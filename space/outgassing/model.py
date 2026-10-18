"""First outgassing screen: predict ASTM E595 pass/fail from what a buyer knows before testing.
Evaluated on unseen product families and unseen brands, against the replicate-test noise floor."""
import json, warnings
import numpy as np, pandas as pd
from scipy import sparse
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.preprocessing import OneHotEncoder
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GroupKFold
from sklearn.metrics import roc_auc_score, accuracy_score, brier_score_loss
warnings.filterwarnings("ignore")

d = pd.read_pickle("og_clean.pkl").reset_index(drop=True)
y = d["pass"].astype(int).values
d["cure_h_b"] = pd.cut(d.cure_h, [-1, 2, 8, 30, 100, 1e6], labels=["<=2h", "2-8h", "8-30h", "30-100h", ">100h"]).astype(str)
d["cure_T_b"] = pd.cut(d.cure_T, [-1, 40, 80, 130, 1e3], labels=["<=40C", "40-80C", "80-130C", ">130C"]).astype(str)
cat_cols_base = ["chem", "use"]
cat_cols_full = ["chem", "use", "MFR", "cure_h_b", "cure_T_b", "devol", "postcure"]

def feats(tr, te, cols, text):
    enc = OneHotEncoder(handle_unknown="ignore", min_frequency=5)
    Xa = enc.fit_transform(d.loc[tr, cols].astype(str)); Xb = enc.transform(d.loc[te, cols].astype(str))
    if text:
        tf = TfidfVectorizer(analyzer="char_wb", ngram_range=(3, 5), min_df=3, sublinear_tf=True, max_features=60000)
        Xa = sparse.hstack([Xa, tf.fit_transform(d.loc[tr, "name"])]).tocsr(); Xb = sparse.hstack([Xb, tf.transform(d.loc[te, "name"])]).tocsr()
    return Xa, Xb

def cv(groups, cols, text, C=1.0):
    p = np.zeros(len(d))
    for tr, te in GroupKFold(5).split(d, y, groups):
        Xa, Xb = feats(tr, te, cols, text)
        m = LogisticRegression(C=C, max_iter=3000).fit(Xa, y[tr])
        p[te] = m.predict_proba(Xb)[:, 1]
    return p

res = {}
for split, groups in [("new_family", d.family.values), ("new_brand", d.brand.values)]:
    for name, cols, text in [("prior: class + usage", cat_cols_base, False),
                             ("+ maker, cure, process flags", cat_cols_full, False),
                             ("+ name text (char n-grams)", cat_cols_full, True)]:
        p = cv(groups, cols, text)
        res[f"{split} | {name}"] = {"auc": round(roc_auc_score(y, p), 3), "acc": round(accuracy_score(y, p > 0.5), 3),
                                    "brier": round(brier_score_loss(y, p), 3)}
        if split == "new_family" and text: p_best = p
for k, v in res.items(): print(k, v)

# triage value: if a buyer only sends the top-ranked candidates to the lab, how many tests per passing material?
d["p"] = p_best
dd = d.sort_values("p", ascending=False)
tri = []
for q in [1.0, 0.5, 0.3, 0.2, 0.1]:
    top = dd.head(int(len(dd) * q))
    tri.append({"test_top_share": q, "pass_rate": round(float(top["pass"].mean()), 3),
                "tests_per_pass": round(1 / top["pass"].mean(), 2),
                "share_of_all_passes_kept": round(float(top["pass"].sum() / d["pass"].sum()), 3)})
print(pd.DataFrame(tri).to_string(index=False))
# calibration by decile
d["dec"] = pd.qcut(d.p, 10, labels=False, duplicates="drop")
cal = d.groupby("dec").agg(pred=("p", "mean"), obs=("pass", "mean"), n=("p", "size")).round(3)
print(cal.to_string())

# alternative criterion many programs use: water-corrected TML (TML - WVR <= 1.0)
has = d.WVR.notna()
alt = ((d.TML - d.WVR.fillna(0)) <= 1.0) & (d.CVCM <= 0.10)
res_alt = {"pass_tml_minus_wvr": round(float(alt[has].mean()), 3), "pass_standard_same_rows": round(float(d.loc[has, "pass"].mean()), 3),
           "nylon_standard": round(float(d.loc[has & (d.chem == "Polyamide / nylon / aramid"), "pass"].mean()), 3),
           "nylon_water_corrected": round(float(alt[has & (d.chem == "Polyamide / nylon / aramid")].mean()), 3)}
print(res_alt)
json.dump({"cv": res, "triage": tri, "calibration": cal.reset_index().to_dict("records"), "water_corrected": res_alt},
          open("model.json", "w"), indent=1, default=float)
