"""Profile NASA GSFC outgassing data (ASTM E595): chemistry classes, process levers, noise floor."""
import re, json
import numpy as np, pandas as pd

d = pd.read_csv("data/nasa_outgassing_db.csv")
d = d.dropna(subset=["TML", "CVCM"]).copy()
d = d[(d.TML < 99) & (d.CVCM < 99)]            # 99.99 / 100 are sentinel values
d["name"] = d["Sample Material"].str.upper().str.strip()
d["use"] = d["Material Usage"].fillna("UNKNOWN").str.upper().str.strip()
d["pass"] = (d.TML <= 1.0) & (d.CVCM <= 0.10)
d["tml_ok"] = d.TML <= 1.0
d["cvcm_ok"] = d.CVCM <= 0.10

CHEM = [  # order matters: first match wins
    ("PFPE / fluorinated oil", r"KRYTOX|BRAYCOTE|FOMBLIN|PFPE|PERFLUORO|FLUORINATED (OIL|GREASE|FLUID)|DEMNUM"),
    ("Hydrocarbon lubricant (MAC/PAO)", r"PENNZANE|NYE ?\w*|RHEOLUBE|APIEZON|PAO\b|MULTIPLY.ALKYL"),
    ("Silicone", r"SILICONE|\bRTV\b|SYLGARD|\bDC ?\d|DOW CORNING|CV-?\d{4}|CV1|NUSIL|SCV-|\bSE ?\d{4}|SILASTIC|POLYSILOXANE|\bPDMS\b"),
    ("Epoxy", r"EPOXY|EPON|HYSOL|\bEA ?\d{3,4}|ECCOBOND|STYCAST|ARALDITE|SCOTCH-?WELD|TRA-?BOND|MASTERBOND|MASTER BOND|\bEP\d"),
    ("Cyanate ester", r"CYANATE|RS-3|EX-1515|954-"),
    ("Polyimide", r"POLYIMIDE|KAPTON|VESPEL|\bPI\b|UPILEX|PYRALIN"),
    ("Fluoropolymer (PTFE/FEP/ETFE)", r"PTFE|TEFLON|\bFEP\b|ETFE|TEFZEL|PFA\b|PVDF|KYNAR|HALAR|VITON|FLUOROELASTOMER|KEL-?F|PCTFE"),
    ("Polyurethane", r"URETHANE|POLYURETHANE|\bPU\b|SOLITHANE|CONATHANE|UROCOAT|ARATHANE"),
    ("Acrylic", r"ACRYLIC|ACRYLATE|\bPMMA\b|PLEXIGLAS|LOCTITE"),
    ("Polyester / PET", r"POLYESTER|MYLAR|\bPET\b|DACRON"),
    ("Polyolefin", r"POLYETHYLENE|POLYPROPYLENE|POLYOLEFIN|\bHDPE\b|\bLDPE\b|UHMW"),
    ("PEEK / PEI / PPS / PSU", r"PEEK|ULTEM|\bPEI\b|\bPPS\b|RYTON|POLYSULFONE|\bPSU\b|POLYETHERIMIDE"),
    ("Polyamide / nylon / aramid", r"NYLON|POLYAMIDE|KEVLAR|NOMEX|ARAMID"),
    ("Phenolic", r"PHENOLIC"),
    ("Rubber (EPDM/NBR/neoprene)", r"EPDM|NITRILE|\bNBR\b|NEOPRENE|BUTYL|RUBBER|ELASTOMER"),
    ("Polysulfide", r"POLYSULFIDE|PR-?\d{3,4}|PRC"),
    ("PVC / vinyl", r"\bPVC\b|VINYL"),
]
def chem(n):
    for lab, pat in CHEM:
        if re.search(pat, n): return lab
    return "Unmapped (trade name only)"
d["chem"] = d["name"].map(chem)

d["devol"] = d["name"].str.contains(r"DEVOL|VACUUM BAKE|VAC BAKE|BAKED|BAKE OUT|BAKEOUT|LOW OUTGAS|LOW VOLATIL|\bLV\b|OUTGASSED", regex=True)
d["postcure"] = d["name"].str.contains(r"POST ?CURE|POSTCURE|\bPC\b", regex=True)
temps = d["name"].str.findall(r"(\d{2,3}) ?(?:DEG ?)?C\b")
d["cure_T"] = temps.map(lambda t: max(map(int, t)) if t else np.nan)
def hours(c):
    if not isinstance(c, str): return np.nan
    m = re.match(r"(\d+(?:\.\d+)?)\s*([MHDW])", c.strip().upper())
    if not m: return np.nan
    v, u = float(m.group(1)), m.group(2)
    return v * {"M": 1/60, "H": 1, "D": 24, "W": 168}[u]
d["cure_h"] = d["Cure"].map(hours)
d["idnum"] = pd.to_numeric(d["ID"].str.extract(r"(\d+)")[0], errors="coerce")
# product family: strip lot/batch/run tokens to group replicate tests
def family(n):
    n = re.sub(r"\b(LOT|BATCH|RUN|S/N|SN|SER|SAMPLE|PART|P/N|DATE)\b.*$", "", n)
    n = re.sub(r"\(.*?\)", "", n)
    return re.sub(r"\s+", " ", n).strip()
d["family"] = d["name"].map(family)
d["brand"] = d["family"].str.split().str[0]

out = {}
out["n"] = int(len(d))
out["pass_rate"] = round(float(d["pass"].mean()), 3)
out["fail_mode"] = {
    "both_fail": round(float((~d.tml_ok & ~d.cvcm_ok).mean()), 3),
    "tml_only": round(float((~d.tml_ok & d.cvcm_ok).mean()), 3),
    "cvcm_only": round(float((d.tml_ok & ~d.cvcm_ok).mean()), 3),
}
out["mapped_share"] = round(float((d.chem != "Unmapped (trade name only)").mean()), 3)
g = d.groupby("chem").agg(n=("pass", "size"), pass_rate=("pass", "mean"),
                         tml_med=("TML", "median"), cvcm_med=("CVCM", "median")).sort_values("n", ascending=False)
out["by_chem"] = g.round(3).reset_index().to_dict("records")
u = d.groupby("use").agg(n=("pass", "size"), pass_rate=("pass", "mean")).query("n>=80").sort_values("pass_rate")
out["by_use"] = u.round(3).reset_index().to_dict("records")

# near misses: fail, but within 2x of both limits -> candidates for process fixes (devolatilization, post-cure)
fail = d[~d["pass"]]
near = fail[(fail.TML <= 2.0) & (fail.CVCM <= 0.20)]
out["near_miss_share_of_fails"] = round(len(near) / len(fail), 3)

# process lever: devolatilized / baked vs not, within the same chemistry class
lev = []
for c, sub in d.groupby("chem"):
    a, b = sub[sub.devol], sub[~sub.devol]
    if len(a) >= 8 and len(b) >= 30:
        lev.append({"chem": c, "n_devol": len(a), "pass_devol": round(a["pass"].mean(), 3), "pass_other": round(b["pass"].mean(), 3),
                    "cvcm_med_devol": float(a.CVCM.median()), "cvcm_med_other": float(b.CVCM.median())})
out["devol_lever"] = lev
out["devol_overall"] = {"n": int(d.devol.sum()), "pass_devol": round(float(d[d.devol]["pass"].mean()), 3), "pass_other": round(float(d[~d.devol]["pass"].mean()), 3)}

# paired: same family tested as-is and devolatilized
pairs = []
for f, sub in d.groupby(d["family"].str.replace(r"\b(DEVOL\w*|VACUUM BAKED|VAC BAKED|BAKED|LOW OUTGAS\w*)\b", "", regex=True).str.replace(r"\s+", " ", regex=True).str.strip()):
    a, b = sub[sub.devol], sub[~sub.devol]
    if len(a) and len(b):
        pairs.append({"family": f, "tml_as_is": float(b.TML.median()), "tml_devol": float(a.TML.median()),
                      "cvcm_as_is": float(b.CVCM.median()), "cvcm_devol": float(a.CVCM.median()),
                      "pass_as_is": float(b["pass"].mean()), "pass_devol": float(a["pass"].mean())})
pp = pd.DataFrame(pairs)
out["devol_pairs"] = {"n_pairs": len(pp), "median_cvcm_ratio": round(float((pp.cvcm_devol + 0.005).div(pp.cvcm_as_is + 0.005).median()), 3) if len(pp) else None,
                      "median_tml_ratio": round(float((pp.tml_devol + 0.005).div(pp.tml_as_is + 0.005).median()), 3) if len(pp) else None,
                      "pass_as_is": round(float(pp.pass_as_is.mean()), 3) if len(pp) else None, "pass_devol": round(float(pp.pass_devol.mean()), 3) if len(pp) else None}

# noise floor: families tested >=3 times; how often does a test agree with the family's majority verdict?
rep = d.groupby("family").filter(lambda s: len(s) >= 3)
maj = rep.groupby("family")["pass"].transform(lambda s: s.mean() >= 0.5)
out["noise_floor"] = {"families": int(rep.family.nunique()), "tests": int(len(rep)),
                      "agree_with_family_majority": round(float((rep["pass"] == maj).mean()), 3),
                      "families_with_mixed_verdicts": round(float(rep.groupby("family")["pass"].agg(lambda s: 0 < s.mean() < 1).mean()), 3),
                      "median_within_family_cvcm_cv": round(float(rep.groupby("family").CVCM.agg(lambda s: s.std() / (s.mean() + 1e-3)).median()), 3)}

# time trend using GSC sample ID as a proxy for test order
dd = d.dropna(subset=["idnum"]).copy()
dd["era"] = pd.qcut(dd.idnum, 5, labels=False)
out["trend"] = dd.groupby("era").agg(id_lo=("idnum", "min"), id_hi=("idnum", "max"), n=("pass", "size"), pass_rate=("pass", "mean")).round(3).reset_index().to_dict("records")

# oils: viscosity/grade vs volatility check (physics sanity)
out["oil_examples"] = d[d.name.str.contains(r"DC 200|KRYTOX 143|FOMBLIN Z|BRAYCOTE 815|PENNZANE|NYE ?(?:SYNTHETIC )?OIL 2001", regex=True)][["name", "TML", "CVCM"]].head(40).to_dict("records")

json.dump(out, open("profile.json", "w"), indent=1, default=float)
d.to_pickle("og_clean.pkl")
for k, v in out.items():
    if k in ("by_chem", "by_use", "trend", "devol_lever"):
        print(k); print(pd.DataFrame(v).to_string(index=False))
    elif k != "oil_examples":
        print(k, v)
