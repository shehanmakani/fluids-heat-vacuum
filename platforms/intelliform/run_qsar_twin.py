import json, os, sys, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, "/home/claude/cheme-nova/intelliform")
import pandas as pd
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from modules import qsar
from modules.digital_twin import simulate_manufacturing

db = pd.read_csv("/home/claude/cheme-nova/intelliform/data/ingredients_db.csv")
card = qsar.initialize_models(db)
print("QSAR tier:", card.active_tier, "n_train", card.n_training)

# Representative structures by chemistry class (the commercial fluids are mixtures; structures are proxies)
REP = {
 "PAO2":  ("hydrogenated 1-decene dimer (C20 branched alkane)", "CCCCCCCCCCC(C)CCCCCCCC", 0),
 "PAO4":  ("hydrogenated 1-decene trimer (C30 branched alkane)", "CCCCCCCCC(CCCCCCCC)CC(C)CCCCCCCCCC", 0),
 "PAO6":  ("hydrogenated 1-decene tetramer (C40 branched alkane)", "CCCCCCCCC(CCCCCCCC)CC(CCCCCCCC)CC(C)CCCCCCCCCC", 0),
 "EF7010":("diester proxy: bis(2-ethylhexyl) sebacate", "CCCCC(CC)COC(=O)CCCCCCCCC(=O)OCC(CC)CCCC", 1),
 "EF3221":("monoester proxy: 2-ethylhexyl oleate", "CCCCCCCC/C=C\\CCCCCCCC(=O)OCC(CC)CCCC", 1),
 "P1976": ("polyol ester proxy: neopentyl glycol dioleate", "CCCCCCCC/C=C\\CCCCCCCC(=O)OCC(C)(C)COC(=O)CCCCCCC/C=C\\CCCCCCCC", 1),
 "NC1000":("natural ester proxy: triolein", "CCCCCCCC/C=C\\CCCCCCCC(=O)OCC(COC(=O)CCCCCCC/C=C\\CCCCCCCC)OC(=O)CCCCCCC/C=C\\CCCCCCCC", 1),
 "HYV2":  ("naphthenic mineral oil proxy: hexyldecalin", "CCCCCCC1CCC2CCCCC2C1", 0),
 "RHC8":  ("renewable isoparaffin proxy: phytane (C20)", "CC(C)CCCC(C)CCCC(C)CCCC(C)CC", 1),
}
rows = []
for k, (desc, smi, rb) in REP.items():
    m = Chem.MolFromSmiles(smi)
    p = qsar.predict_properties(smi)
    rows.append({"stock": k, "proxy": desc, "formula": CalcMolFormula(m), "datasheet_readily_biodeg": rb,
                 "qsar_biodeg": round(p.biodegradability, 1), "qsar_ecotox": round(p.ecotoxicity, 1),
                 "confidence": p.confidence, "used_ml": p.used_ml})
q = pd.DataFrame(rows); print(q.to_string(index=False))
# agreement: QSAR biodeg >= 70 counts as "readily"
q["qsar_readily"] = (q.qsar_biodeg >= 70).astype(int)
agree = (q.qsar_readily == q.datasheet_readily_biodeg).mean()
print("agreement with data sheets:", agree)

# Digital twin: 1,000 kg production batch of each recommended blend at ChemRich NJ (simple blending)
names = {"PAO2": "Synfluid PAO 2", "PAO4": "Synfluid PAO 4", "PAO6": "Synfluid PAO 6", "EF7010": "Priolube EF 7010",
         "P1976": "Priolube 1976", "RHC8": "Renewable isoparaffin", "EF3221": "Priolube EF 3221"}
twdb = pd.DataFrame([{"Ingredient": names[k], "Function": "Dielectric Base Fluid", "Vertical": "industrial"} for k in names])
pareto = json.load(open("intelliform_pareto.json"))
picks = {"Lead A (round 3)": {"PAO4": 40, "P1976": 35, "PAO2": 25}, "Lead G (round 3)": {"RHC8": 62, "P1976": 38},
         "IntelliForm TOPSIS, standard": pareto["standard_150C"]["recommended"]["blend"],
         "IntelliForm TOPSIS, premium": pareto["premium_210C"]["recommended"]["blend"]}
tw = {}
for nm, b in picks.items():
    blend = {names[k]: float(v) for k, v in b.items()}
    r = simulate_manufacturing(blend, twdb, batch_kg=1000.0, manufacturing_process="blending", vertical="industrial")
    d = r.__dict__ if hasattr(r, "__dict__") else r
    keep = {k: v for k, v in d.items() if isinstance(v, (int, float, str)) }
    tw[nm] = {"blend": b, **keep, "risks": [getattr(x, "__dict__", x) for x in (d.get("scale_up_risks") or d.get("risks") or [])][:4]}
    print(nm, {k: keep[k] for k in list(keep)[:12]})
json.dump({"qsar": rows, "agreement": agree, "qsar_tier": card.active_tier, "twin": tw}, open("intelliform_qsar_twin.json", "w"), indent=1, default=str)
