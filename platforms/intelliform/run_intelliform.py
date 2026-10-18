"""Run ThermalForm through IntelliForm v1.6: NSGA-III + TOPSIS, QSAR, digital twin."""
import json, os, sys, time
os.environ["THERMALFORM_PATH"] = "/home/claude/fluids-heat-vacuum/thermalform"
sys.path.insert(0, "/home/claude/cheme-nova/intelliform")
from modules.thermalform_vertical import run_thermalform_pareto
out = {}
for name, fmin in [("standard_150C", 150.0), ("premium_210C", 210.0)]:
    t = time.time()
    res, details, ctx = run_thermalform_pareto(flash_min=fmin)
    rec = res.recommended
    print(f"\n== {name}: {res.n_solutions} Pareto blends, {time.time()-t:.0f}s, backend {res.backend}")
    for s in sorted(res.frontier, key=lambda s: -s.perf_score)[:12]:
        print(f"  {s.blend}  pump ${s.cost_per_kg:,.0f}/MW-yr  bio {s.bio_pct}%  perf {s.perf_score}")
    d = next(x for x in details if x["blend"] == rec.blend)
    print("  TOPSIS recommended:", rec.blend, {k: d[k] for k in ("dT_vs_bench", "pump_usd_vs_bench", "dk_worst", "flash_c", "bio_pct")})
    out[name] = {"n_front": res.n_solutions, "recommended": d,
                 "front": [next(x for x in details if x["blend"] == s.blend) for s in res.frontier]}
json.dump(out, open("/home/claude/integration/intelliform_pareto.json", "w"), indent=1, default=float)
