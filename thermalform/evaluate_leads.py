"""Reproduce the round 3 lead-blend scorecard from the published recipes.

The recipes came out of ChemeNova's IntelliForm optimizer, which is not part of this repo.
This script evaluates them with the open property, heat-sink, loop, dielectric and OCP models.
Run: python evaluate_leads.py
"""
import json
import numpy as np
from thermalform.properties import load_base_stocks, evaluate, fluid_at
from thermalform.physics import HeatSink, forced
from thermalform.loop import loop
from thermalform.ocp import check, fom1
from thermalform.dielectric import fluid_effect

bs = load_base_stocks()
IDX = {k: i for i, k in enumerate(bs.ids)}
Q, T_IN, T_BASE_MAX = 1000.0, 40.0, 70.0      # 1 kW chip, 40 C supply, 70 C heat-sink base limit
FLOW, SINK = 6 / 60000, HeatSink()            # 6 L/min per heat sink, 1.5 mm fin gap

def mix(**fr):
    w = np.zeros(bs.n)
    for k, v in fr.items():
        w[IDX[k]] = v
    return w / w.sum()

FLUIDS = {
    "PAO benchmark": mix(PAO4=0.88, PAO2=0.12),
    "Lead A (all-round)": mix(PAO4=0.40, P1976=0.35, PAO2=0.25),
    "Lead G (biodegradable)": mix(RHC8=0.62, P1976=0.38),
    "Lead P (premium flash)": mix(PAO4=0.81, EF7010=0.19),
    "All-ester (needs GHz Dk data)": mix(EF7010=0.52, P1976=0.48),
}

bench = FLUIDS["PAO benchmark"]
rb, lb = forced(bench, bs, FLOW, SINK, Q, T_IN), loop(bench, bs)
rows = []
for name, w in FLUIDS.items():
    r, lp, oc, p = forced(w, bs, FLOW, SINK, Q, T_IN), loop(w, bs), check(w, bs), evaluate(w, bs)
    rows.append({
        "fluid": name,
        "chip_dT_vs_bench_C": round((r["R_fluid"] - rb["R_fluid"]) * Q, 3),
        "max_chip_W_at_70C": round((T_BASE_MAX - T_IN) / r["R_fluid"]),
        "pump_kW_per_MW": round(lp["pump_kw"], 2),
        "pump_usd_yr_vs_bench": round(lp["pump_usd_yr"] - lb["pump_usd_yr"]),
        "dk_central": round(oc["dk_ghz"], 3), "dk_worst": round(oc["dk_ghz_worst"], 3),
        "flash_C_model": round(p["flash_c"]), "biodegradable_frac": round(p["biodegradable_frac"], 2),
        "trace_ohm_worst": round(fluid_effect(oc["dk_ghz_worst"])["z_fluid"], 2),
    })
print(f"{'fluid':32s} {'dT C':>7s} {'Qmax W':>7s} {'pump kW':>8s} {'$/yr vs PAO':>12s} {'Dk worst':>9s} {'flash':>6s}")
for x in rows:
    print(f"{x['fluid']:32s} {x['chip_dT_vs_bench_C']:7.2f} {x['max_chip_W_at_70C']:7d} {x['pump_kW_per_MW']:8.2f} {x['pump_usd_yr_vs_bench']:12,d} {x['dk_worst']:9.3f} {x['flash_C_model']:6d}")
print(f"\nCheck: PAO 4 OCP FOM1 at 50 C = {fom1(fluid_at(mix(PAO4=1), bs, 50.0)):.1f} (Evonik publishes 46)")
json.dump(rows, open("results/leads.json", "w"), indent=1)
