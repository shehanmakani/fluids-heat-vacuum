"""Physics simulator of the ASTM E595 outgassing test (125 C, 24 h, vacuum, collector at 25 C).

The E595 specimen compartment is an effusion (Knudsen) cell: vapor leaves through one small exit
port toward a cooled collector plate. For each volatile species we integrate
  - effusion out of the compartment: partial pressure (Raoult on remaining mole fractions) x Hertz-Knudsen
    flux x effective port area,
  - capture on the collector (fraction ETA) and re-evaporation from the 25 C deposit (Raoult within deposit).
TML = mass lost / initial mass; CVCM = deposit mass at 24 h / initial mass.

Inputs per species: mass fraction, molar mass, and two vapor-pressure points (or one point + dH).
No parameter is fitted to TML. ETA (collector capture) is set from the two diffusion-pump fluids whose
condensable fraction is unambiguous (DC 705, polyphenyl ether).
"""
import json, math
import numpy as np


R = 8.314
T_SAMPLE, T_COLL, HOURS = 398.15, 298.15, 24.0
M0 = 300e-6                         # kg, nominal specimen mass (about 300 mg)
D_PORT = 6.3e-3                     # m, exit-port diameter (as commonly cited for E595; verify against the standard)
A_SURF = 1.5e-4                     # m^2, exposed specimen surface (boat), sensitivity tested
A_PORT = math.pi * D_PORT ** 2 / 4
A_EFF = A_PORT * A_SURF / (A_PORT + A_SURF)   # series resistance: surface evaporation then port effusion
A_DEP = 5e-4                        # m^2, deposit spot on collector
ETA = 0.85                          # fraction of effused vapor that reaches the collector

def vp_fun(sp):
    """Clausius-Clapeyron through (T1,P1),(T2,P2) in C and Pa, or (T1,P1)+dH."""
    (t1, p1) = sp["vp"][0]
    if len(sp["vp"]) > 1:
        t2, p2 = sp["vp"][1]
        b = math.log(p2 / p1) / (1 / (t1 + 273.15) - 1 / (t2 + 273.15))
    else:
        b = sp["dH"] / R
    return lambda T: p1 * math.exp(b * (1 / (t1 + 273.15) - 1 / T)), b * R

def hk(P, M, T):            # Hertz-Knudsen mass flux, kg m^-2 s^-1
    return P * math.sqrt(M / (2 * math.pi * R * T))

def simulate(species, m0=M0, a_eff=A_EFF, eta=ETA, a_dep=A_DEP, hours=HOURS, t_sample=T_SAMPLE, t_coll=T_COLL):
    n = len(species)
    vps = [vp_fun(s)[0] for s in species]
    Ms = np.array([s["M"] for s in species])
    w = np.array([s["w"] for s in species]) * m0
    vol = np.array([s.get("volatile", True) for s in species])
    vp_s = np.array([vps[i](t_sample) for i in range(n)]); vp_c = np.array([vps[i](t_coll) for i in range(n)])
    k_s = vp_s * np.sqrt(Ms / (2 * math.pi * R * t_sample)) * a_eff * vol   # kg/s at mole fraction 1
    k_c = vp_c * np.sqrt(Ms / (2 * math.pi * R * t_coll)) * a_dep
    liq, dep = w.copy(), np.zeros(n)
    dt = 20.0                                       # s, explicit step; fluxes vary slowly relative to this
    for _ in range(int(hours * 3600 / dt)):
        ml = liq / Ms; xl = ml / max(ml.sum(), 1e-30)
        md = dep / Ms; xd = md / max(md.sum(), 1e-30)
        out = np.minimum(k_s * xl * dt, liq)
        rev = np.minimum(k_c * xd * dt, dep)
        liq = liq - out; dep = dep + eta * out - rev
    return {"TML": 100 * (m0 - liq.sum()) / m0, "CVCM": 100 * dep.sum() / m0,
            "per_species": [{"name": s["name"], "lost_pct": 100 * (w[i] - liq[i]) / m0, "deposit_pct": 100 * dep[i] / m0} for i, s in enumerate(species)]}

# ---- validation fluids: single, well-characterized molecules or narrow fluids with published vapor pressure ----
TORR = 133.322
FLUIDS = {
    "DC 704 (tetraphenyl tetramethyl trisiloxane)": dict(M=0.4848, vp=[(25, 2e-8 * TORR), (210, 0.5 * TORR)],
        nasa=[(99.85, 11.30)], vp_src="Duniway: 2e-8 Torr at 25 C; bp 210 C at 0.5 Torr"),
    "DC 705 (pentaphenyl trimethyl trisiloxane)": dict(M=0.5469, vp=[(25, 2.4e-9 * TORR), (243, 0.5 * TORR)],
        nasa=[(97.20, 81.00)], vp_src="Duniway: 2.4e-9 Torr at 25 C; bp 243 C at 0.5 Torr"),
    "Polyphenyl ether (5P4E class; NASA entry is Convalex-10)": dict(M=0.4465, vp=[(25, 4e-10 * TORR), (275, 0.5 * TORR)],
        nasa=[(46.40, 40.00)], vp_src="Santovac 5 (5P4E), Duniway: 4e-10 Torr at 25 C; bp 275 C at 0.5 Torr. Convalex-10 is a related polyphenyl ether, not identical"),
    "MAC (Pennzane / Nye 2001 family)": dict(M=0.91, vp=[(25, 3.1e-10), (150, 2.3e-3)],
        nasa=[(0.47, 0.23), (1.42, 0.61), (1.33, 0.55), (0.47, 0.08), (0.61, 0.07), (0.93, 0.14), (0.21, 0.05)],
        vp_src="NASA CR-2009-215637 (Pennzane X2000)"),
    "PFPE (Fomblin Z25 / Brayco 815Z)": dict(M=9.5, vp=[(25, 1.2e-11), (150, 1.1e-5)],
        nasa=[(0.03, 0.00), (0.06, 0.02)], vp_src="NASA CR-2009-215637 (Fomblin Z25); molar mass assumed 9.5 kg/mol"),
}

def run_validation(**kw):
    rows = []
    for name, f in FLUIDS.items():
        r = simulate([{"name": name, "w": 1.0, "M": f["M"], "vp": f["vp"]}], **kw)
        tm = np.median([a for a, b in f["nasa"]]); cm = np.median([b for a, b in f["nasa"]])
        vp, dH = vp_fun(f)
        rows.append({"fluid": name, "P125_Pa": vp(T_SAMPLE), "dH_kJ": dH / 1000, "M": f["M"],
                     "TML_pred": r["TML"], "CVCM_pred": r["CVCM"], "TML_nasa_median": float(tm), "CVCM_nasa_median": float(cm),
                     "TML_nasa_range": [min(a for a, b in f["nasa"]), max(a for a, b in f["nasa"])],
                     "CVCM_nasa_range": [min(b for a, b in f["nasa"]), max(b for a, b in f["nasa"])], "n_tests": len(f["nasa"]), "vp_src": f["vp_src"]})
    return rows

def pass_threshold(limit_pct, condensable=False, eta=ETA, a_eff=A_EFF):
    """Largest P125*sqrt(M) (Pa*(kg/mol)^0.5) a pure low-volatility fluid can have and still meet the limit,
    in the small-loss regime where loss is linear in time."""
    allowed = limit_pct / 100 * M0 / (eta if condensable else 1.0)
    J = allowed / (a_eff * HOURS * 3600)
    return J * math.sqrt(2 * math.pi * R * T_SAMPLE)

if __name__ == "__main__":
    out = {"geometry": {"A_port_cm2": A_PORT * 1e4, "A_surface_cm2": A_SURF * 1e4, "A_eff_cm2": A_EFF * 1e4, "eta": ETA, "A_dep_cm2": A_DEP * 1e4}}
    val = run_validation(); out["validation"] = val
    for r in val:
        print(f"{r['fluid'][:48]:48s} P125 {r['P125_Pa']:.1e} Pa  TML {r['TML_pred']:6.2f} (NASA {r['TML_nasa_median']:6.2f})  CVCM {r['CVCM_pred']:6.2f} (NASA {r['CVCM_nasa_median']:6.2f})")
    # sensitivity: port diameter 4-8 mm, surface 0.5-3 cm2, eta 0.6-1.0
    sens = []
    for d in (4e-3, 6.3e-3, 8e-3):
        for a in (0.5e-4, 1.5e-4, 3e-4):
            ap = math.pi * d * d / 4; ae = ap * a / (ap + a)
            v = run_validation(a_eff=ae)
            sens.append({"port_mm": d * 1e3, "surf_cm2": a * 1e4, "A_eff_cm2": ae * 1e4, "TML": {r["fluid"][:12]: round(r["TML_pred"], 3) for r in v}})
    out["sensitivity"] = sens
    # log error on TML for the low-volatility fluids where the test is not saturated
    out["tml_log10_error"] = {r["fluid"][:30]: round(math.log10(r["TML_pred"] / r["TML_nasa_median"]), 2) for r in val}
    # design rule
    out["design_rule"] = {"TML_1pct_P125_sqrtM": pass_threshold(1.0), "CVCM_0p1pct_P125_sqrtM": pass_threshold(0.1, condensable=True),
                          "units": "Pa*(kg/mol)^0.5 at 125 C"}
    # predictions for the lubricant candidates from the space-lubricant model
    sihc = {"name": "penta-silahydrocarbon", "w": 1, "M": 0.91, "vp": [(25, 3.0e-12)], "dH": vp_fun(FLUIDS["MAC (Pennzane / Nye 2001 family)"])[1]}
    target = {"name": "OrbitChem target oil", "w": 1, "M": 0.91, "vp": [(25, 1.0e-11)], "dH": sihc["dH"]}
    out["lubricant_predictions"] = {s["name"]: {k: round(v, 4) for k, v in simulate([s]).items() if k != "per_species"} for s in (sihc, target)}
    # diagnosis demo: a passing base oil carrying a volatile tail, and what devolatilization does
    base = {"name": "base (MAC-like)", "M": 0.91, "vp": [(25, 3.1e-11), (150, 2.3e-4)]}   # a 10x less volatile, distilled MAC grade (illustrative)
    tail = {"name": "volatile tail", "M": 0.63, "vp": [(25, 3.1e-9), (150, 2.3e-2)]}      # 100x more volatile homolog (illustrative)
    demo = []
    for f in (0.0, 0.005, 0.01, 0.02, 0.04):
        r = simulate([{**base, "w": 1 - f}, {**tail, "w": f}])
        demo.append({"tail_pct": f * 100, "TML": round(r["TML"], 3), "CVCM": round(r["CVCM"], 3)})
    out["tail_demo"] = demo
    # devolatilization: 300 mg/cm2 film, 150 C, open evaporation (no port), t hours, then E595
    devol = []
    for hrs in (0, 2, 8, 24):
        if hrs == 0:
            comp = [{**base, "w": 0.98}, {**tail, "w": 0.02}]
        else:
            r = simulate([{**base, "w": 0.98}, {**tail, "w": 0.02}], a_eff=1e-4, hours=hrs, t_sample=423.15, eta=0.0)
            lost = {p["name"]: p["lost_pct"] / 100 for p in r["per_species"]}
            wb, wt = 0.98 - lost["base (MAC-like)"], 0.02 - lost["volatile tail"]
            comp = [{**base, "w": wb / (wb + wt)}, {**tail, "w": wt / (wb + wt)}]
        e = simulate(comp)
        devol.append({"devol_hours_150C": hrs, "tail_left_pct": round(comp[1]["w"] * 100, 3), "TML": round(e["TML"], 3), "CVCM": round(e["CVCM"], 3),
                      "pass": e["TML"] <= 1 and e["CVCM"] <= 0.1})
    out["devol_demo"] = devol
    print(json.dumps({k: out[k] for k in ("tml_log10_error", "design_rule", "lubricant_predictions", "tail_demo", "devol_demo")}, indent=1, default=float))
    json.dump(out, open("e595sim.json", "w"), indent=1, default=float)
