"""Operating-window model for vacuum space lubricant base oils.

Cold limit: temperature where kinematic viscosity reaches NU_COLD (startup torque proxy),
            from a Walther (ASTM D341) fit through two published viscosities.
Hot limit:  temperature where Langmuir free-molecular evaporation from 1 cm^2 of free
            surface removes LOSS_MG over LIFE_YR (unsealed, worst case), using a
            Clausius-Clapeyron fit through published vapor pressures.
Blends: Refutas viscosity mixing (mass basis); Raoult's law on mole fractions for vapor
pressure (ideal solution assumption).
All inputs are literature values; see DATA notes.
"""
import math, json
import numpy as np

R = 8.314
NU_COLD = 10000.0      # cSt
LOSS_MG, AREA_CM2, LIFE_YR = 5.0, 1.0, 15.0

def walther(t1, n1, t2, n2):
    z = lambda n: math.log10(math.log10(n + 0.7))
    a1, a2 = math.log10(t1 + 273.15), math.log10(t2 + 273.15)
    B = (z(n1) - z(n2)) / (a2 - a1); A = z(n1) + B * a1
    return lambda T: 10 ** (10 ** (A - B * math.log10(T + 273.15))) - 0.7

def clausius(t1, p1, t2, p2):
    b = (math.log(p1) - math.log(p2)) / (1 / (t2 + 273.15) - 1 / (t1 + 273.15))
    a = math.log(p1) + b / (t1 + 273.15)
    return lambda T: math.exp(a - b / (T + 273.15)), b

DATA = {
 # PFPE: viscosity Brayco 815Z data sheet (-40, 38, 99 C); vapor pressure Fomblin Z25 (NASA CR-2009-215637)
 "PFPE (Brayco 815Z / Fomblin Z25)": dict(visc=((-40, 6500), (38, 148), (99, 45)), vp=((25, 1.2e-11), (150, 1.1e-5)),
                                          M=9.5, pour=-72, pfas=True, rho=1.85),
 # MAC: Pennzane 2001A viscosity 38/99 C (NASA TM-2001-211196); vapor pressure X2000 = 2001A base (NASA CR-2009-215637)
 "MAC (Pennzane 2001A)": dict(visc=((-20, 8500), (38, 108), (99, 15)), vp=((25, 3.1e-10), (150, 2.3e-3)), M=0.91, pour=-55, pfas=False, rho=0.84),
 # Penta-silahydrocarbon MJD991029: viscosity 38/99 C and VP at 25 C (NASA TM-2001-211196); VP slope and MW ASSUMED = MAC X2000
 "Silahydrocarbon (penta, MJD991029)": dict(visc=((-20, 3500), (38, 77), (99, 13)), vp=((25, 3.0e-12), None), M=0.91, pour=None, pfas=False, rho=0.84,
                                            assumed="VP temperature slope and molar mass assumed equal to MAC X2000"),
 # -20 C viscosities for MAC and silahydrocarbon are NASA's extrapolated values (TM-2001-211196)
}

def build(name, d, ref_b=None):
    v = sorted(d["visc"])
    # cold-limit accuracy matters most: fit through the two coldest published points
    nu = walther(v[0][0], v[0][1], v[1][0], v[1][1])
    if d["vp"][1] is not None:
        pv, b = clausius(d["vp"][0][0], d["vp"][0][1], d["vp"][1][0], d["vp"][1][1])
    else:
        t1, p1 = d["vp"][0]; b = ref_b
        a = math.log(p1) + b / (t1 + 273.15)
        pv = lambda T, a=a, b=b: math.exp(a - b / (T + 273.15))
    return nu, pv, b

def evap_mg(pv, M, T):
    Tk = T + 273.15
    w = pv(T) * math.sqrt(M / (2 * math.pi * R * Tk))        # kg/m2/s
    return w * AREA_CM2 * 1e-4 * LIFE_YR * 3.156e7 * 1e6     # mg over life

def solve(f, lo, hi):
    for _ in range(200):
        m = (lo + hi) / 2
        if f(m) > 0: hi = m
        else: lo = m
    return (lo + hi) / 2

def window(nu, pv, M):
    tc = solve(lambda T: NU_COLD - nu(T), -120, 100)   # nu falls with T
    th = solve(lambda T: evap_mg(pv, M, T) - LOSS_MG, -50, 300)
    return tc, th

# ---- blends: MAC + silahydrocarbon (both PFAS-free) ----
def vbn(n): return 14.534 * math.log(math.log(n + 0.8)) + 10.975
def vbn_inv(x): return math.exp(math.exp((x - 10.975) / 14.534)) - 0.8

out = {"assumptions": dict(nu_cold_cst=NU_COLD, loss_mg=LOSS_MG, area_cm2=AREA_CM2, life_yr=LIFE_YR), "fluids": {}, "blends": []}
models = {}
b_mac = None
for n, d in DATA.items():
    nu, pv, b = build(n, d, ref_b=b_mac)
    if n.startswith("MAC"): b_mac = b
    models[n] = (nu, pv, d["M"])
for n, d in DATA.items():   # rebuild silahydrocarbon now that b_mac exists
    nu, pv, b = build(n, d, ref_b=b_mac); models[n] = (nu, pv, d["M"])
    tc, th = window(nu, pv, d["M"])
    out["fluids"][n] = dict(T_cold=tc, T_hot=th, width=th - tc, pour=d["pour"], pfas=d["pfas"],
                            nu_m40=nu(-40), nu_m20=nu(-20), nu_20=nu(20), nu_100=nu(100),
                            vp25=pv(25), vp100=pv(100), dH_kJ=b * R / 1000, note=d.get("assumed", ""))

numac, pvmac, Mmac = models["MAC (Pennzane 2001A)"]
nusi, pvsi, Msi = models["Silahydrocarbon (penta, MJD991029)"]
for x in np.linspace(0, 1, 11):           # mass fraction silahydrocarbon
    nu = lambda T, x=x: vbn_inv((1 - x) * vbn(numac(T)) + x * vbn(nusi(T)))
    nm, ns = (1 - x) / Mmac, x / Msi
    xm, xs = nm / (nm + ns), ns / (nm + ns)
    pv = lambda T, xm=xm, xs=xs: xm * pvmac(T) + xs * pvsi(T)
    tc, th = window(nu, pv, Mmac)
    out["blends"].append(dict(x_sihc=float(x), T_cold=tc, T_hot=th))

json.dump(out, open("spacelub_results.json", "w"), indent=1, default=float)
for n, r in out["fluids"].items():
    print(f"{n:38s} cold {r['T_cold']:6.1f} C  hot {r['T_hot']:6.1f} C  width {r['width']:5.1f}  pour {r['pour']}  nu(-40) {r['nu_m40']:8.0f}  vp25 {r['vp25']:.1e} dH {r['dH_kJ']:.0f} kJ/mol")
for b in out["blends"]:
    print(f"MAC + {b['x_sihc']*100:3.0f}% SiHC: cold {b['T_cold']:6.1f}  hot {b['T_hot']:6.1f}")
