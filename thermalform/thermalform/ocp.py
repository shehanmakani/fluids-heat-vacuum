"""Open Compute Project single-phase immersion fluid checks.

Figures of merit and limits from the OCP Base Specification for Immersion Fluids
(2022-12-01); limits referenced to 25 C in that revision, 50 C announced for later:
  FOM1 natural convection   k * (beta * cp * rho^2 / (mu * k))^0.2813   Tier 1 > 35, Tier 2 > 45
  FOM2 forced, developing   k * (rho * cp / k)^(1/3)                     > 19
  FOM3 viscosity            mu < 0.015 Pa.s (pump check at 0 C)
  Dielectric constant       <= 2.3, 20 MHz - 40 GHz
  Flash point               > 150 C (base spec: closed cup, D93; requirements doc rev 2: COC, D92)
  Pour point                < -30 C
OCP ACS Immersion Requirements rev 2 adds: acidity (synthetic esters <= 0.03 mg KOH/g unused),
'mineral oils should be avoided'.
"""
from __future__ import annotations

import numpy as np

from .properties import blend_dk, blend_flash, blend_pour_worst, fluid_at

DK_LIMIT = 2.3
FLASH_LIMIT = 150.0
POUR_LIMIT = -30.0
FOM1_T1, FOM1_T2, FOM2_MIN, FOM3_MAX = 35.0, 45.0, 19.0, 0.015


def fom1(p):
    return p["k"] * (p["beta"] * p["cp"] * p["rho"] ** 2 / (p["mu"] * p["k"])) ** 0.2813


def fom2(p):
    return p["k"] * (p["rho"] * p["cp"] / p["k"]) ** (1 / 3)


def check(w, bs, chemistry_classes=None):
    w = np.asarray(w, float) / np.sum(w)
    out = {}
    for T in (25.0, 50.0):
        p = fluid_at(w, bs, T)
        out[f"fom1_{int(T)}"] = fom1(p)
        out[f"fom2_{int(T)}"] = fom2(p)
        out[f"fom3_{int(T)}"] = p["mu"]
    out["mu_0c"] = fluid_at(w, bs, 0.0)["mu"]
    out["dk_ghz"] = blend_dk(w, bs, "ghz")
    out["dk_ghz_worst"] = blend_dk(w, bs, "worst")
    out["flash_c"] = blend_flash(w, bs)
    out["pour_worst_c"] = blend_pour_worst(w, bs)
    tier = "Tier 2" if out["fom1_25"] > FOM1_T2 else ("Tier 1" if out["fom1_25"] > FOM1_T1 else "below Tier 1")
    out["fom1_tier_25"] = tier
    out["fom1_tier_50"] = "Tier 2" if out["fom1_50"] > FOM1_T2 else ("Tier 1" if out["fom1_50"] > FOM1_T1 else "below Tier 1")
    out["pass_fom2"] = out["fom2_25"] > FOM2_MIN
    out["pass_fom3_25"] = out["fom3_25"] < FOM3_MAX
    out["pass_fom3_50"] = out["fom3_50"] < FOM3_MAX
    out["pass_dk_central"] = out["dk_ghz"] <= DK_LIMIT
    out["pass_dk_worst"] = out["dk_ghz_worst"] <= DK_LIMIT
    out["pass_flash"] = out["flash_c"] > FLASH_LIMIT
    out["pass_pour"] = (out["pour_worst_c"] < POUR_LIMIT) if not np.isnan(out["pour_worst_c"]) else None
    mo = sum(w[i] for i in range(bs.n) if bs.chemistry[i] == "mineral oil")
    out["mineral_oil_frac"] = float(mo)
    return out
