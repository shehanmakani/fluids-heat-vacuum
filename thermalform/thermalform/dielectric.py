"""Signal-integrity effect of the immersion fluid on outer-layer PCB traces.

Inner-layer striplines are fully embedded in laminate and do not see the fluid.
Outer-layer microstrips (and connector/cable dielectrics) have the fluid, instead of
air, above them, which raises the effective permittivity, lowers impedance and slows
propagation.

Microstrip model: Hammerstad-Jensen closed forms for the air-covered line; the
fluid cover is handled with the standard filling-factor approximation
    e_eff,fluid = q * e_r + (1 - q) * e_fluid,  q = (e_eff,air - 1) / (e_r - 1)
which is first-order but adequate for ranking fluids.
"""
from __future__ import annotations

import math

C_IN_PER_NS = 11.8028     # speed of light, inches per ns
ETA0 = 376.73             # ohm


def eeff_air(er, u):
    """Hammerstad-Jensen effective permittivity, u = w/h."""
    a = 1 + (1 / 49) * math.log((u ** 4 + (u / 52) ** 2) / (u ** 4 + 0.432)) + (1 / 18.7) * math.log(1 + (u / 18.1) ** 3)
    b = 0.564 * ((er - 0.9) / (er + 3)) ** 0.053
    return (er + 1) / 2 + (er - 1) / 2 * (1 + 10 / u) ** (-a * b)


def z0_vacuum(u):
    """Hammerstad-Jensen characteristic impedance of the line with e_r = 1."""
    f = 6 + (2 * math.pi - 6) * math.exp(-((30.666 / u) ** 0.7528))
    return ETA0 / (2 * math.pi) * math.log(f / u + math.sqrt(1 + (2 / u) ** 2))


def z0(er, u, e_cover=1.0):
    ea = eeff_air(er, u)
    q = (ea - 1) / (er - 1)
    ee = q * er + (1 - q) * e_cover
    return z0_vacuum(u) / math.sqrt(ee), ee, q


def width_for(z_target, er, e_cover=1.0):
    lo, hi = 0.05, 20.0
    for _ in range(80):
        m = math.sqrt(lo * hi)
        z = z0(er, m, e_cover)[0]
        lo, hi = (m, hi) if z > z_target else (lo, m)
    return math.sqrt(lo * hi)


def fluid_effect(e_fluid, er=3.5, z_design=50.0, tol=0.10):
    """A trace designed for z_design in air, then immersed in a fluid of permittivity e_fluid."""
    u = width_for(z_design, er)
    z_air, ee_air, q = z0(er, u, 1.0)
    z_f, ee_f, _ = z0(er, u, e_fluid)
    gamma = abs(z_f - z_design) / (z_f + z_design)
    return {
        "w_over_h": u, "fill_q": q, "z_fluid": z_f, "dz_pct": (z_f / z_design - 1) * 100,
        "return_loss_db": -20 * math.log10(gamma) if gamma > 0 else float("inf"),
        "delay_ps_per_in_air": math.sqrt(ee_air) / C_IN_PER_NS * 1000,
        "delay_ps_per_in_fluid": math.sqrt(ee_f) / C_IN_PER_NS * 1000,
        "skew_ps_per_in": (math.sqrt(ee_f) - math.sqrt(ee_air)) / C_IN_PER_NS * 1000,
        "within_tol": abs(z_f / z_design - 1) <= tol,
    }


def redesign_width_change(e_fluid, er=3.5, z_design=50.0):
    """Width change an OEM must make to restore z_design in the fluid (%)."""
    return (width_for(z_design, er, e_fluid) / width_for(z_design, er, 1.0) - 1) * 100


def max_polar_fraction(e_nonpolar, e_polar, limit=2.3):
    """Largest volume fraction of a polar stock in a non-polar one (Lichtenecker) under a Dk limit."""
    if e_polar <= limit:
        return 1.0
    return (math.log(limit) - math.log(e_nonpolar)) / (math.log(e_polar) - math.log(e_nonpolar))
