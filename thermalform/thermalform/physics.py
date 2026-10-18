"""Heat-sink physics for single-phase immersion cooling.

A plate-fin heat sink on a 1 kW accelerator, immersed in the fluid. Two modes:

forced  - a fan/pump pushes fluid through the fin channels (as in the Intel/Submer
          forced-convection heat sink, ~18 W of fan power per ~800 W chip).
          Heat transfer: Teertstra et al. composite model for laminar developing
          flow in plate-fin channels. Pressure drop: Muzychka-Yovanovich apparent
          friction for developing flow between parallel plates + Kays-London
          entrance/exit losses.
natural - open bath, no fan; buoyancy drives flow up the vertical fin channels.
          Heat transfer: Bar-Cohen & Rohsenow composite correlation for symmetric
          isothermal parallel plates.

Fin efficiency: straight rectangular fin, adiabatic tip.
All properties are evaluated at the film temperature, found by iteration.
These are standard textbook correlations; absolute numbers depend on the
assumed geometry, but fluid-to-fluid ratios are robust to it.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, replace

import numpy as np
from scipy.optimize import brentq

from .properties import BaseStocks, fluid_at

G = 9.81


@dataclass(frozen=True)
class HeatSink:
    W: float = 0.090        # m, width across fins
    L: float = 0.090        # m, length along flow
    H: float = 0.030        # m, fin height
    t: float = 0.0006       # m, fin thickness
    b: float = 0.0015       # m, gap between fins
    k_fin: float = 390.0    # W/m-K, copper
    eta_mover: float = 0.35 # wire-to-fluid efficiency of a small immersion fan/pump

    @property
    def n_fins(self) -> int:
        return int((self.W + self.b) // (self.b + self.t))

    @property
    def n_ch(self) -> int:
        return self.n_fins - 1


def _fin_area_conductance(h, hs: HeatSink):
    m = math.sqrt(2 * h / (hs.k_fin * hs.t))
    eta = math.tanh(m * hs.H) / (m * hs.H)
    a_fin = hs.n_fins * 2 * hs.H * hs.L
    a_base = max(hs.W - hs.n_fins * hs.t, 0) * hs.L
    return h * (eta * a_fin + a_base), eta


def forced(w, bs: BaseStocks, vdot: float, hs: HeatSink = HeatSink(), Q: float = 1000.0, T_in: float = 40.0):
    """Forced convection through the fin channels at volumetric flow vdot (m3/s)."""
    return forced_props(lambda T: fluid_at(w, bs, T), vdot, hs, Q, T_in)


def forced_props(props_fn, vdot: float, hs: HeatSink = HeatSink(), Q: float = 1000.0, T_in: float = 40.0):
    """Same model, driven by any props_fn(T) -> property dict (used for sensitivity studies)."""
    T_film = T_in + 15.0
    for _ in range(8):
        p = props_fn(T_film)
        u = vdot / (hs.n_ch * hs.b * hs.H)
        re_b = u * hs.b / p["nu"]
        re_s = re_b * hs.b / hs.L
        nu_b = ((re_s * p["Pr"] / 2) ** -3 +
                (0.664 * math.sqrt(re_s) * p["Pr"] ** (1 / 3) * math.sqrt(1 + 3.65 / math.sqrt(re_s))) ** -3) ** (-1 / 3)
        h = nu_b * p["k"] / hs.b
        hA, eta = _fin_area_conductance(h, hs)
        mdot = p["rho"] * vdot
        r_conv = 1 / hA
        r_cal = 1 / (2 * mdot * p["cp"])
        T_base = T_in + Q * (r_conv + r_cal)
        T_fluid_mean = T_in + Q / (2 * mdot * p["cp"])
        T_new = min(0.5 * (T_base + T_fluid_mean), 150.0)
        if abs(T_new - T_film) < 0.01:
            break
        T_film = T_new
    dh = 2 * hs.b * hs.H / (hs.b + hs.H)
    re_dh = u * dh / p["nu"]
    l_star = hs.L / (dh * re_dh)
    f_re = math.sqrt((3.44 / math.sqrt(l_star)) ** 2 + 24.0 ** 2)   # Fanning
    dyn = 0.5 * p["rho"] * u ** 2
    sigma = hs.b / (hs.b + hs.t)
    dp = 4 * (f_re / re_dh) * (hs.L / dh) * dyn + (0.42 * (1 - sigma ** 2) + (1 - sigma) ** 2) * dyn
    p_hyd = vdot * dp
    return {
        "R_conv": r_conv, "R_cal": r_cal, "R_fluid": r_conv + r_cal, "T_base": T_base, "T_film": T_film,
        "h": h, "fin_eff": eta, "u": u, "Re_Dh": re_dh, "dp_pa": dp, "P_hyd_w": p_hyd,
        "P_elec_w": p_hyd / hs.eta_mover, "vdot_lpm": vdot * 60000, "nu_film_cst": p["nu"] * 1e6,
    }


def forced_at_power(w, bs, P_elec: float, hs: HeatSink = HeatSink(), Q=1000.0, T_in=40.0):
    """Flow that consumes exactly P_elec watts of fan/pump power, and the result."""
    f = lambda lv: forced(w, bs, 10 ** lv, hs, Q, T_in)["P_elec_w"] - P_elec
    lv = brentq(f, -5.5, -2.0, xtol=1e-6)
    return forced(w, bs, 10 ** lv, hs, Q, T_in)


def forced_power_for_R(w, bs, R_target: float, hs: HeatSink = HeatSink(), Q=1000.0, T_in=40.0):
    """Fan/pump power needed for fluid-side resistance R_target (K/W). NaN if unreachable."""
    f = lambda lv: forced(w, bs, 10 ** lv, hs, Q, T_in)["R_fluid"] - R_target
    try:
        lv = brentq(f, -5.5, -2.0, xtol=1e-7)
    except ValueError:
        return float("nan"), None
    r = forced(w, bs, 10 ** lv, hs, Q, T_in)
    return r["P_elec_w"], r


def natural(w, bs, hs: HeatSink = HeatSink(b=0.004), Q: float = 1000.0, T_bath: float = 40.0):
    """Open-bath natural convection; fins vertical, channel height = hs.L."""
    return natural_props(lambda T: fluid_at(w, bs, T), hs, Q, T_bath)


def natural_props(props_fn, hs: HeatSink = HeatSink(b=0.004), Q: float = 1000.0, T_bath: float = 40.0):
    dT = 30.0
    for _ in range(60):
        T_film = T_bath + dT / 2
        p = props_fn(T_film)
        alpha = p["k"] / (p["rho"] * p["cp"])
        ra_b = G * p["beta"] * dT * hs.b ** 3 / (p["nu"] * alpha)
        el = ra_b * hs.b / hs.L
        nu_b = (576.0 / el ** 2 + 2.873 / math.sqrt(el)) ** -0.5
        h = nu_b * p["k"] / hs.b
        hA, eta = _fin_area_conductance(h, hs)
        dT_new = Q / hA
        if abs(dT_new - dT) < 0.01:
            break
        dT = 0.5 * dT + 0.5 * dT_new
    return {"R_nat": dT / Q, "dT_base": dT, "h": h, "fin_eff": eta, "El": el}


def best_gap(fn, gaps=None):
    """Scan fin gap; fn(b) -> resistance. Returns (b_opt, R_opt, curve)."""
    gaps = gaps if gaps is not None else np.round(np.arange(0.5, 6.01, 0.1), 2) * 1e-3
    curve = []
    for b in gaps:
        try:
            curve.append((b, fn(b)))
        except (ValueError, ZeroDivisionError, OverflowError):
            continue
    b_opt, r_opt = min(curve, key=lambda c: c[1])
    return b_opt, r_opt, curve
