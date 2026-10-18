"""Facility-level cooling loop for 1 MW of single-phase immersion.

    tank (1 MW IT) --> pump --> plate heat exchanger (CDU) --> back to tank
                                      ^ facility water (dry cooler / chiller)

Design rule held constant across fluids: the tank receives fluid at T_SUPPLY and
returns it DT_TANK warmer, so the volumetric flow each fluid needs is
Q / (rho * cp * DT_TANK). Fluids with higher volumetric heat capacity need less flow.

Plate heat exchanger: chevron plates, VDI Heat Atlas form of the Martin (1996)
correlation for both friction and heat transfer, with the (mu/mu_w)^(1/6) wall
correction on the oil side. Counterflow epsilon-NTU.
Piping: Darcy-Weisbach (64/Re laminar, Haaland turbulent) plus summed minor losses.

Outputs per fluid: loop flow, pressure drop, pump power, annual pump cost, and the
warmest facility-water inlet that still delivers T_SUPPLY to the tank.

Geometry, loss coefficients and pump efficiency are stated assumptions; absolute
values depend on them, fluid-to-fluid differences much less so.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

from .properties import fluid_at

Q_IT = 1.0e6          # W
T_SUPPLY = 40.0       # C, fluid to tank (and to the heat sinks)
DT_TANK = 10.0        # K, tank temperature rise
PRICE_KWH = 0.1629    # USD/kWh, EIA New Jersey all-sector average, 2024
HOURS = 8760


@dataclass(frozen=True)
class PlateHX:
    width: float = 0.45       # m, effective plate width
    length: float = 1.20      # m, effective port-to-port length
    b: float = 0.0025         # m, pressing depth (channel gap)
    phi_deg: float = 60.0     # chevron angle vs. flow direction
    enlarge: float = 1.17     # surface enlargement factor
    n_plates: int = 121       # -> 60 channels per side
    t_plate: float = 0.0005   # m
    k_plate: float = 16.0     # W/m-K, stainless
    water_kg_s: float = 30.0  # facility water flow (about 8 K rise at 1 MW)

    @property
    def n_ch(self):
        return (self.n_plates - 1) // 2

    @property
    def dh(self):
        return 2 * self.b / self.enlarge

    @property
    def area(self):            # total heat-transfer area, m2
        return (self.n_plates - 2) * self.width * self.length * self.enlarge


@dataclass(frozen=True)
class Piping:
    d: float = 0.15           # m, inner diameter (6 in)
    length: float = 60.0      # m, supply + return
    k_minor: float = 12.0     # elbows, valves, strainer, tank manifold
    eta_pump: float = 0.65
    rough: float = 4.5e-5     # m, commercial steel


# ---------------- water properties ----------------
def water(T):
    mu = 2.414e-5 * 10 ** (247.8 / (T + 133.15))
    rho = 1000.0 - 0.0178 * abs(T - 4.0) ** 1.7
    k = 0.5706 + 0.00175 * T - 6.9e-6 * T ** 2
    cp = 4180.0
    return {"rho": rho, "mu": mu, "nu": mu / rho, "k": k, "cp": cp, "Pr": mu * cp / k}


# ---------------- Martin / VDI chevron correlation ----------------
def martin_xi(re, phi):
    if re < 2000:
        x0, x1 = 64.0 / re, 597.0 / re + 3.85
    else:
        x0, x1 = (1.8 * math.log10(re) - 1.5) ** -2, 39.0 / re ** 0.289
    c = math.cos(phi)
    inv = c / math.sqrt(0.045 * math.tan(phi) + 0.09 * math.sin(phi) + x0 / c) + (1 - c) / math.sqrt(3.8 * x1)
    return inv ** -2


def martin_side(vdot, p, hx: PlateHX, visc_ratio=1.0):
    """Heat transfer coefficient and pressure drop for one side of the PHE."""
    phi = math.radians(hx.phi_deg)
    u = vdot / (hx.n_ch * hx.width * hx.b)
    re = u * hx.dh / p["nu"]
    xi = martin_xi(re, phi)
    nu_ = 0.122 * p["Pr"] ** (1 / 3) * visc_ratio ** (1 / 6) * (xi * re ** 2 * math.sin(2 * phi)) ** 0.374
    h = nu_ * p["k"] / hx.dh
    dp = xi * hx.length / hx.dh * p["rho"] * u ** 2 / 2
    return {"h": h, "dp": dp, "Re": re, "u": u, "xi": xi, "Nu": nu_}


def counterflow_eps(ntu, cr):
    if abs(1 - cr) < 1e-9:
        return ntu / (1 + ntu)
    e = math.exp(-ntu * (1 - cr))
    return (1 - e) / (1 - cr * e)


def pipe_dp(vdot, p, pp: Piping):
    a = math.pi * pp.d ** 2 / 4
    u = vdot / a
    re = u * pp.d / p["nu"]
    if re < 2300:
        f = 64.0 / re
    else:
        f = (-1.8 * math.log10((pp.rough / pp.d / 3.7) ** 1.11 + 6.9 / re)) ** -2
    dyn = p["rho"] * u ** 2 / 2
    return f * pp.length / pp.d * dyn + pp.k_minor * dyn, re


def loop(w, bs, hx: PlateHX = PlateHX(), pp: Piping = Piping(), Q=Q_IT, T_sup=T_SUPPLY, dT=DT_TANK):
    T_ret = T_sup + dT
    T_mean = T_sup + dT / 2
    pf = fluid_at(w, bs, T_mean)
    vdot = Q / (pf["rho"] * pf["cp"] * dT)
    m_f = pf["rho"] * vdot
    C_h = m_f * pf["cp"]
    # water side: guess inlet, iterate on wall temperature for the viscosity correction
    T_w_in = T_sup - 5.0
    for _ in range(30):
        T_w_out = T_w_in + Q / (hx.water_kg_s * 4180.0)
        T_wm = 0.5 * (T_w_in + T_w_out)
        pw = water(T_wm)
        vw = hx.water_kg_s / pw["rho"]
        T_wall = 0.5 * (T_mean + T_wm)
        mu_wall_oil = fluid_at(w, bs, T_wall)["mu"]
        sf = martin_side(vdot, pf, hx, visc_ratio=pf["mu"] / mu_wall_oil)
        sw = martin_side(vw, pw, hx, visc_ratio=pw["mu"] / water(T_wall)["mu"])
        U = 1.0 / (1 / sf["h"] + hx.t_plate / hx.k_plate + 1 / sw["h"])
        C_c = hx.water_kg_s * pw["cp"]
        cmin, cmax = min(C_h, C_c), max(C_h, C_c)
        ntu = U * hx.area / cmin
        eps = counterflow_eps(ntu, cmin / cmax)
        new = T_ret - Q / (eps * cmin)
        if abs(new - T_w_in) < 1e-4:
            break
        T_w_in = new
    dp_pipe, re_pipe = pipe_dp(vdot, pf, pp)
    dp = sf["dp"] + dp_pipe
    p_pump = vdot * dp / pp.eta_pump
    return {
        "flow_lpm": vdot * 60000, "flow_m3s": vdot, "T_water_in_max": T_w_in,
        "approach_c": T_sup - T_w_in, "U": U, "h_fluid": sf["h"], "h_water": sw["h"],
        "Re_hx": sf["Re"], "Re_pipe": re_pipe, "dp_hx_kpa": sf["dp"] / 1e3, "dp_pipe_kpa": dp_pipe / 1e3,
        "dp_kpa": dp / 1e3, "pump_kw": p_pump / 1e3, "pump_pct_it": p_pump / Q * 100,
        "pump_usd_yr": p_pump / 1e3 * HOURS * PRICE_KWH, "eps": eps, "ntu": ntu,
    }
