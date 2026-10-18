"""Blend property models for single-phase immersion cooling fluids.

Blends are MASS fractions of base stocks. Pure-component properties come from
published data sheets (see data/base_stocks.csv, every value tagged with its basis).

v2 adds temperature dependence so heat transfer can be evaluated at the real
film temperature on a hot heat sink instead of at 40 C:
  - kinematic viscosity: Walther equation (ASTM D341) fitted to KV40 and KV100
  - density: linear thermal expansion
  - thermal conductivity, heat capacity: linear in T, anchored at each data
    sheet's reference temperature, slopes from the Synfluid PAO data
  - water saturation and an oxidation-exposure flag for the chemistry analyses
"""
from __future__ import annotations

import csv
import math
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

DATA = Path(__file__).resolve().parent.parent / "data"

# Reference temperatures of the published k and cp values (deg C), by stock id.
K_REF_T = {"PAO2": 38, "PAO4": 38, "PAO6": 38, "EF7010": 30, "EF3221": 30, "P1976": 30,
           "NC1000": 20, "HYV2": 40, "RHC8": 40}
CP_REF_T = {"PAO2": 40, "PAO4": 40, "PAO6": 40, "EF7010": 30, "EF3221": 30, "P1976": 30,
            "NC1000": 60, "HYV2": 40, "RHC8": 40}
RHO_REF_T = {"PAO2": 15.6, "PAO4": 15.6, "PAO6": 15.6, "NC1000": 20, "HYV2": 15}  # others 20
# Slopes. dk/dT from Synfluid PAO data 38->93 C (mean of three grades);
# dcp/dT from Synfluid PAO data 40->120 C (mean of three grades); thermal expansion typical of oils.
DK_DT = -(0.005 / 55 + 0.007 / 55 + 0.004 / 55) / 3      # W/m-K per K  (~ -9.7e-5)
DCP_DT = ((2.504 - 2.203) + (2.510 - 2.143) + (2.327 - 2.028)) / 3 / 80   # kJ/kg-K per K (~ +0.0040)
BETA = {"PAO": 7.0e-4, "mineral oil": 7.0e-4, "synthetic ester": 7.5e-4, "natural ester": 7.2e-4,
        "renewable hydrocarbon": 8.4e-4}  # 1/K; renewable HC from EC-120 density data (795.3 -> 761.8 kg/m3, 20 -> 70 C)


@dataclass
class BaseStocks:
    ids: list[str]
    names: list[str]
    chemistry: list[str]
    kv40: np.ndarray
    kv100: np.ndarray
    flash: np.ndarray
    pour: np.ndarray
    rho: np.ndarray
    k: np.ndarray
    cp: np.ndarray
    bdv: np.ndarray
    biodeg: np.ndarray
    biobased: np.ndarray
    price: np.ndarray
    wsat: np.ndarray          # water saturation at 20 C, ppm
    unsat: np.ndarray         # 1 = unsaturated (natural ester, oxidation-prone in air)
    dk_lowf: np.ndarray       # relative permittivity, low frequency (<= 100 kHz)
    dk_ghz: np.ndarray        # central estimate, 20-40 GHz
    dk_ghz_worst: np.ndarray  # worst case, 20-40 GHz (polar stocks: low-frequency value)
    flags: dict[str, list[str]]
    walther_a: np.ndarray = field(init=False)
    walther_b: np.ndarray = field(init=False)
    beta: np.ndarray = field(init=False)

    def __post_init__(self):
        # Walther / ASTM D341: log10(log10(nu + 0.7)) = A - B log10(T_K)
        t1, t2 = math.log10(313.15), math.log10(373.15)
        z1 = np.log10(np.log10(self.kv40 + 0.7))
        z2 = np.log10(np.log10(self.kv100 + 0.7))
        self.walther_b = (z1 - z2) / (t2 - t1)
        self.walther_a = z1 + self.walther_b * t1
        self.beta = np.array([BETA[c] for c in self.chemistry])

    @property
    def n(self) -> int:
        return len(self.ids)

    # ----- pure-component properties at temperature T (deg C) -----
    def nu_T(self, T):
        z = self.walther_a - self.walther_b * math.log10(T + 273.15)
        return 10 ** (10 ** z) - 0.7

    def rho_T(self, T):
        tref = np.array([RHO_REF_T.get(i, 20.0) for i in self.ids])
        return self.rho * (1 - self.beta * (T - tref))

    def k_T(self, T):
        tref = np.array([K_REF_T[i] for i in self.ids])
        return self.k + DK_DT * (T - tref)

    def cp_T(self, T):
        tref = np.array([CP_REF_T[i] for i in self.ids])
        return self.cp + DCP_DT * (T - tref)


def _f(x: str) -> float:
    return float(x) if x not in ("", None) else float("nan")


def load_base_stocks(path: Path | None = None, include_only: bool = True) -> BaseStocks:
    path = path or DATA / "base_stocks.csv"
    rows = list(csv.DictReader(open(path, newline="", encoding="utf-8")))
    if include_only:
        rows = [r for r in rows if r["include"] == "1"]
    flags = {}
    for r in rows:
        flags[r["id"]] = [f"{c}: {r[c]}" for c in
                          ("kv100_basis", "flash_basis", "pour_basis", "density_basis", "k_basis", "cp_basis", "price_basis",
                           "water_sat_basis", "dk_ghz_basis")
                          if "ESTIMATE" in r.get(c, "") or "PLACEHOLDER" in r.get(c, "")]
    a = lambda c: np.array([_f(r[c]) for r in rows])
    return BaseStocks(
        ids=[r["id"] for r in rows], names=[r["name"] for r in rows],
        chemistry=[r["chemistry"] for r in rows],
        kv40=a("kv40_cst"), kv100=a("kv100_cst"), flash=a("flash_coc_c"), pour=a("pour_c"),
        rho=a("density_kg_l"), k=a("k_w_mk"), cp=a("cp_kj_kgk"), bdv=a("bdv_kv"),
        biodeg=a("readily_biodegradable"), biobased=a("bio_based"), price=a("price_usd_l"),
        wsat=a("water_sat_ppm_20c"), unsat=a("unsaturated"),
        dk_lowf=a("dk_lowf"), dk_ghz=a("dk_ghz_central"), dk_ghz_worst=a("dk_ghz_worst"), flags=flags,
    )


def vol_fractions(w, bs, T=20.0):
    v = w / bs.rho_T(T)
    return v / v.sum()


# ---------- viscosity: Refutas blending number (mass basis) ----------
def _vbn(nu):
    return 14.534 * np.log(np.log(nu + 0.8)) + 10.975


def _vbn_inv(vbn):
    return math.exp(math.exp((vbn - 10.975) / 14.534)) - 0.8


def blend_nu(w, bs, T):
    return _vbn_inv(float(np.dot(w, _vbn(bs.nu_T(T)))))


def blend_kv40(w, bs):
    return _vbn_inv(float(np.dot(w, _vbn(bs.kv40))))


def blend_kv100(w, bs):
    return _vbn_inv(float(np.dot(w, _vbn(bs.kv100))))


# ---------- flash point: Hu-Burns blending index (volume basis) ----------
def _flash_index(t_c):
    return 10 ** (-6.1188 + 2414.0 / (t_c + 273.15 - 42.6))


def blend_flash(w, bs):
    bi = float(np.dot(vol_fractions(w, bs), _flash_index(bs.flash)))
    return 2414.0 / (math.log10(bi) + 6.1188) + 42.6 - 273.15


# ---------- mixture thermophysical properties at T ----------
def blend_rho(w, bs, T=20.0):
    return 1.0 / float(np.dot(w, 1.0 / bs.rho_T(T)))


def blend_cp(w, bs, T=40.0):
    return float(np.dot(w, bs.cp_T(T)))


def blend_k(w, bs, T=40.0):
    """Filippov equation (mass fractions), generalised pairwise."""
    kk = bs.k_T(T)
    lin = float(np.dot(w, kk))
    non = 0.0
    for i in range(bs.n):
        if w[i] == 0:
            continue
        for j in range(i + 1, bs.n):
            non += 0.72 * w[i] * w[j] * abs(kk[i] - kk[j])
    return lin - non


def fluid_at(w, bs, T):
    """All transport properties of blend w at temperature T (SI units)."""
    rho = blend_rho(w, bs, T) * 1000.0
    nu = blend_nu(w, bs, T) * 1e-6
    k = blend_k(w, bs, T)
    cp = blend_cp(w, bs, T) * 1000.0
    beta = float(np.dot(vol_fractions(w, bs, T), bs.beta))
    return {"rho": rho, "nu": nu, "mu": nu * rho, "k": k, "cp": cp, "Pr": nu * rho * cp / k, "beta": beta}


# ---------- screening-level indicators ----------
def blend_pour_worst(w, bs, tol=1e-3):
    vals = bs.pour[(w > tol)]
    vals = vals[~np.isnan(vals)]
    return float(vals.max()) if len(vals) else float("nan")


def blend_bdv_floor(w, bs, tol=1e-3):
    present = (w > tol) & ~np.isnan(bs.bdv)
    return float(bs.bdv[present].min()) if present.any() else float("nan")


def blend_water_sat(w, bs):
    """Water saturation of the blend at 20 C, ppm. Mass-weighted mixing is an
    approximation; real blends of polar esters and non-polar PAO can deviate."""
    return float(np.dot(w, bs.wsat))


def relative_saturation(water_ppm, wsat_ppm):
    return water_ppm / wsat_ppm


def oxidation_exposure(w, bs):
    """Mass fraction of unsaturated (natural-ester) stock: the part that polymerises in air."""
    return float(np.dot(w, bs.unsat))


def blend_dk(w, bs, which="ghz", T=25.0):
    """Lichtenecker logarithmic mixing on volume fractions: ln e = sum(phi_i ln e_i)."""
    e = {"ghz": bs.dk_ghz, "worst": bs.dk_ghz_worst, "lowf": bs.dk_lowf}[which]
    return float(math.exp(np.dot(vol_fractions(w, bs, T), np.log(e))))


def cost_per_l(w, bs):
    return float(np.dot(vol_fractions(w, bs), bs.price))


def mouromtseff_fom(kv_cst, rho_kg_l, k, cp_kj):
    rho = rho_kg_l * 1000.0
    mu = kv_cst * 1e-6 * rho
    return rho ** 0.8 * k ** 0.67 * (cp_kj * 1000.0) ** 0.33 / mu ** 0.47


def evaluate(w, bs, T=40.0) -> dict:
    """Screening properties of blend w (round-1 compatible, plus v2 fields)."""
    w = np.asarray(w, dtype=float)
    w = w / w.sum()
    kv40, kv100 = blend_kv40(w, bs), blend_kv100(w, bs)
    rho, k, cp = blend_rho(w, bs, T), blend_k(w, bs, T), blend_cp(w, bs, T)
    return {
        "kv40_cst": kv40,
        "kv100_cst": kv100,
        "kv20_cst": blend_nu(w, bs, 20.0),   # cold-start / cold-aisle viscosity
        "flash_c": blend_flash(w, bs),
        "density_kg_l": rho,
        "k_w_mk": k,
        "cp_kj_kgk": cp,
        "fom_turbulent": mouromtseff_fom(blend_nu(w, bs, T), rho, k, cp),
        "pour_worst_c": blend_pour_worst(w, bs),
        "bdv_floor_kv": blend_bdv_floor(w, bs),
        "water_sat_ppm": blend_water_sat(w, bs),
        "oxidation_exposure": oxidation_exposure(w, bs),
        "dk_ghz": blend_dk(w, bs, "ghz"),
        "dk_ghz_worst": blend_dk(w, bs, "worst"),
        "biodegradable_frac": float(np.dot(w, bs.biodeg)),
        "biobased_frac": float(np.dot(w, bs.biobased)),
        "cost_usd_l": cost_per_l(w, bs),
    }
