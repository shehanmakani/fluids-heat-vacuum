"""Sanity tests for the blend and heat-sink models. Run: python -m pytest -q"""
import math

import numpy as np
import pytest

from thermalform.physics import HeatSink, forced, natural
from thermalform.properties import (RHO_REF_T, blend_flash, blend_k, blend_kv40, blend_nu, blend_rho,
                                    evaluate, load_base_stocks)

bs = load_base_stocks()


def one(i):
    w = np.zeros(bs.n)
    w[i] = 1.0
    return w


def mix(**fr):
    w = np.zeros(bs.n)
    for k, v in fr.items():
        w[bs.ids.index(k)] = v
    return w


@pytest.mark.parametrize("i", range(bs.n))
def test_pure_component_round_trips(i):
    w = one(i)
    assert blend_kv40(w, bs) == pytest.approx(bs.kv40[i], rel=1e-9)
    assert blend_flash(w, bs) == pytest.approx(bs.flash[i], abs=1e-6)
    # Walther fit reproduces both anchor viscosities
    assert blend_nu(w, bs, 40.0) == pytest.approx(bs.kv40[i], rel=1e-6)
    assert blend_nu(w, bs, 100.0) == pytest.approx(bs.kv100[i], rel=1e-6)
    tref = RHO_REF_T.get(bs.ids[i], 20.0)
    assert blend_rho(w, bs, tref) == pytest.approx(bs.rho[i], rel=1e-9)


def test_viscosity_falls_with_temperature():
    for i in range(bs.n):
        assert bs.nu_T(70)[i] < bs.nu_T(40)[i] < bs.nu_T(20)[i]


def test_blend_viscosity_between_components():
    w = mix(PAO2=0.5, PAO6=0.5)
    kv = blend_kv40(w, bs)
    assert bs.kv40[0] < kv < bs.kv40[2]
    assert kv < (bs.kv40[0] + bs.kv40[2]) / 2       # log-type blending sits below the mean


def test_flash_dominated_by_lighter_component():
    w = mix(PAO2=0.5, NC1000=0.5)
    fp = blend_flash(w, bs)
    assert 155 < fp < (155 + 325) / 2


def test_filippov_below_linear():
    w = mix(PAO2=0.5, NC1000=0.5)
    lin = 0.5 * bs.k_T(40)[0] + 0.5 * bs.k_T(40)[bs.ids.index("NC1000")]
    assert blend_k(w, bs, 40) < lin


def test_forced_more_flow_lowers_resistance_and_costs_power():
    w = mix(PAO4=1.0)
    a, b = forced(w, bs, 2 / 60000), forced(w, bs, 8 / 60000)
    assert b["R_fluid"] < a["R_fluid"]
    assert b["P_elec_w"] > a["P_elec_w"]
    assert a["Re_Dh"] < 2300 and b["Re_Dh"] < 2300   # laminar, correlations valid


def test_energy_balance():
    w = mix(PAO4=1.0)
    r = forced(w, bs, 6 / 60000)
    assert r["T_base"] > 40.0
    assert r["fin_eff"] < 1.0


def test_natural_worse_than_forced():
    w = mix(PAO4=1.0)
    assert natural(w, bs)["R_nat"] > forced(w, bs, 6 / 60000)["R_fluid"]


def test_evaluate_runs_for_every_stock():
    for i in range(bs.n):
        p = evaluate(one(i), bs)
        assert all(not (isinstance(v, float) and math.isinf(v)) for v in p.values())


# ---------------- round 3: loop, dielectric, OCP ----------------
from thermalform.dielectric import eeff_air, fluid_effect, max_polar_fraction, width_for, z0
from thermalform.loop import loop, martin_xi, water
from thermalform.ocp import check, fom1
from thermalform.properties import blend_dk, fluid_at


def test_microstrip_matches_textbook_fr4():
    # 50 ohm on e_r 4.3 needs w/h of about 1.9-2.0 (standard calculators)
    assert 1.85 < width_for(50, 4.3) < 2.05


def test_fluid_cover_lowers_impedance_monotonically():
    zs = [fluid_effect(e)["z_fluid"] for e in (1.0, 2.1, 2.3, 3.0)]
    assert zs[0] == pytest.approx(50, abs=0.01)
    assert zs == sorted(zs, reverse=True)


def test_lichtenecker_limit_fraction_reproduces_limit():
    f = max_polar_fraction(2.08, 3.0, 2.3)
    assert math.exp(f * math.log(3.0) + (1 - f) * math.log(2.08)) == pytest.approx(2.3, rel=1e-9)


def test_hydrocarbon_dk_is_frequency_flat_and_esters_bounded():
    for i, c in enumerate(bs.chemistry):
        assert bs.dk_ghz[i] <= bs.dk_ghz_worst[i] + 1e-12
        if c in ("PAO", "mineral oil", "renewable hydrocarbon"):
            assert bs.dk_ghz[i] == bs.dk_ghz_worst[i]


def test_martin_friction_increases_with_chevron_angle():
    assert martin_xi(500, math.radians(30)) < martin_xi(500, math.radians(60))


def test_water_viscosity_reference():
    assert water(20)["mu"] == pytest.approx(1.0e-3, rel=0.03)


def test_loop_energy_balance_and_pump():
    r = loop(mix(PAO4=1.0), bs)
    assert r["T_water_in_max"] < 40.0
    assert 0 < r["pump_pct_it"] < 5
    # higher volumetric heat capacity -> less flow needed
    assert loop(mix(NC1000=1.0), bs)["flow_lpm"] < loop(mix(PAO2=1.0), bs)["flow_lpm"]


def test_fom1_matches_evonik_published_pao4():
    assert fom1(fluid_at(mix(PAO4=1.0), bs, 50.0)) == pytest.approx(46.0, rel=0.03)


def test_ocp_check_flags_ester_dk():
    c = check(mix(EF7010=0.52, P1976=0.48), bs)
    assert not c["pass_dk_worst"]
    assert check(mix(PAO4=1.0), bs)["pass_dk_worst"]
