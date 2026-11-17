"""Unit tests for nf_properties.py. Run:  python -m pytest -v"""
import math
import pytest
from nf_properties import *

WATER = BaseFluid("water", 996.5, 4181.0, 0.613, 8.55e-4)


def cuo(phi):
    return Particle("CuO", phi, 6350.0, 535.6, 69.0, "spherical")


# --- Spec Section 5 test vector --------------------------------------------
def test_spec_vector():
    s = case_summary(WATER, [cuo(0.002)], 10000)
    assert s["rho_NF"] == pytest.approx(1007.21, rel=1e-4)
    assert s["cp_NF"] == pytest.approx(4135.0, rel=1e-3)
    assert s["k_NF"] == pytest.approx(0.61659, rel=1e-4)
    assert s["mu_NF"] == pytest.approx(8.593e-4, rel=1e-4)
    assert s["Pr_NF"] == pytest.approx(5.763, rel=1e-3)
    assert s["alpha_BF"] == pytest.approx(1.4713e-7, rel=1e-3)
    assert s["alpha_p"] == pytest.approx(2.0287e-5, rel=1e-3)
    assert s["alpha_r"] == pytest.approx(137.9, rel=1e-3)
    assert s["V"] == pytest.approx(0.3709, rel=1e-3)
    assert s["mdot"] == pytest.approx(0.1552, rel=2e-3)
    assert s["f"] == pytest.approx(0.030565, rel=1e-4)
    assert s["Nu_BF"] == pytest.approx(72.4, rel=3e-3)
    assert s["eta_eq4"] == pytest.approx(1.210, rel=2e-3)
    assert s["Nu_eq4"] == pytest.approx(87.6, rel=3e-3)
    assert s["Nu_eq3"] == pytest.approx(82.5, rel=3e-3)


# --- phi_ref justification: split-batch consistency -------------------------
@pytest.mark.parametrize("shape", list(SHAPES))
def test_split_batch_total_phi(shape):
    """Two identical batches at phi/2 each must equal one particle at phi."""
    phi = 0.01
    one = Particle("X", phi, 6350.0, 535.6, 69.0, shape)
    a = Particle("X", phi / 2, 6350.0, 535.6, 69.0, shape)
    b = Particle("X", phi / 2, 6350.0, 535.6, 69.0, shape)
    r1 = nanofluid_properties(WATER, [one])
    r2 = nanofluid_properties(WATER, [a, b])
    assert r2.k == pytest.approx(r1.k, rel=1e-12)
    assert r2.mu == pytest.approx(r1.mu, rel=1e-12)
    assert r2.rho == pytest.approx(r1.rho, rel=1e-12)
    assert r2.cp == pytest.approx(r1.cp, rel=1e-12)
    assert r2.alpha_r == pytest.approx(r1.alpha_r, rel=1e-12)


def test_split_batch_own_phi_fails():
    """Documents why 'own' fraction is rejected: it halves the enhancement."""
    phi = 0.01
    one = nanofluid_properties(WATER, [cuo(phi)])
    split = nanofluid_properties(WATER, [cuo(phi / 2), cuo(phi / 2)], phi_ref_mode="own")
    enh_one = one.k / WATER.k - 1
    enh_split = split.k / WATER.k - 1
    assert enh_split == pytest.approx(0.5 * enh_one, rel=0.02)


def test_ternary_weights_and_limits():
    ps = [Particle("a", 0.003, 6350, 535.6, 69, "spherical"),
          Particle("b", 0.003, 3970, 765, 40, "cylinder"),
          Particle("c", 0.004, 2220, 710, 3000, "platelet")]
    r = nanofluid_properties(WATER, ps)
    assert r.phi == pytest.approx(0.01)
    assert r.n_eff == pytest.approx(0.3 * 3 + 0.3 * 4.9 + 0.4 * 5.7)
    assert min(r.k_components) <= r.k <= max(r.k_components)


def test_input_validation():
    with pytest.raises(ValueError):
        nanofluid_properties(WATER, [])
    with pytest.raises(ValueError):
        nanofluid_properties(WATER, [cuo(0.001)] * 4)
    with pytest.raises(ValueError):
        nanofluid_properties(WATER, [cuo(0.0)])


def test_limit_phi_to_zero():
    r = nanofluid_properties(WATER, [cuo(1e-9)])
    assert r.k == pytest.approx(WATER.k, rel=1e-6)
    assert r.mu == pytest.approx(WATER.mu, rel=1e-6)


# --- Correlations and geometry ----------------------------------------------
def test_geometry_constants():
    assert math.pi / 4 * D_I ** 2 == pytest.approx(FLOW_AREA, rel=1e-4)
    assert math.pi * D_O * 2.0 / 8 == pytest.approx(SECTOR_AREA_OUTER, rel=1e-4)
    assert math.pi * D_I * 2.0 == pytest.approx(HEATED_AREA_INNER, rel=1e-4)
    assert FLUX_RATIO == pytest.approx(1.087, abs=1e-3)


def test_wall_distance_y1():
    """Spec Sec. 4: y1 (y+=1) = 37.2 um at Re 1e4 and 16.6 um at Re 2.5e4."""
    for Re, expect in [(10000, 37.2e-6), (25000, 16.6e-6)]:
        y1 = D_I / (Re * math.sqrt(friction_factor(Re) / 8))
        assert y1 == pytest.approx(expect, rel=1.5e-2)  # spec rounds 16.8 um to 16.6 um


def test_reynolds_roundtrip():
    nf = nanofluid_properties(WATER, [cuo(0.002)])
    V = inlet_velocity(18000, nf.rho, nf.mu)
    assert reynolds(nf.rho, V, nf.mu) == pytest.approx(18000)


# --- Flux scaling -----------------------------------------------------------
def test_flux_scaling_water():
    pat = placeholder_pattern()
    assert sum(pat) / 8 == pytest.approx(1.0)
    q10 = scale_fluxes(pat, 10000, WATER.k, WATER.Pr)
    q25 = scale_fluxes(pat, 25000, WATER.k, WATER.Pr)
    assert 8.0e3 < sum(q10) / 8 < 1.0e4       # ~9 kW/m2
    assert 1.8e4 < sum(q25) / 8 < 2.3e4       # ~20 kW/m2


def test_placeholder_pattern_shape():
    p = placeholder_pattern()
    assert p[0] == pytest.approx(p[7])          # Q1 = Q8
    assert p[1] == pytest.approx(p[6])          # Q2 = Q7
    assert p[0] == pytest.approx(2.828, abs=2e-3)
    assert p[1] == pytest.approx(1.172, abs=2e-3)
    assert all(abs(v) < 1e-9 for v in p[2:6])   # Q3..Q6 = 0


# --- Temperature tables -----------------------------------------------------
def test_property_table():
    funcs = dict(rho=lambda T: 1000 - 0.3 * (T - 300), cp=lambda T: 2000 + 3 * (T - 300),
                 k=lambda T: 0.13 - 1e-4 * (T - 300), mu=lambda T: 0.05 * math.exp(-0.03 * (T - 300)))
    tbl = build_property_table(funcs, [cuo(0.002)], 300, 340, 20)
    assert len(tbl) == 20
    assert tbl[0]["T"] == 300 and tbl[-1]["T"] == 340
    assert tbl[-1]["mu"] < tbl[0]["mu"]
    with pytest.raises(ValueError):
        build_property_table(funcs, [cuo(0.002)], 300, 340, 10)


def test_reference_T():
    assert estimate_reference_T(300, 10) == pytest.approx(307)
