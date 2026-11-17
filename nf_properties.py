"""
nf_properties.py
Single-phase ternary nanofluid property calculator (Phase 1).

Follows Tables 3 and 4 of Upadhyay, Chandra and Sarkar, Appl. Therm. Eng. 190 (2021)
116469, as summarised in the build specification, extended to up to three particles.
Pure Python, no ANSYS needed. SI units throughout. Volume fractions are fractions.

Working assumptions (pending supervisor):
  * phi_ref = TOTAL phi for every particle (see tests: split-batch consistency).
  * alpha_i uses the PARTICLE conductivity k_p (matches the spec test vector).
  * Particle properties are constant (not temperature dependent).
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Sequence

# ----------------------------------------------------------------------------
# Constants
# ----------------------------------------------------------------------------
D_I = 0.023                      # inner diameter [m]
D_O = 0.025                      # outer diameter [m]
FLOW_AREA = 4.1548e-4            # [m2], pi/4 * D_I^2
SECTOR_AREA_OUTER = 0.019635     # [m2] heated outer area per 45 deg sector
HEATED_AREA_INNER = 0.14451      # [m2]
FLUX_RATIO = D_O / D_I           # 1.087: mean inner flux = FLUX_RATIO * mean(q_k)

# Shape table (spec Section 5): name -> (sphericity, n, A1, A2)
SHAPES: Dict[str, Dict[str, float]] = {
    "spherical": dict(sphericity=1.00, n=3.0, A1=2.5, A2=6.2),
    "brick":     dict(sphericity=0.81, n=3.7, A1=1.9, A2=471.4),
    "cylinder":  dict(sphericity=0.62, n=4.9, A1=13.5, A2=904.4),
    "platelet":  dict(sphericity=0.52, n=5.7, A1=37.1, A2=612.6),
}

# Paper's water-based correlation constants (Eq. 4)
PAPER_C = dict(C1=0.93, C2=-0.059, C3=-0.026, C4=0.357)


# ----------------------------------------------------------------------------
# Data classes
# ----------------------------------------------------------------------------
@dataclass
class BaseFluid:
    name: str
    rho: float   # kg/m3
    cp: float    # J/kg/K
    k: float     # W/m/K
    mu: float    # Pa.s

    @property
    def alpha(self) -> float:
        return self.k / (self.rho * self.cp)

    @property
    def Pr(self) -> float:
        return self.mu * self.cp / self.k


@dataclass
class Particle:
    name: str
    phi: float          # volume fraction contributed by this particle
    rho: float
    cp: float
    k: float
    shape: str = "spherical"

    @property
    def n(self) -> float:
        return SHAPES[self.shape]["n"]

    @property
    def A1(self) -> float:
        return SHAPES[self.shape]["A1"]

    @property
    def A2(self) -> float:
        return SHAPES[self.shape]["A2"]

    @property
    def alpha(self) -> float:
        """Particle thermal diffusivity k_p / (rho_p cp_p)."""
        return self.k / (self.rho * self.cp)


@dataclass
class NanofluidProps:
    phi: float
    rho: float
    cp: float
    k: float
    mu: float
    Pr: float
    alpha_r: float
    n_eff: float
    alpha_p: float
    alpha_bf: float
    k_components: List[float] = field(default_factory=list)
    mu_components: List[float] = field(default_factory=list)


# ----------------------------------------------------------------------------
# Single-particle building blocks
# ----------------------------------------------------------------------------
def hamilton_crosser(k_bf: float, k_p: float, n: float, phi_ref: float) -> float:
    """Hamilton-Crosser conductivity of a single-particle suspension at phi_ref."""
    num = k_p + (n - 1.0) * k_bf - (n - 1.0) * phi_ref * (k_bf - k_p)
    den = k_p + (n - 1.0) * k_bf + phi_ref * (k_bf - k_p)
    return k_bf * num / den


def viscosity_single(mu_bf: float, A1: float, A2: float, phi_ref: float) -> float:
    """Single-particle viscosity: mu_BF (1 + A1 phi + A2 phi^2)."""
    return mu_bf * (1.0 + A1 * phi_ref + A2 * phi_ref ** 2)


# ----------------------------------------------------------------------------
# Mixture properties
# ----------------------------------------------------------------------------
def nanofluid_properties(base: BaseFluid,
                         particles: Sequence[Particle],
                         phi_ref_mode: str = "total") -> NanofluidProps:
    """
    Compute nanofluid properties for 1 to 3 particles.

    phi_ref_mode:
      "total" (default) - phi_ref = total phi for every particle
      "own"             - phi_ref = particle's own phi_i (kept only for comparison;
                          fails the split-batch consistency test)
    """
    if not 1 <= len(particles) <= 3:
        raise ValueError("Between 1 and 3 particles are supported")
    if phi_ref_mode not in ("total", "own"):
        raise ValueError("phi_ref_mode must be 'total' or 'own'")

    phi = sum(p.phi for p in particles)
    if phi <= 0:
        raise ValueError("Total volume fraction must be > 0")

    rho = (1.0 - phi) * base.rho + sum(p.phi * p.rho for p in particles)
    cp = ((1.0 - phi) * base.rho * base.cp
          + sum(p.phi * p.rho * p.cp for p in particles)) / rho

    k_parts, mu_parts = [], []
    for p in particles:
        ref = phi if phi_ref_mode == "total" else p.phi
        k_parts.append(hamilton_crosser(base.k, p.k, p.n, ref))
        mu_parts.append(viscosity_single(base.mu, p.A1, p.A2, ref))
    w = [p.phi / phi for p in particles]
    k = sum(wi * ki for wi, ki in zip(w, k_parts))
    mu = sum(wi * mi for wi, mi in zip(w, mu_parts))

    alpha_p = sum(wi * p.alpha for wi, p in zip(w, particles))
    n_eff = sum(wi * p.n for wi, p in zip(w, particles))
    alpha_r = alpha_p / base.alpha

    return NanofluidProps(phi=phi, rho=rho, cp=cp, k=k, mu=mu, Pr=mu * cp / k,
                          alpha_r=alpha_r, n_eff=n_eff, alpha_p=alpha_p,
                          alpha_bf=base.alpha, k_components=k_parts,
                          mu_components=mu_parts)


# ----------------------------------------------------------------------------
# Flow quantities
# ----------------------------------------------------------------------------
def inlet_velocity(Re: float, rho: float, mu: float, D: float = D_I) -> float:
    return Re * mu / (rho * D)


def mass_flow(rho: float, V: float, area: float = FLOW_AREA) -> float:
    return rho * V * area


def reynolds(rho: float, V: float, mu: float, D: float = D_I) -> float:
    return rho * V * D / mu


# ----------------------------------------------------------------------------
# Correlations
# ----------------------------------------------------------------------------
def friction_factor(Re: float) -> float:
    """Smooth-pipe Darcy factor (Petukhov): f = (1.84 log10 Re - 1.64)^-2 (spec Sec. 4)."""
    return (1.84 * math.log10(Re) - 1.64) ** -2


def gnielinski(Re: float, Pr: float) -> float:
    """Gnielinski Nusselt number with the spec's friction factor."""
    f = friction_factor(Re)
    f8 = f / 8.0
    return f8 * (Re - 1000.0) * Pr / (1.0 + 12.7 * math.sqrt(f8) * (Pr ** (2.0 / 3.0) - 1.0))


def paper_eq3(Re: float, Pr: float) -> float:
    """Paper Eq. 3 (as quoted in the spec): Nu = 0.074 Re^0.72 Pr^0.22."""
    return 0.074 * Re ** 0.72 * Pr ** 0.22


def eta_paper(alpha_r: float, phi: float, n_eff: float,
              C1: float = PAPER_C["C1"], C2: float = PAPER_C["C2"],
              C3: float = PAPER_C["C3"], C4: float = PAPER_C["C4"]) -> float:
    """Correction factor eta = C1 alpha_r^C2 phi^C3 n^C4 (paper Eq. 4 form)."""
    return C1 * alpha_r ** C2 * phi ** C3 * n_eff ** C4


def paper_eq4(Re: float, Pr: float, alpha_r: float, phi: float, n_eff: float) -> float:
    """Nu_NF from Eq. 4 = eta * Gnielinski(Re_NF, Pr_NF)."""
    return eta_paper(alpha_r, phi, n_eff) * gnielinski(Re, Pr)


# ----------------------------------------------------------------------------
# Heat flux scaling
# ----------------------------------------------------------------------------
def scale_fluxes(pattern: Sequence[float], Re: float, base_or_nf_k: float,
                 Pr: float, dT_target: float = 5.0, D: float = D_I) -> List[float]:
    """
    Scale a flux pattern (any positive/zero numbers, ratios kept) so that the mean
    wall-to-bulk temperature difference is about dT_target.

      h        = Nu_Gnielinski * k / D
      mean q_k = dT_target * h / FLUX_RATIO       (FLUX_RATIO = D_o/D_i = 1.087)
    The pattern is normalised to mean 1 first.
    """
    h = gnielinski(Re, Pr) * base_or_nf_k / D
    mean_q = dT_target * h / FLUX_RATIO
    m = sum(pattern) / len(pattern)
    if m <= 0:
        raise ValueError("Flux pattern must have a positive mean")
    return [mean_q * q / m for q in pattern]


def placeholder_pattern(theta0_deg: float = 0.0) -> List[float]:
    """
    Placeholder 8-sector pattern (PENDING supervisor's real values): clipped cosine
    peaking between Q8 and Q1, normalised to mean 1. Sector k centre = (k-1)*45 deg.
    The peak sits at -22.5 deg (midway between Q8 at 315 and Q1 at 0).
    """
    vals = []
    for k in range(1, 9):
        theta = (k - 1) * 45.0 + theta0_deg
        vals.append(max(0.0, math.cos(math.radians(theta + 22.5))))
    m = sum(vals) / 8.0
    return [v / m for v in vals]


# ----------------------------------------------------------------------------
# Temperature-dependent tables (oils, salts)
# ----------------------------------------------------------------------------
def build_property_table(base_funcs: Dict[str, Callable[[float], float]],
                         particles: Sequence[Particle],
                         T_start: float, T_end: float, n_points: int = 20,
                         phi_ref_mode: str = "total") -> List[Dict[str, float]]:
    """
    Evaluate nanofluid properties on a temperature grid.

    base_funcs: dict with callables 'rho','cp','k','mu' taking T [K].
    Returns a list of dicts: T, rho, cp, k, mu (ready for Fluent piecewise-linear tables).
    Use T_start = T_in and T_end = T_in + expected rise + 10 K (spec Section 5).
    """
    if n_points < 20:
        raise ValueError("Spec requires at least 20 temperature points")
    rows = []
    for i in range(n_points):
        T = T_start + (T_end - T_start) * i / (n_points - 1)
        b = BaseFluid("tbl", base_funcs["rho"](T), base_funcs["cp"](T),
                      base_funcs["k"](T), base_funcs["mu"](T))
        nf = nanofluid_properties(b, particles, phi_ref_mode)
        rows.append(dict(T=T, rho=nf.rho, cp=nf.cp, k=nf.k, mu=nf.mu))
    return rows


def estimate_reference_T(T_in: float, bulk_rise: float, fraction: float = 0.7) -> float:
    """Estimate bulk temperature in the evaluation window: T_in + 70% of bulk rise."""
    return T_in + fraction * bulk_rise


# ----------------------------------------------------------------------------
# Convenience: full case summary
# ----------------------------------------------------------------------------
def case_summary(base: BaseFluid, particles: Sequence[Particle], Re: float) -> dict:
    nf = nanofluid_properties(base, particles)
    V = inlet_velocity(Re, nf.rho, nf.mu)
    return dict(
        phi=nf.phi, rho_NF=nf.rho, cp_NF=nf.cp, k_NF=nf.k, mu_NF=nf.mu, Pr_NF=nf.Pr,
        alpha_BF=nf.alpha_bf, alpha_p=nf.alpha_p, alpha_r=nf.alpha_r, n_eff=nf.n_eff,
        V=V, mdot=mass_flow(nf.rho, V), f=friction_factor(Re),
        Nu_BF=gnielinski(Re, nf.Pr),
        eta_eq4=eta_paper(nf.alpha_r, nf.phi, nf.n_eff),
        Nu_eq4=paper_eq4(Re, nf.Pr, nf.alpha_r, nf.phi, nf.n_eff),
        Nu_eq3=paper_eq3(Re, nf.Pr),
    )


if __name__ == "__main__":
    water = BaseFluid("water", 996.5, 4181.0, 0.613, 8.55e-4)
    cuo = Particle("CuO", 0.002, 6350.0, 535.6, 69.0, "spherical")
    for key, val in case_summary(water, [cuo], 10000).items():
        print(f"{key:10s} {val:.6g}")
