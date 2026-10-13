"""Flavor factors for P -> gamma* gamma with the off-shell photon singlet.

Convention: the pseudoscalar unit fields and eta mixing are exactly those in
derive_radiative_sm_flavor_algebra.py.  The t2 flavor trace is
tr(Q {Q, T_P}); no ALP spurion is present.

Run with: .venv/bin/python scripts/derive_radiative_tff_flavor.py
"""
import sympy as sp

from derive_radiative_sm_flavor_algebra import (
    ALONI_CT, ALONI_TH, E1, ETA8, PI0, Q, SQ2, TH, CT,
)


def tff_flavor(pseudoscalar: sp.Matrix) -> sp.Expr:
    """Flavor trace for the t2 two-photon structure."""
    return sp.simplify(sp.trace(Q * (Q * pseudoscalar + pseudoscalar * Q)))


def main() -> None:
    eta = CT * ETA8 + TH * E1
    eta_prime = -TH * ETA8 + CT * E1
    factors = {
        "pi0": tff_flavor(PI0),
        "eta": tff_flavor(eta),
        "eta_prime": tff_flavor(eta_prime),
    }
    expected = {
        "pi0": sp.sqrt(2) / 3,
        "eta": (CT + 4 * TH) / (3 * sp.sqrt(3)),
        "eta_prime": (-TH + 4 * CT) / (3 * sp.sqrt(3)),
    }
    for name in factors:
        assert sp.simplify(factors[name] - expected[name]) == 0
        print(f"{name}: {sp.simplify(factors[name])}")

    numerical = {name: sp.N(value.subs({CT: ALONI_CT, TH: ALONI_TH}), 12)
                 for name, value in factors.items()}
    assert numerical["pi0"] == sp.N(sp.sqrt(2) / 3, 12)
    print("Aloni mixing (sin(theta)=-1/3):")
    for name, value in numerical.items():
        print(f"  {name}: {value}")


if __name__ == "__main__":
    main()
