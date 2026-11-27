"""Leading WZW coefficient of a -> gamma P_1 P_2 with finite ALP-meson mixing.

Implements the boxed formula of
memory/alp/photon-two-pseudoscalar/finite-mixing/alp-lo-finite-mixing-wzw-pp-amplitude.md,
with the tadpole source and the mass operator of
memory/alp/photon/finite-mixing/alp-lo-finite-mixing-photon-coefficient.md, the pole
equations of memory/alp/meson-mixing/finite-mixing/spectrum-and-poles/pole-equations/alp-lo-exact-pole-equations.md and the
K0 doublet of memory/alp/meson-mixing/finite-mixing/spectrum-and-poles/pole-equations/alp-lo-k0-doublet-exact-pole.md. Neutral
kaons in the final state are the flavour states K0 (d sbar) and K0bar
(s dbar), and K_S = p K0 + q K0bar, K_L = p K0 - q K0bar with the
conjugated coefficients for final-state bras, as in
memory/alp/photon-two-pseudoscalar/odd-photon-neutral-kaon-selection-rules.md.
epsilon_K, p and q are those of three_meson_width.py.

Units: masses, f_a, f_pi, B_0 and the source matrices hat_S[0], hat_P[1] in
GeV, m_0_squared in GeV^2, e dimensionless. The returned coefficient is in
GeV^-3. Conventions: metric (+,-,-,-), epsilon^{0123} = +1, common WZW
orientation = -1 (already inserted), Q = diag(2/3, -1/3, -1/3) in (u, d, s)
order, N_c = 3, N_f = 3, U = exp(i Phi / f_pi), tr(lambda_a lambda_b) = 2
delta_ab, epsilon = f_pi/f_a. In Phi the field K0 sits at the (d, s) entry
with coefficient sqrt(2), so (lambda_6 phi_6 + lambda_7 phi_7) gives
K0 = (phi_6 - i phi_7)/sqrt(2) and K0bar = (phi_6 + i phi_7)/sqrt(2); the
outgoing generator lambda^dagger of a final K0 is sqrt(2) at the (s, d)
entry.

Physical pole masses are taken from ``calculate_masses`` of
calculate_alp_meson_masses_finite.py. Its labels 'K0_1', 'K0_2' (the two
real poles of the K0 doublet) are used only inside this module.

Minimal invocation from the repository root::

    .venv/bin/python scripts/implementation/alp_photon_pp_wzw_finite_mixing.py
"""
import math
import sys
from itertools import permutations
from pathlib import Path

import numpy as np

_IMPLEMENTATION_DIRECTORY = str(Path(__file__).resolve().parent)
if _IMPLEMENTATION_DIRECTORY not in sys.path:
    # calculate_alp_meson_masses_finite imports alda_neutral_inputs by bare name.
    sys.path.insert(0, _IMPLEMENTATION_DIRECTORY)

from calculate_alp_meson_masses_finite import SNAPSHOT, calculate_masses  # noqa: E402
from three_meson_width import EPSILON_K, _kaon_mixing_coefficients  # noqa: E402

NUMBER_OF_COLOURS = 3
NUMBER_OF_FLAVOURS = 3
QUARK_CHARGE_MATRIX = np.diag([2/3, -1/3, -1/3]).astype(complex)

_SQRT_TWO = math.sqrt(2.)
# Outgoing generators lambda^dagger, sqrt(2) times a unit matrix, as in
# alpaca_chiral.three_pseudoscalar._pseudoscalars.
CHARGED_OUTGOING_GENERATORS = {
    'pi+': _SQRT_TWO*np.array([[0, 0, 0], [1, 0, 0], [0, 0, 0]], dtype=complex),
    'pi-': _SQRT_TWO*np.array([[0, 1, 0], [0, 0, 0], [0, 0, 0]], dtype=complex),
    'K+': _SQRT_TWO*np.array([[0, 0, 0], [0, 0, 0], [1, 0, 0]], dtype=complex),
    'K-': _SQRT_TWO*np.array([[0, 0, 1], [0, 0, 0], [0, 0, 0]], dtype=complex),
}
FLAVOUR_KAON_OUTGOING_GENERATORS = {
    'K0': _SQRT_TWO*np.array([[0, 0, 0], [0, 0, 0], [0, 1, 0]], dtype=complex),
    'K0bar': _SQRT_TWO*np.array([[0, 0, 0], [0, 0, 1], [0, 0, 0]], dtype=complex),
}

# Neutral generators with tr(lambda_a lambda_b) = 2 delta_ab, in the order of
# calculate_masses: lambda_3, lambda_8, lambda_0, lambda_6, lambda_7.
NEUTRAL_GELL_MANN_GENERATORS = np.array([
    np.diag([1, -1, 0]),
    np.diag([1, 1, -2])/math.sqrt(3),
    np.eye(3)*math.sqrt(2/3),
    [[0, 0, 0], [0, 0, 1], [0, 1, 0]],
    [[0, 0, 0], [0, 0, -1j], [0, 1j, 0]],
], dtype=complex)

NEUTRAL_DAUGHTER_LABELS = ('pi0', 'eta', 'eta_prime', 'K0', 'K0bar', 'K_L', 'K_S')
DAUGHTER_LABELS = tuple(CHARGED_OUTGOING_GENERATORS) + NEUTRAL_DAUGHTER_LABELS
# Unmixed generators that can have a nonzero tr(U C). K0_2 is excluded: it is
# the decoupled combination of the K0 doublet, with tr(U_{K0_2} C) = 0.
_COUPLED_LABELS = ('pi0', 'eta', 'eta_prime', 'K0_1')
_KAON_POLE_LABELS = ('K0_1', 'K0_2')


def tadpole_source(hat_P_1, hat_Theta_1, B_0, m_0_squared):
    r"""Return $\mathcal C=B_0\widehat P_1-\frac{m_0^2}{N_f}\widehat\Theta_1\mathbf 1$ in GeV^2."""
    return (B_0*np.asarray(hat_P_1, dtype=complex)
            - m_0_squared/NUMBER_OF_FLAVOURS*hat_Theta_1*np.eye(3))


def unmixed_neutral_generators(quark_mass_term, m_0_squared, source):
    r"""Return the unmixed squared masses and generators of the neutral mesons.

    The generators $U_k$ are the eigen-matrices of the $s$-independent part of
    $\mathbb M_s(X)=\{B_0m_q,X\}+\frac{m_0^2}{N_f}\operatorname{tr}(X)\mathbf 1-sX$,
    $\mathbb M_0(U_k)=\mu_kU_k$, normalized to $\operatorname{tr}(U_kU_l)=2\delta_{kl}$:

    - 'pi0', 'eta', 'eta_prime': eigenvectors of the diagonal block in
      ascending order of $\mu_k$; the sign is fixed by a positive overlap with
      $\lambda_3$, $\lambda_8$, $\lambda_0$, respectively.
    - 'K0_1', 'K0_2' (internal labels of the two real poles of the K0
      doublet), both with $\mu=B_0(m_d+m_s)$. They are the coupled and
      decoupled combinations of Step 2 of the K0 doublet record, with $(d,s)$
      entries $(\mathcal C)_{ds}/|(\mathcal C)_{ds}|$ and
      $i(\mathcal C)_{ds}/|(\mathcal C)_{ds}|$:
      'K0_1' $=(\operatorname{tr}(\lambda_6\mathcal C)\lambda_6+\operatorname{tr}(\lambda_7\mathcal C)\lambda_7)/\sqrt{\operatorname{tr}(\lambda_6\mathcal C)^2+\operatorname{tr}(\lambda_7\mathcal C)^2}$
      and 'K0_2' $=(\operatorname{tr}(\lambda_7\mathcal C)\lambda_6-\operatorname{tr}(\lambda_6\mathcal C)\lambda_7)/\sqrt{\operatorname{tr}(\lambda_6\mathcal C)^2+\operatorname{tr}(\lambda_7\mathcal C)^2}$,
      which is minus the normalized generator
      $-\operatorname{tr}(\lambda_7\mathcal C)\lambda_6+\operatorname{tr}(\lambda_6\mathcal C)\lambda_7$
      of the record (proportional to it, as the record states). When
      $(\mathcal C)_{ds}=0$ both combinations decouple (record, Step 2) and
      $\lambda_6$, $\lambda_7$ are used, as in ``calculate_masses``. The
      flavour and $K_{L,S}$ legs do not depend on these signs.

    ``quark_mass_term`` is $B_0m_q$ (diagonal, GeV^2), ``source`` is
    $\mathcal C$. Returns two dicts keyed by label: squared masses (GeV^2)
    and 3x3 generators.
    """
    diagonal = NEUTRAL_GELL_MANN_GENERATORS[:3]
    block = np.array([[np.real(np.trace(first @ (quark_mass_term @ second + second @ quark_mass_term))/2
                               + m_0_squared/NUMBER_OF_FLAVOURS*np.trace(first).real*np.trace(second).real/2)
                       for second in diagonal] for first in diagonal])
    squared_masses, vectors = np.linalg.eigh(block)
    unmixed_squared_masses, generators = {}, {}
    for column, label in enumerate(('pi0', 'eta', 'eta_prime')):
        generator = np.tensordot(vectors[:, column], diagonal, axes=1)
        if np.trace(generator @ diagonal[column]).real < 0:
            generator = -generator
        unmixed_squared_masses[label] = squared_masses[column]
        generators[label] = generator
    trace_lambda_6, trace_lambda_7 = np.real([np.trace(NEUTRAL_GELL_MANN_GENERATORS[3] @ source),
                                              np.trace(NEUTRAL_GELL_MANN_GENERATORS[4] @ source)])
    kaon_norm = math.hypot(trace_lambda_6, trace_lambda_7)
    first_kaon, second_kaon = (trace_lambda_6/kaon_norm, trace_lambda_7/kaon_norm) if kaon_norm else (1., 0.)
    generators['K0_1'] = first_kaon*NEUTRAL_GELL_MANN_GENERATORS[3] + second_kaon*NEUTRAL_GELL_MANN_GENERATORS[4]
    generators['K0_2'] = second_kaon*NEUTRAL_GELL_MANN_GENERATORS[3] - first_kaon*NEUTRAL_GELL_MANN_GENERATORS[4]
    kaon_squared_mass = np.real(quark_mass_term[1, 1] + quark_mass_term[2, 2])
    unmixed_squared_masses['K0_1'] = unmixed_squared_masses['K0_2'] = kaon_squared_mass
    return unmixed_squared_masses, generators


def alp_leg(alp_squared_mass, mixing_ratio, unmixed_squared_masses, generators, source):
    r"""Return the ALP leg $2\epsilon X/\sqrt{1+2\epsilon^2\operatorname{tr}X^2}$, $X=\mathbb M_{m_a^2}^{-1}(\mathcal C)$.

    $X=\sum_k U_k\operatorname{tr}(U_k\mathcal C)/(2(\mu_k-m_a^2))$ is the
    eigenbasis form of $\mathbb M_s^{-1}$ of the exact pole-equation record.
    No ALP-meson degeneracy is assumed, so no denominator vanishes.
    """
    weights = {label: np.trace(generators[label] @ source).real
               / (unmixed_squared_masses[label] - alp_squared_mass) for label in _COUPLED_LABELS}
    matrix = sum(weights[label]*generators[label] for label in _COUPLED_LABELS)
    return mixing_ratio*matrix/math.sqrt(1 + mixing_ratio**2*sum(weight**2 for weight in weights.values()))


def neutral_pole_leg(label, pole_squared_mass, alp_squared_mass, mixing_ratio,
                     unmixed_squared_masses, generators, source):
    r"""Return the leg of a real neutral meson pole, with the sign of the boxed record expression.

    ``label`` is one of 'pi0', 'eta', 'eta_prime', 'K0_1', 'K0_2'. For a
    pole $j$ at $s_j=m_j^2$ the record's leg is
    $2\epsilon X/\sqrt{1+2\epsilon^2\operatorname{tr}X^2}$ with
    $X=\mathbb M_{s_j}^{-1}(\mathcal C)=\sum_kU_k\operatorname{tr}(U_k\mathcal C)/(2(\mu_k-s_j))$;
    its ALP component $1/\sqrt{1+2\epsilon^2\operatorname{tr}X^2}$ is
    positive. The term $k=j$ has the denominator $\mu_j-s_j$, which the exact
    meson pole equation (Step 5 of the exact pole-equation record),
    $(\mu_j-s_j)D_j=\epsilon^2\operatorname{tr}(U_j\mathcal C)^2$ with
    $D_j=M^2_{\hat a\hat a,\rm can}-s_j-\epsilon^2\sum_{k\ne j}\operatorname{tr}(U_k\mathcal C)^2/(\mu_k-s_j)$
    and $M^2_{\hat a\hat a,\rm can}=m_a^2+\epsilon^2\sum_k\operatorname{tr}(U_k\mathcal C)^2/(\mu_k-m_a^2)$
    (Step 3), replaces. Multiplying numerator and normalization by
    $|\epsilon\operatorname{tr}(U_j\mathcal C)|$ gives the same leg as

    $$
    \operatorname{sign}(\operatorname{tr}(U_j\mathcal C))\,
    \frac{D_jU_j+\epsilon^2\operatorname{tr}(U_j\mathcal C)\sum_{k\ne j}\frac{\operatorname{tr}(U_k\mathcal C)}{\mu_k-s_j}U_k}
    {\sqrt{D_j^2+\epsilon^2\operatorname{tr}(U_j\mathcal C)^2
    +\sum_{k\ne j}\Big(\frac{\epsilon^2\operatorname{tr}(U_j\mathcal C)\operatorname{tr}(U_k\mathcal C)}{\mu_k-s_j}\Big)^2}},
    $$

    which contains no difference $\mu_j-s_j$ and keeps full precision for
    any $f_a$. For $\epsilon\to0$ it tends to
    $\operatorname{sign}(\operatorname{tr}(U_j\mathcal C)D_j)\,U_j$.
    For $\operatorname{tr}(U_j\mathcal C)=0$ the record's expression is
    singular and the record gives the leg as the pole's own unit generator,
    without fixing its sign; the expression above then gives
    $\operatorname{sign}(D_j)\,U_j$ (the sign function is taken as $+1$ at
    zero). This happens for 'K0_1' when $(\mathcal C)_{ds}=0$.

    'K0_2' is the decoupled pole of the K0 doublet record (Step 2): its pole
    is $s=\mu_K$ exactly, with zero ALP component, for any sources. The
    expression above does not apply to it ($D_j$ contains
    $1/(\mu_{K0_1}-s_j)=1/0$), and its leg is the record's finite
    expression, the unit generator $U_{K0_2}$ of ``unmixed_neutral_generators``.
    """
    if label == 'K0_2':
        return generators['K0_2']
    couplings = {other: np.trace(generators[other] @ source).real for other in _COUPLED_LABELS}
    own_coupling = couplings[label]
    others = [other for other in _COUPLED_LABELS if other != label]
    alp_pole_sum = sum(couplings[other]**2/(unmixed_squared_masses[other] - alp_squared_mass)
                       for other in _COUPLED_LABELS)
    pole_sum = sum(couplings[other]**2/(unmixed_squared_masses[other] - pole_squared_mass) for other in others)
    reduced_alp_entry = alp_squared_mass - pole_squared_mass + mixing_ratio**2*(alp_pole_sum - pole_sum)
    weights = {other: mixing_ratio**2*own_coupling*couplings[other]
               / (unmixed_squared_masses[other] - pole_squared_mass) for other in others}
    matrix = reduced_alp_entry*generators[label] + sum(weights[other]*generators[other] for other in others)
    normalization = math.sqrt(reduced_alp_entry**2 + (mixing_ratio*own_coupling)**2
                              + sum(weight**2 for weight in weights.values()))
    sign = -1. if own_coupling < 0 else 1.
    return sign*matrix/normalization


def flavour_kaon_leg(outgoing_generator, kaon_pole_legs):
    r"""Return the leg of a final K0 or K0bar from the legs of the two real kaon poles.

    ``outgoing_generator`` is $\lambda^\dagger$ of the flavour state ($\sqrt2$
    at the $(s,d)$ entry for K0, at $(d,s)$ for K0bar). ``kaon_pole_legs`` are
    the legs $L_1$, $L_2$ of the coupled and decoupled poles ('K0_1', 'K0_2').

    Definition (the flavour state at finite mixing): with $|1\rangle$,
    $|2\rangle$ the two kaon pole states and $d_k$ the $(\lambda_6,\lambda_7)$
    part of $L_k$, $d_k=\sum_{n=6,7}\lambda_n\operatorname{tr}(\lambda_nL_k)/2$,
    the flavour state is
    $|K^0\rangle=\sum_k|k\rangle\langle k|K^0\rangle_0/\sqrt{\operatorname{tr}(d_k^2)/2}$,
    where $\langle k|K^0\rangle_0$ is the overlap with the unmixed flavour
    state, $\langle K^0|k\rangle_0=\operatorname{tr}(\lambda^\dagger d_k)/2$.
    By the K0 doublet record, $d_1$ is along the coupled generator and
    $d_2=L_2$ is the decoupled one, so the coefficients are the
    decomposition of the unmixed $K^0$ onto the two doublet generators,
    multiplied by the sign of each pole state. The outgoing leg is then

    $$
    L_{K^0}=\sum_{k=1,2}\frac{\operatorname{tr}(\lambda^\dagger d_k)/2}{\sqrt{\operatorname{tr}(d_k^2)/2}}\,L_k ,
    $$

    which does not change under $L_k\to-L_k$, tends to $\lambda^\dagger$ for
    $\epsilon\to0$, and keeps $\langle\bar K^0|K^0\rangle=0$ and unit norm.
    """
    leg = np.zeros((3, 3), dtype=complex)
    for pole_leg in kaon_pole_legs:
        doublet_part = sum(np.trace(generator @ pole_leg)/2*generator
                           for generator in NEUTRAL_GELL_MANN_GENERATORS[3:])
        overlap = np.trace(outgoing_generator @ doublet_part)/2
        leg += overlap/math.sqrt(np.trace(doublet_part @ doublet_part).real/2)*pole_leg
    return leg


def physical_kaon_legs(kaon_0_leg, kaon_0_bar_leg, epsilon_K=EPSILON_K):
    r"""Return the outgoing legs of $K_S$ and $K_L$ from the $K^0$ and $\bar K^0$ legs.

    $|K_S\rangle=p|K^0\rangle+q|\bar K^0\rangle$,
    $|K_L\rangle=p|K^0\rangle-q|\bar K^0\rangle$, with
    $p=(1+\epsilon_K)/\sqrt{2(1+|\epsilon_K|^2)}$,
    $q=(1-\epsilon_K)/\sqrt{2(1+|\epsilon_K|^2)}$ from three_meson_width.py.
    For final-state bras the coefficients are complex-conjugated (odd-photon
    selection-rule record): $L_{K_S}=p^*L_{K^0}+q^*L_{\bar K^0}$,
    $L_{K_L}=p^*L_{K^0}-q^*L_{\bar K^0}$. Returns ``{'K_S': ..., 'K_L': ...}``.
    """
    p_coefficient, q_coefficient = _kaon_mixing_coefficients(epsilon_K)
    return {'K_S': np.conj(p_coefficient)*kaon_0_leg + np.conj(q_coefficient)*kaon_0_bar_leg,
            'K_L': np.conj(p_coefficient)*kaon_0_leg - np.conj(q_coefficient)*kaon_0_bar_leg}


def alternating_trace(first, second, third):
    r"""Return $\sum_{\sigma\in S_3}\operatorname{sgn}(\sigma)\operatorname{tr}(Q X_{\sigma(1)}X_{\sigma(2)}X_{\sigma(3)})$."""
    legs = (first, second, third)
    result = 0j
    for permutation in permutations(range(3)):
        inversions = sum(permutation[i] > permutation[j] for i in range(3) for j in range(i + 1, 3))
        result += (-1)**inversions*np.trace(
            QUARK_CHARGE_MATRIX @ legs[permutation[0]] @ legs[permutation[1]] @ legs[permutation[2]])
    return complex(result)


def _all_legs(m_a0, f_a, hat_S, hat_P, hat_Theta, f_pi, B_0, m_0_squared, epsilon_K):
    """Return the legs of the ALP and of every daughter label, the internal kaon pole legs, and the pole masses."""
    scalar_sources = list(hat_S) + [np.zeros((3, 3), dtype=complex)]*max(0, 3 - len(hat_S))
    pseudoscalar_sources = list(hat_P) + [np.zeros((3, 3), dtype=complex)]*max(0, 2 - len(hat_P))
    gluonic_sources = list(hat_Theta) + [0.]*max(0, 2 - len(hat_Theta))
    pole_masses = calculate_masses(
        f_a, m_a0, list(scalar_sources), list(pseudoscalar_sources), list(gluonic_sources),
        B_0=B_0, f_pi=f_pi, m_0_squared=m_0_squared)
    mixing_ratio = f_pi/f_a
    quark_mass_term = B_0*np.asarray(scalar_sources[0], dtype=complex)
    source = tadpole_source(pseudoscalar_sources[1], gluonic_sources[1], B_0, m_0_squared)
    unmixed_squared_masses, generators = unmixed_neutral_generators(quark_mass_term, m_0_squared, source)
    alp_squared_mass = pole_masses['a']**2
    pole_legs = {label: neutral_pole_leg(label, pole_masses[label]**2, alp_squared_mass, mixing_ratio,
                                         unmixed_squared_masses, generators, source)
                 for label in ('pi0', 'eta', 'eta_prime') + _KAON_POLE_LABELS}
    legs = {'a': alp_leg(alp_squared_mass, mixing_ratio, unmixed_squared_masses, generators, source),
            **CHARGED_OUTGOING_GENERATORS,
            **{label: pole_legs[label] for label in ('pi0', 'eta', 'eta_prime')}}
    for label, outgoing in FLAVOUR_KAON_OUTGOING_GENERATORS.items():
        legs[label] = flavour_kaon_leg(outgoing, [pole_legs[pole] for pole in _KAON_POLE_LABELS])
    legs.update(physical_kaon_legs(legs['K0'], legs['K0bar'], epsilon_K))
    return legs, pole_legs, pole_masses


def external_legs(meson_1, meson_2, m_a0, f_a, hat_S, hat_P, hat_Theta, *, f_pi, B_0, m_0_squared,
                  epsilon_K=EPSILON_K):
    r"""Return the legs (ALP, meson_1, meson_2) and the physical masses.

    Inputs as in ``wzw_coefficient``. Charged daughters have the leg
    $\lambda^\dagger$ of ``CHARGED_OUTGOING_GENERATORS``; 'pi0', 'eta',
    'eta_prime' the pole leg of ``neutral_pole_leg``; 'K0', 'K0bar' the leg
    of ``flavour_kaon_leg``; 'K_S', 'K_L' that of ``physical_kaon_legs``.
    The masses (GeV) are those of ``calculate_masses`` for 'a', 'pi0',
    'eta', 'eta_prime' and the charged mesons. The neutral kaon flavour
    states and $K_{L,S}$ are combinations of two poles of different mass at
    finite mixing, so no mass is returned for them.
    """
    for label in (meson_1, meson_2):
        if label not in DAUGHTER_LABELS:
            raise ValueError(f'unknown daughter {label!r}; use one of {DAUGHTER_LABELS} '
                             '(the internal kaon pole labels K0_1, K0_2 are not accepted)')
    legs, _, pole_masses = _all_legs(m_a0, f_a, hat_S, hat_P, hat_Theta, f_pi, B_0, m_0_squared, epsilon_K)
    masses = {label: mass for label, mass in pole_masses.items() if label not in _KAON_POLE_LABELS}
    return [legs['a'], legs[meson_1], legs[meson_2]], masses


def wzw_coefficient(meson_1, meson_2, m_a0, f_a, hat_S, hat_P, hat_Theta, *,
                    f_pi, B_0, m_0_squared, e, epsilon_K=EPSILON_K):
    r"""Return the LO WZW coefficient $\mathcal A$ of a -> gamma P_1 P_2 at finite mixing, in GeV^-3.

    $\mathcal M=\epsilon_{\mu\nu\rho\sigma}\varepsilon_\gamma^\mu p_1^\nu p_2^\rho k^\sigma\,\mathcal A$ with

    $$
    \mathcal A=\frac{N_ce}{24\pi^2f_\pi^3}\sum_{\sigma\in S_3}\operatorname{sgn}(\sigma)
    \operatorname{tr}\big(Q\,L_{\sigma(1)}L_{\sigma(2)}L_{\sigma(3)}\big),
    $$

    $(L_1,L_2,L_3)$ the legs of (ALP, meson_1, meson_2). A neutral pole leg
    at physical mass $m$ is
    $2\epsilon\,\mathbb M_{m^2}^{-1}(\mathcal C)/\sqrt{1+2\epsilon^2\operatorname{tr}[(\mathbb M_{m^2}^{-1}(\mathcal C))^2]}$
    with the sign of the boxed record (see ``alp_leg`` and
    ``neutral_pole_leg``); a charged leg is its outgoing generator
    $\lambda^\dagger$; K0 and K0bar legs are combinations of the two real
    kaon pole legs (``flavour_kaon_leg``), and K_S, K_L legs are
    $p^*L_{K^0}\pm q^*L_{\bar K^0}$ (``physical_kaon_legs``). For three
    neutral pole legs the prefactor times $(2\epsilon)^3$ is
    $N_ce/(3\pi^2f_a^3)$, the boxed form.
    $\mathcal C=B_0\widehat P_1-\frac{m_0^2}{N_f}\widehat\Theta_1\mathbf 1$
    and $\mathbb M_s(X)=\{B_0m_q,X\}+\frac{m_0^2}{N_f}\operatorname{tr}(X)\mathbf 1-sX$,
    with $m_q=$ ``hat_S[0]``, $\widehat P_1=$ ``hat_P[1]``,
    $\widehat\Theta_1=$ ``hat_Theta[1]``. This $\mathcal C$ is the source
    whose $-\epsilon\operatorname{tr}(\lambda_i\mathcal C)$ is the ALP-meson
    mass entry of ``calculate_masses``.

    ``hat_S``, ``hat_P``, ``hat_Theta`` are lists of derivatives at
    $a/f_a=0$ indexed by order, as in ``calculate_masses``; missing entries
    are zero. hat_S[2] and $m_{a,0}$ = ``m_a0`` enter only through the
    physical masses. ``epsilon_K`` (complex) defaults to ``EPSILON_K`` of
    three_meson_width.py and enters only the K_S and K_L legs.

    Daughter labels: 'pi+', 'pi-', 'K+', 'K-', 'pi0', 'eta', 'eta_prime',
    'K0', 'K0bar', 'K_L', 'K_S'. Other labels, including the internal kaon
    pole labels 'K0_1', 'K0_2', raise ValueError.

    LO, tree level, leading WZW term, no ALP-meson degeneracy. No kinematic
    threshold is applied.
    """
    legs, _ = external_legs(meson_1, meson_2, m_a0, f_a, hat_S, hat_P, hat_Theta,
                            f_pi=f_pi, B_0=B_0, m_0_squared=m_0_squared, epsilon_K=epsilon_K)
    return NUMBER_OF_COLOURS*e/(24*math.pi**2*f_pi**3)*alternating_trace(*legs)


if __name__ == '__main__':
    import json
    snapshot = json.loads(SNAPSHOT.read_text())
    example_B_0 = 2.5
    quark_masses = np.diag(snapshot['B0mq_gev2'])/example_B_0
    constants = dict(f_pi=snapshot['fpi_gev'], B_0=example_B_0, m_0_squared=snapshot['m0_squared_gev2'],
                     e=math.sqrt(4*math.pi/137.035999))
    m_a0, hat_Theta_1 = 1.5, .4
    print('Inputs: m_a0 =', m_a0, 'GeV, B_0 =', example_B_0, 'GeV, f_pi =', constants['f_pi'],
          'GeV, m_0_squared =', constants['m_0_squared'], 'GeV^2, e =', constants['e'],
          ', epsilon_K =', EPSILON_K)
    print('hat_S = [diag', list(np.diag(quark_masses)), 'GeV, 0, 0]; hat_Theta = [0,', hat_Theta_1, ']')
    for name, axial_minus_pseudoscalar in (
            ('diagonal', np.diag([1., .5, -.3]).astype(complex)),
            ('flavour-changing', np.array([[1., 0., 0.], [0., .5, .2+.3j], [0., .2-.3j, -.3]]))):
        hat_P_1 = -(quark_masses @ axial_minus_pseudoscalar + axial_minus_pseudoscalar @ quark_masses)
        print(f'\nhat_P = [0, -{{hat_S[0], c^AP}}] with {name} c^AP =\n', axial_minus_pseudoscalar)
        for f_a in (10., 1e4, 1e8):
            scale = NUMBER_OF_COLOURS*constants['e']/(12*math.pi**2*f_a*constants['f_pi']**2)
            print(f'  f_a = {f_a:g} GeV, N_c e/(12 pi^2 f_a f_pi^2) = {scale:.6e} GeV^-3')
            for pair in (('pi+', 'pi-'), ('K0', 'pi0'), ('K0bar', 'pi0'), ('K0', 'K0bar'), ('K_L', 'K_S'),
                         ('K_S', 'pi0'), ('K_L', 'pi0'), ('K_S', 'K_S'), ('pi0', 'eta')):
                value = wzw_coefficient(*pair, m_a0, f_a, [quark_masses], [0*quark_masses, hat_P_1],
                                        [0., hat_Theta_1], **constants)
                print(f'    {pair}: A = {value:.6e} GeV^-3, A/scale = {value/scale:.6e}')
