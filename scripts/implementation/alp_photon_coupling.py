"""Vector-photon contribution to the no-mixing ALP photon coefficient.

Natural units, masses and f_a in GeV, epsilon^0123=+1, metric (+---).
The coefficient convention is L = alpha_em/(4*pi) c_gamma a/f_a F Ftilde.
Couplings c_q^A and c_q^P are already in the no-mixing basis, in u,d,s
order. No meson mixing transformation or additional WZW term is applied.

The only script dependency is the explicitly reused alp_vv_width.py fit
reader. It reconstructs the neutral mass generators and signed one-trace
source residues, including z_V, from odd_preliminary_fit.json. The default
is its central_single_trace_negative provisional solution.

Invocation from the repository root::

    .venv/bin/python scripts/implementation/alp_photon_coupling.py \
        --axial 1 0 0 --pseudoscalar 0 0 0 \
        --alp-mass 1 --alp-decay-constant 1000 --c-gamma 0

This is the original-eleven-operator, tree-level, leading-electromagnetic
contribution with two real photons. Internal vectors are evaluated at
q_gamma^2=0, with no widths, loops, resummed mixing or extra odd operators.
The ALP mass must be in the domain in which this resonance model is used;
the implementation does not impose an additional phenomenological cutoff.
"""
import argparse
import json
import math

import numpy as np

if __package__:
    from .alp_vv_width import read_preliminary_fit_parameters, DEFAULT_FIT_LABEL
else:
    from alp_vv_width import read_preliminary_fit_parameters, DEFAULT_FIT_LABEL


def _finite_real(value, name):
    if isinstance(value, (bool, complex)) or not np.isscalar(value):
        raise ValueError(name + ' must be a finite real number')
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(name + ' must be a finite real number')
    return value


def _quark_couplings(values, name):
    couplings = np.asarray(values)
    if (couplings.shape != (3,) or np.iscomplexobj(couplings)
            or not np.all(np.isfinite(couplings))):
        raise ValueError(name + ' must contain three finite real entries (u,d,s)')
    return couplings.astype(float)


def alp_photon_coupling(axial_quark_couplings, pseudoscalar_quark_couplings,
                       alp_mass, alp_decay_constant, *, c_gamma=None,
                       fit_parameters=None):
    """Return dimensionless delta_c_gamma and its VJP/VVP contributions.

    An optional local ``c_gamma`` is added to delta_c_gamma to form
    ``c_gamma_eff``. It must already be in the same no-mixing basis.
    ``fit_parameters`` may be the result of the reused provisional fit
    reader; omitting it selects its documented default.

    The scalar G is defined by i M = i e^2 E G with
    E=epsilon(q1,eps1,q2,eps2). Thus delta_c_gamma=-4*pi^2*f_a*G.
    Its units are GeV^-1. The returned amplitude coefficient includes e^2;
    the coefficient of a F Ftilde is alpha_em*delta_c_gamma/(4*pi*f_a).

    Both photon placements in VJP and both vector assignments in VVP are
    included. No statistical uncertainty is attached to the provisional
    coefficient representative.
    """
    alp_mass = _finite_real(alp_mass, 'alp_mass')
    alp_decay_constant = _finite_real(alp_decay_constant, 'alp_decay_constant')
    if alp_mass < 0 or alp_decay_constant <= 0:
        raise ValueError('alp_mass must be nonnegative and alp_decay_constant positive')
    axial = _quark_couplings(axial_quark_couplings, 'axial_quark_couplings')
    pseudoscalar = _quark_couplings(
        pseudoscalar_quark_couplings, 'pseudoscalar_quark_couplings')
    fitted = fit_parameters if fit_parameters is not None else read_preliminary_fit_parameters()
    derivative_source = np.diag(2*axial/alp_decay_constant)
    pseudoscalar_matrix = np.diag(pseudoscalar)
    mass_source = -2/alp_decay_constant*(
        pseudoscalar_matrix@fitted.quark_mass_matrix
        + fitted.quark_mass_matrix@pseudoscalar_matrix)
    electric_charges = np.diag([2/3, -1/3, -1/3])
    odd = fitted.odd_coefficients
    neutral_vectors = ('rho0', 'omega', 'phi')
    real_photon_response = {
        name: fitted.photon_vector_residues[name]/fitted.vector_masses[name]**2
        for name in neutral_vectors}

    one_vector_contributions = {}
    two_vector_contributions = {}
    for name in neutral_vectors:
        generator_charge = fitted.vector_generators[name]@electric_charges
        derivative_trace = float(np.trace(generator_charge@derivative_source))
        mass_trace = float(np.trace(generator_charge@mass_source))
        ordered_vjp_numerator = 4/fitted.operator_mass*(
            alp_mass**2/2*(odd['c1']+odd['c2']-odd['c5'])*derivative_trace
            - 2*odd['c3']*mass_trace)
        # Each placement gives -2 sqrt(2) times the ordered numerator.
        one_vector_contributions[name] = (
            -4*math.sqrt(2)*real_photon_response[name]*ordered_vjp_numerator)
        for second_name in neutral_vectors:
            generator_product = (fitted.vector_generators[name]
                                 @fitted.vector_generators[second_name])
            derivative_pair_trace = float(np.trace(generator_product@derivative_source))
            mass_pair_trace = float(np.trace(generator_product@mass_source))
            vvp_numerator = (
                4*alp_mass**2*(odd['d1']-odd['d3'])*derivative_pair_trace
                - 16*odd['d2']*mass_pair_trace)
            # This ordered double sum includes both distinct vector assignments.
            two_vector_contributions[name+':'+second_name] = (
                2*real_photon_response[name]*real_photon_response[second_name]
                * vvp_numerator)
    vjp_scalar = math.fsum(one_vector_contributions.values())
    vvp_scalar = math.fsum(two_vector_contributions.values())
    coefficient_conversion = -4*math.pi**2*alp_decay_constant
    vjp_coefficient = coefficient_conversion*vjp_scalar
    vvp_coefficient = coefficient_conversion*vvp_scalar
    delta_c_gamma = vjp_coefficient+vvp_coefficient
    result = {
        'fit_label': fitted.label,
        'delta_c_gamma_vjp': vjp_coefficient,
        'delta_c_gamma_vvp': vvp_coefficient,
        'delta_c_gamma': delta_c_gamma,
        'resonance_scalar_GeV_inv': vjp_scalar+vvp_scalar,
        'resonance_amplitude_coefficient_GeV_inv': (
            4*math.pi*fitted.alpha_em*(vjp_scalar+vvp_scalar)),
        'resonance_a_FdualF_coefficient_GeV_inv': (
            fitted.alpha_em/(4*math.pi*alp_decay_constant)*delta_c_gamma),
    }
    if c_gamma is not None:
        c_gamma = _finite_real(c_gamma, 'c_gamma')
        result['c_gamma_local'] = c_gamma
        result['c_gamma_eff'] = c_gamma+delta_c_gamma
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--axial', type=float, nargs=3, required=True, metavar=('CU', 'CD', 'CS'))
    parser.add_argument('--pseudoscalar', type=float, nargs=3, required=True,
                        metavar=('CU', 'CD', 'CS'))
    parser.add_argument('--alp-mass', type=float, required=True, help='GeV')
    parser.add_argument('--alp-decay-constant', type=float, required=True, help='GeV')
    parser.add_argument('--c-gamma', type=float)
    parser.add_argument('--fit-path')
    parser.add_argument('--fit-label', default=DEFAULT_FIT_LABEL)
    arguments = parser.parse_args()
    fitted = read_preliminary_fit_parameters(
        arguments.fit_path, fit_label=arguments.fit_label)
    result = alp_photon_coupling(
        arguments.axial, arguments.pseudoscalar,
        arguments.alp_mass, arguments.alp_decay_constant,
        c_gamma=arguments.c_gamma, fit_parameters=fitted)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
