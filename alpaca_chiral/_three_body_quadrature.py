"""Converged Gauss-Legendre integration of three-body Dalitz domains."""
import numpy as np
from scipy.integrate import quad
from scipy.special import roots_legendre


def _pole_boundaries(lower, upper, poles):
    """Return physical endpoints and finite-width pole subdivisions in GeV²."""
    points = {lower, upper}
    for pole in poles:
        for offset in (-8., -2., 0., 2., 8.):
            point = pole.real+offset*abs(pole.imag)
            if lower < point < upper:
                points.add(point)
    return np.array(sorted(points))


def integrate_three_body(fixed_squared, limits, lower, upper, total_squared,
                         invariant_poles, *, order=16, epsrel=1e-5):
    """Integrate a coherent squared amplitude over its two Dalitz invariants.

    ``fixed_squared(s)`` returns an array-compatible function of the inner
    invariant. ``limits(s)`` gives its physical endpoints; pole groups denote
    s,t,u with u=total_squared-s-t. All squared invariants use GeV². Four
    orders are compared on a sine-squared outer map and pole-separated inner
    intervals. Unresolved outer grids use adaptive outer integration with
    converged inner rules and an adaptive inner fallback.
    No physical propagator approximation or pole regularization is made.
    """
    rules = [roots_legendre(order*factor) for factor in (1, 2, 4, 8)]
    polynomial_nodes, polynomial_weights = roots_legendre(3)
    s_poles, t_poles, u_poles = invariant_poles

    def integrate_inner(s):
        """Integrate all pole-separated inner intervals at fixed outer s."""
        inner_lower, inner_upper = limits(s)
        if inner_upper <= inner_lower:
            return 0.
        squared = fixed_squared(s)
        if not t_poles and not u_poles:
            # Both present squared amplitudes are at most quadratic in t
            # without crossed propagators. Three nodes integrate them exactly.
            half_range = (inner_upper-inner_lower)/2
            coordinates = (inner_upper+inner_lower)/2+half_range*polynomial_nodes
            return float(half_range*(squared(coordinates)@polynomial_weights))
        poles = (*t_poles, *(total_squared-s-pole.conjugate() for pole in u_poles))
        boundaries = _pole_boundaries(inner_lower, inner_upper, poles)
        centers = (boundaries[1:]+boundaries[:-1])/2
        half_ranges = np.diff(boundaries)/2
        previous = None
        for nodes, weights in rules:
            coordinates = centers[:, None]+half_ranges[:, None]*nodes
            integral = float(np.sum(half_ranges*(squared(coordinates)@weights)))
            if previous is not None and abs(integral-previous) <= epsrel/4*abs(integral):
                return integral
            previous = integral
        return quad(lambda t: float(squared(t)), inner_lower, inner_upper,
                    points=boundaries[1:-1], epsabs=0., epsrel=epsrel/4,
                    limit=max(200, len(boundaries)*4))[0]

    outer_boundaries = _pole_boundaries(lower, upper, s_poles)
    previous = None
    for nodes, weights in rules:
        angles = (nodes+1)*np.pi/4
        intervals = np.diff(outer_boundaries)
        s_values = (outer_boundaries[:-1, None]+intervals[:, None]*np.sin(angles)**2).ravel()
        outer_weights = (intervals[:, None]*weights*np.pi/4*np.sin(2*angles)).ravel()
        inner_lower, inner_upper = limits(s_values)
        boundaries = [inner_lower, inner_upper]
        for reflected, poles in ((False, t_poles), (True, u_poles)):
            for pole in poles:
                center = total_squared-s_values-pole.real if reflected else np.full_like(s_values, pole.real)
                for offset in (-8., -2., 0., 2., 8.):
                    boundaries.append(np.clip(center+offset*abs(pole.imag), inner_lower, inner_upper))
        boundaries = np.sort(np.array(boundaries).T, axis=1)
        half_ranges = np.diff(boundaries, axis=1)/2
        inner_nodes, inner_weights = ((nodes, weights) if t_poles or u_poles
                                      else (polynomial_nodes, polynomial_weights))
        coordinates = ((boundaries[:, 1:]+boundaries[:, :-1])[:, :, None]/2
                       +half_ranges[:, :, None]*inner_nodes)
        squared = fixed_squared(s_values[:, None, None])(coordinates)
        inner_integrals = np.sum(half_ranges*(squared@inner_weights), axis=1)
        integral = float(outer_weights@inner_integrals)
        if previous is not None and abs(integral-previous) <= epsrel/4*abs(integral):
            return integral
        previous = integral

    points = outer_boundaries[1:-1]
    return quad(integrate_inner, lower, upper, points=points,
                epsabs=0., epsrel=epsrel, limit=max(200, len(points)*4))[0]
