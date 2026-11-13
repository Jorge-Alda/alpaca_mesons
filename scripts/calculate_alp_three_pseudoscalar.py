"""Tree-level LO a -> three pseudoscalars with finite ALP-meson mixing.

Sources are derivatives at x=a/f_a=0. Masses, f_a, f_pi, B_0 and mass
sources use GeV, m_0_squared and Mandelstam invariants GeV^2, widths GeV.
Require S[0:5], P[0:4], Theta[0:4]; supply zero matrices/scalars for
vanishing derivatives. Aligned vacuum and charge-conserving sources assumed.
Labels are those of calculate_alp_meson_masses, with K0_1, K0_2 real pole
states. The labels of mixed poles use the same maximum-overlap convention.
An accessible internal pole with nonzero residue requires physics beyond the
stable-particle tree calculation and raises ValueError.

Example:
    PYTHONPATH=scripts .venv/bin/python -c 'import numpy as np; from calculate_alp_three_pseudoscalar import decay_width; z=np.zeros((3,3)); print(decay_width(1., .6, [np.diag([.003,.005,.1]),z,z,z,z], [z,z,z,z], [0.,0.,0.,0.], ("pi0","pi+","pi-"), f_pi=.092, B_0=2.5, m_0_squared=.75))'
"""
from collections import Counter
from itertools import combinations, permutations
from math import factorial
import numpy as np
from alpaca_chiral._three_body_quadrature import integrate_three_body
from scipy.optimize import linear_sum_assignment
from calculate_alp_meson_masses import calculate_masses, _matrix


def _anticommutator(a, b):
    return a @ b + b @ a


def _commutator(a, b):
    return a @ b - b @ a


def _basis():
    unit = np.eye(3, dtype=complex)
    def pair(i, j):
        m = np.zeros((3, 3), complex); m[i, j] = 1
        return m
    u, v = pair(0, 1), pair(0, 2)
    d = pair(1, 2)
    neutral = [np.diag([1, -1, 0]), np.diag([1, 1, -2])/np.sqrt(3),
               unit*np.sqrt(2/3), d+d.T, -1j*(d-d.T)]
    return np.array(neutral + [u+u.T, -1j*(u-u.T), v+v.T, -1j*(v-v.T)]), {
        'pi+': np.sqrt(2)*u.T, 'pi-': np.sqrt(2)*u,
        'K+': np.sqrt(2)*v.T, 'K-': np.sqrt(2)*v}


def prepare_amplitude(f_a, m_a0, hat_S, hat_P, hat_Theta, mesons, *, f_pi, B_0, m_0_squared):
    """Return a prepared amplitude object; call it with s12, s23 in GeV^2.

    The output is a complex amplitude in the convention S=1+iT. A separate
    decay_width function integrates its modulus squared over physical phase space.
    """
    if len(hat_S) < 5 or len(hat_P) < 4 or len(hat_Theta) < 4:
        raise ValueError('finite mixing requires S[0:5], P[0:4], Theta[0:4]')
    if len(mesons) != 3:
        raise ValueError('provide three meson labels')
    masses = calculate_masses(f_a, m_a0, hat_S, hat_P, hat_Theta,
                              f_pi=f_pi, B_0=B_0, m_0_squared=m_0_squared)
    for name in mesons:
        if name not in masses or name == 'a':
            raise ValueError(f'unknown meson label: {name}')
    scalar = [_matrix(x, f'hat_S[{i}]') for i, x in enumerate(hat_S[:5])]
    pseudo = [_matrix(x, f'hat_P[{i}]') for i, x in enumerate(hat_P[:4])]
    theta_input = np.asarray(hat_Theta[:4])
    if not np.all(np.isreal(theta_input)):
        raise ValueError('Theta derivatives must be finite and real')
    theta = np.asarray(theta_input.real, dtype=float)
    if not np.all(np.isfinite(theta)):
        raise ValueError('Theta derivatives must be finite and real')
    for matrix in scalar + pseudo:
        if not np.allclose(matrix[0, 1:], 0, atol=1e-14, rtol=0):
            raise ValueError('source matrices must conserve electric charge')
    generators, charged = _basis()
    ratio = f_pi/f_a
    chi = 2*B_0*scalar[0]
    matrix = np.zeros((10, 10))
    for i, first in enumerate(generators):
        for j, second in enumerate(generators):
            matrix[i,j] = (np.trace(_anticommutator(first,second)@chi)/4 +
                           m_0_squared*np.trace(first)*np.trace(second)/6).real
        matrix[i,9] = matrix[9,i] = (-ratio*B_0*np.trace(first@pseudo[1]) +
                                       ratio*m_0_squared*theta[1]*np.trace(first)/3).real
    matrix[9,9] = m_a0**2 + ratio**2*(2*m_0_squared*theta[1]**2/3-B_0*np.trace(scalar[2]).real)
    neutral_indices = [0,1,2,3,4,9]
    neutral_matrix = matrix[np.ix_(neutral_indices, neutral_indices)]
    _, diagonal_vectors = np.linalg.eigh(neutral_matrix[:3,:3])
    reference = np.zeros((6,6)); reference[:3,:3] = diagonal_vectors
    coupling = neutral_matrix[3:5,5]; norm = np.linalg.norm(coupling)
    if norm:
        reference[3:5,3] = coupling/norm
        reference[3:5,4] = [-reference[4,3],reference[3,3]]
    else:
        reference[3:5,3:5] = np.eye(2)
    reference[5,5] = 1
    poles, vectors = np.linalg.eigh(neutral_matrix)
    rows, cols = linear_sum_assignment(-np.abs(reference.T@vectors)**2)
    labels = ('pi0','eta','eta_prime','K0_1','K0_2','a')
    columns = {labels[r]: vectors[:,c] for r,c in zip(rows,cols)}
    def direction(name):
        if name in charged:
            return (charged[name], 0.)
        column = columns[name]
        return (sum((column[i]*generators[i] for i in range(5)), np.zeros((3,3),complex)), column[5])
    external = [direction('a')] + [direction(name) for name in mesons]
    internal = [(generator,0.) for generator in generators] + [(np.zeros((3,3),complex),1.)]
    def z1(leg):
        return np.trace(leg[0])/f_pi + 2*theta[1]*leg[1]/f_a
    def z2(a,b):
        return 2*theta[2]*a[1]*b[1]/f_a**2
    def z3(a,b,c):
        return 2*theta[3]*a[1]*b[1]*c[1]/f_a**3
    plus = [4*B_0*scalar[n]/f_a**n for n in range(5)]
    minus = [4j*B_0*pseudo[n]/f_a**n for n in range(4)]
    def cubic(a,b,c):
        x,y,z = a,b,c
        meson = x[1]*_anticommutator(y[0],z[0])+y[1]*_anticommutator(z[0],x[0])+z[1]*_anticommutator(x[0],y[0])
        two = x[1]*y[1]*z[0]+x[1]*z[1]*y[0]+y[1]*z[1]*x[0]
        value = -np.trace(plus[1]@meson)/8 + f_pi*np.trace(minus[2]@two)/(4j) + f_pi**2*x[1]*y[1]*z[1]*np.trace(plus[3])/4
        return value - m_0_squared*f_pi**2/6*(z1(x)*z2(y,z)+z1(y)*z2(x,z)+z1(z)*z2(x,y))
    matrices = [x[0] for x in external]; overlaps = [x[1] for x in external]
    contact = sum(np.trace(plus[0]@matrices[i]@matrices[j]@matrices[k]@matrices[l]) for i,j,k,l in permutations(range(4)))/(96*f_pi**2)
    for a in range(4):
        other = [i for i in range(4) if i != a]
        contact -= overlaps[a]*sum(np.trace(minus[1]@matrices[i]@matrices[j]@matrices[k]) for i,j,k in permutations(other))/(24j*f_pi)
    for a,b in combinations(range(4),2):
        i,j = [x for x in range(4) if x not in (a,b)]
        contact -= overlaps[a]*overlaps[b]*np.trace(plus[2]@_anticommutator(matrices[i],matrices[j]))/8
    for i in range(4):
        contact += f_pi*np.prod([overlaps[j] for j in range(4) if j != i])*np.trace(minus[3]@matrices[i])/(4j)
    contact += f_pi**2*np.prod(overlaps)*np.trace(plus[4])/4
    contact -= m_0_squared*f_pi**2/6*(sum(z1(external[i])*z3(*[external[j] for j in range(4) if j != i]) for i in range(4)) + sum(z2(external[a],external[b])*z2(external[c],external[d]) for a,b,c,d in [(0,1,2,3),(0,2,1,3),(0,3,1,2)]))
    kinetic = {}
    for i,j,k,l in permutations(range(4)):
        pair = tuple(sorted((j,l)))
        kinetic[pair] = kinetic.get(pair,0) - np.trace(_commutator(matrices[i],matrices[j])@_commutator(matrices[k],matrices[l]))/(48*f_pi**2)
    values, rotation = np.linalg.eigh(matrix)
    groups = []
    for index, pole in enumerate(values):
        if groups and pole == values[groups[-1][0]]:
            groups[-1].append(index)
        else:
            groups.append([index])
    grouped_poles = np.array([values[group[0]] for group in groups])
    exchanges = []
    for spectator, left, right in ((1,2,3),(2,1,3),(3,1,2)):
        first = np.array([cubic(external[0],external[spectator],leg) for leg in internal])
        second = np.array([cubic(leg,external[left],external[right]) for leg in internal])
        residues = (first@rotation)*(rotation.T@second)
        exchanges.append((spectator, np.array([sum(residues[group]) for group in groups])))
    return ThreePseudoscalarAmplitude(masses, tuple(mesons), contact, kinetic, grouped_poles, exchanges)


class ThreePseudoscalarAmplitude:
    def __init__(self, masses, mesons, contact, kinetic, poles, exchanges):
        self.masses, self.mesons, self.contact = masses, mesons, contact
        self.kinetic, self.poles, self.exchanges = kinetic, poles, exchanges

    def __call__(self, s12, s23):
        masses2 = [self.masses['a']**2]+[self.masses[x]**2 for x in self.mesons]
        s13 = sum(masses2)-s12-s23
        invariants = {(1,2):s12, (2,3):s23, (1,3):s13}
        products = {pair:(value-masses2[pair[0]]-masses2[pair[1]])/2 for pair,value in invariants.items()}
        for spectator in (1,2,3):
            other = tuple(x for x in (1,2,3) if x != spectator)
            products[(0,spectator)] = -(masses2[0]+masses2[spectator]-invariants[other])/2
        value = self.contact+sum(coef*products[pair] for pair,coef in self.kinetic.items())
        for spectator, residues in self.exchanges:
            other = tuple(x for x in (1,2,3) if x != spectator)
            for pole, residue in zip(self.poles, residues):
                if residue == 0:
                    continue
                den = pole-invariants[other]
                if np.any(den == 0):
                    raise ValueError('internal pole lies at this invariant; stable tree amplitude diverges')
                value = value + residue/den
        return value


def decay_width(f_a, m_a0, hat_S, hat_P, hat_Theta, mesons, *, f_pi, B_0, m_0_squared, epsrel=1e-5):
    """Integrate the squared amplitude over the full Dalitz region, in GeV."""
    masses = calculate_masses(f_a,m_a0,hat_S,hat_P,hat_Theta,f_pi=f_pi,B_0=B_0,m_0_squared=m_0_squared)
    if len(mesons)!=3 or any(name not in masses or name=='a' for name in mesons):
        raise ValueError('provide three valid meson labels')
    charge = {'pi+':1,'pi-':-1,'K+':1,'K-':-1}
    if sum(charge.get(name,0) for name in mesons) != 0:
        return 0.0
    parent = masses['a']; m1,m2,m3 = (masses[name] for name in mesons)
    if parent <= m1+m2+m3:
        return 0.0
    amp = prepare_amplitude(f_a,m_a0,hat_S,hat_P,hat_Theta,mesons,f_pi=f_pi,B_0=B_0,m_0_squared=m_0_squared)
    for spectator, residues in amp.exchanges:
        other = [i for i in (1,2,3) if i != spectator]
        daughter = [m1,m2,m3]
        low = (daughter[other[0]-1]+daughter[other[1]-1])**2
        high = (parent-daughter[spectator-1])**2
        groups = []
        for index, pole in enumerate(amp.poles):
            if groups and pole == amp.poles[groups[-1][0]]:
                groups[-1].append(index)
            else:
                groups.append([index])
        for group in groups:
            pole = amp.poles[group[0]]
            total = sum(residues[group])
            if low <= pole <= high and total != 0:
                raise ValueError('accessible internal pole: stable tree-level width diverges')
    if amp.contact == 0 and all(v == 0 for v in amp.kinetic.values()) and all(np.all(r == 0) for _,r in amp.exchanges):
        return 0.0
    def kallen(x, y, z):
        return (x-y-z)**2-4*y*z
    def limits(s):
        root = np.sqrt(np.maximum(0., kallen(s,m1*m1,m2*m2)*kallen(parent*parent,s,m3*m3)))
        midpoint = m2*m2+m3*m3+(s+m2*m2-m1*m1)*(parent*parent-s-m3*m3)/(2*s)
        half = root/(2*s)
        return midpoint-half, midpoint+half
    def fixed_squared(s):
        return lambda t: np.abs(amp(s,t))**2
    pole_groups = {}
    for spectator, residues in amp.exchanges:
        pole_groups[spectator] = tuple(complex(pole) for pole, residue in zip(amp.poles, residues) if residue != 0)
    total_squared = parent**2+m1**2+m2**2+m3**2
    integral = integrate_three_body(fixed_squared, limits, (m1+m2)**2, (parent-m3)**2,
                                    total_squared, (pole_groups[3], pole_groups[1], pole_groups[2]),
                                    epsrel=epsrel)
    symmetry = np.prod([factorial(count) for count in Counter(mesons).values()])
    return integral/(256*np.pi**3*parent**3*symmetry)
