"""One raw-observable likelihood for the preliminary, GV and BaBar datasets.

Run: .venv/bin/python scripts/fit_babar_kkpi_global_vector_sector.py --output data/vector_sector/babar_kkpi_global_vector_sector_fit.json
Canonical epsilon^0123=+1; GeV. The 73 observed summaries appear once.
Preliminary coefficients initialize the search and are not likelihood priors.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import least_squares
from scipy.stats import chi2
import sympy as sp

sys.path.insert(0,str(Path(__file__).resolve().parent))
import calibrate_odd_low_energy as calibration
import fit_vector_gv as gv_fit
from fit_babar_kkpi_g_coefficients import load_inputs,F,MV,MK,MPI,ALPHA,CHANNELS,triangle,propagator as reference_propagator

ROOT=Path(__file__).resolve().parents[1]
ODD_NAMES=('c1_minus_c2_plus_c5','c1_minus_c2_minus_c5_plus_2c6',
    'c1_plus_c2_minus_c5','c3','c4','d1','d2','d3',
    'g1_plus_2g2_minus_g3','g2','g4','2g4_plus_g5')
EVEN_NAMES=('m_rho_GeV','m_Kstar_GeV','Delta_V_GeV2','kappa_V','F_V_GeV','z_V_GeV_inverse')
FREE_TOTALS=('rho0','rhoplus','omega','phi','etap')
PARAMETER_NAMES=tuple('scaled_'+name for name in EVEN_NAMES)+('G_V_GeV',)+tuple(
    'total_width_pull_'+name for name in FREE_TOTALS)+ODD_NAMES+(
    'gT_times_CT_dimensionless','gT_times_DT_GeV','beta_first_kkpi0','beta_first_kskpi',
    'beta_second_kkpi0','beta_second_kskpi','normalization_pull_kkpi0','normalization_pull_kskpi')
EVEN_SCALES=np.array([.001,.001,.01,.02,.01,.05])
ODD_START=12
TENSOR_START=24
BETA_START=26
NORMALIZATION_START=30
RELATIVE_SYSTEMATICS={'kkpi0':.065,'kskpi':.051}


def running_propagator(name,invariant,beta,mass,shared_width,reference):
    """Rescale the reference running width to the shared on-shell width.

    Masses and widths are in GeV, and the invariant is in GeV squared. The
    zero-excitation propagator is reused only within this call, so changes of
    parameters or quadrature grids cannot reuse a result from an earlier call.
    """
    ground=reference_propagator(name,invariant,reference,0.)
    denominator=1/ground
    on_shell=1/reference_propagator(name,mass**2,reference,0.)
    width=-denominator.imag/mass*shared_width/(-on_shell.imag/mass)
    excited=reference_propagator(name,invariant,reference,beta)-ground
    return 1/(mass**2-invariant-1j*mass*width)+excited


def odd_representative(coordinates):
    """Represent the 12 measured combinations with c5=c7=d4=g3=0."""
    a,b,p,c3,c4,d1,d2,d3,g_a,g2,g4,g_mass=coordinates
    return dict(c1=(a+p)/2,c2=(p-a)/2,c3=c3,c4=c4,c5=0.,
        c6=(b-a)/2,c7=0.,d1=d1,d2=d2,d3=d3,d4=0.,
        g1=g_a-2*g2,g2=g2,g3=0.,g4=g4,g5=g_mass-2*g4)


def original_odd_matrix():
    """Linear representative map to the eleven original VJP/VVP coefficients."""
    matrix=np.zeros((11,12))
    for index in range(12):
        representative=odd_representative(np.eye(12)[index])
        matrix[:,index]=[representative[key] for key in [*(f'c{i}' for i in range(1,8)),*(f'd{i}' for i in range(1,5))]]
    return matrix


def exact_odd_map():
    """Compute the original sixteen-coefficient map and structural nulls exactly."""
    c=sp.symbols('c1:8');d=sp.symbols('d1:5');g=sp.symbols('g1:6')
    values=sp.Matrix([c[0]-c[1]+c[4],c[0]-c[1]-c[4]+2*c[5],
        c[0]+c[1]-c[4],c[2],c[3],d[0],d[1],d[2],
        g[0]+2*g[1]-g[2],g[1],g[3],2*g[3]+g[4]])
    matrix=values.jacobian(c+d+g)
    return dict(original_order=list(map(str,c+d+g)),coordinate_order=list(ODD_NAMES),
        exact_rank=int(matrix.rank()),exact_null_basis=[[str(v) for v in vector] for vector in matrix.nullspace()],
        coordinate_matrix=[[str(v) for v in row] for row in matrix.tolist()])


class EnergyBinBlocks:
    """Quadrature orders and fixed energy-bin slices; grids are built on use."""
    def __init__(self,data,orders,block_size):
        self.data=data
        self.orders=tuple(orders)
        self.block_size=block_size

    def blocks(self):
        for start in range(0,max(len(self.data[channel]['edges']) for channel in CHANNELS),self.block_size):
            yield {channel:{name:(value[start:start+self.block_size]
                if isinstance(value,np.ndarray) and value.ndim and len(value)==len(self.data[channel]['edges'])
                else value) for name,value in self.data[channel].items()} for channel in CHANNELS}


def phase_points(data,energy_order=4,pair_order=32,angle_order=32):
    """Bin-averaged Dalitz weights in nb, with exact phi/Kstar pole mappings.

    Three KK-invariant segments resolve the narrow phi. Both angular halves
    use arctangent changes of variable around the two Kstar poles. Fixed
    reference pole positions affect numerical convergence only.
    """
    energy_nodes,energy_weights=leggauss(energy_order)
    pair_nodes,pair_weights=leggauss(pair_order)
    angles,angle_weights=leggauss(angle_order)
    result={}
    for channel in CHANNELS:
        masses=(.493677,.493677,.1349768) if channel=='kkpi0' else (.497611,.493677,.13957039)
        first,second,third=masses
        lower,upper=data[channel]['edges'].T
        energies=(lower[:,None]+upper[:,None])/2+(upper-lower)[:,None]*energy_nodes/2
        q=energies**2
        s_lower=(first+second)**2
        s_upper=(energies-third)**2
        phi_pole=1.019461**2;phi_scale=1.019461*.004249
        boundaries=np.stack((np.full_like(s_upper,s_lower),np.clip(np.full_like(s_upper,phi_pole-10*phi_scale),s_lower,s_upper),np.clip(np.full_like(s_upper,phi_pole+10*phi_scale),s_lower,s_upper),s_upper),axis=-1)
        segment_lower=boundaries[...,:-1];segment_upper=boundaries[...,1:]
        pair_values=(segment_lower+segment_upper)[...,None]/2+(segment_upper-segment_lower)[...,None]*pair_nodes/2
        pair_measure=(segment_upper-segment_lower)[...,None]/2*np.broadcast_to(pair_weights,pair_values.shape)
        low_angle=np.arctan((segment_lower[...,1]-phi_pole)/phi_scale)
        high_angle=np.arctan((segment_upper[...,1]-phi_pole)/phi_scale)
        phi_angle=(low_angle+high_angle)[...,None]/2+(high_angle-low_angle)[...,None]*pair_nodes/2
        pair_values[...,1,:]=phi_pole+phi_scale*np.tan(phi_angle)
        pair_measure[...,1,:]=phi_scale/np.cos(phi_angle)**2*(high_angle-low_angle)[...,None]/2*pair_weights
        pair_s=pair_values.reshape(*s_upper.shape,3*pair_order)
        pair_measure=pair_measure.reshape(*s_upper.shape,3*pair_order)
        outer=np.maximum(triangle(q[...,None],pair_s,third**2),0)
        inner=np.maximum(triangle(pair_s,first**2,second**2),0)
        span=np.sqrt(outer*inner)/pair_s
        pair_s13_center=(first**2+third**2+(pair_s+first**2-second**2)
            *(q[...,None]-pair_s-third**2)/(2*pair_s))
        t_center=q[...,None]+sum(mass**2 for mass in masses)-pair_s-pair_s13_center
        t_lower=t_center-span/2;t_upper=t_center+span/2
        invariant_sum=q[...,None]+sum(mass**2 for mass in masses)-pair_s
        # An exact change of integration variable resolves both Kstar poles.
        # Fixed reference values determine the numerical mapping only.
        pole=.8917**2;pole_scale=.8917*.05
        centers=np.stack((np.minimum(pole,invariant_sum-pole),np.maximum(pole,invariant_sum-pole)),axis=-1)
        limits_lower=np.stack((t_lower,t_center),axis=-1)
        limits_upper=np.stack((t_center,t_upper),axis=-1)
        angle_lower=np.arctan((limits_lower-centers)/pole_scale)
        angle_upper=np.arctan((limits_upper-centers)/pole_scale)
        mapped_angle=(angle_lower+angle_upper)[...,None]/2+(angle_upper-angle_lower)[...,None]*angles/2
        mapped_t=centers[...,None]+pole_scale*np.tan(mapped_angle)
        t_weights=pole_scale/np.cos(mapped_angle)**2*(angle_upper-angle_lower)[...,None]/2*angle_weights
        pair_t=mapped_t.reshape(*pair_s.shape,2*angle_order)
        t_weights=t_weights.reshape(*pair_s.shape,2*angle_order)
        cosine=2*(pair_t-t_center[...,None])/np.maximum(span[...,None],1e-30)
        current_squared=(outer*inner/(16*pair_s))[...,None]*np.maximum(1-cosine**2,0)
        measure=pair_measure[...,None]*t_weights
        weight=(389379.3656*ALPHA**2/(48*np.pi*q[...,None,None]**3)*current_squared*measure
            *energy_weights[None,:,None,None]/2)
        if channel=='kskpi':weight*=2
        result[channel]=dict(q=q[...,None,None],s=pair_s[...,None],t=pair_t,
            weight=weight*np.exp(-2*q[...,None,None]**3/2.5**6),shape=pair_t.shape)
    return result


class GlobalLikelihood:
    """All experimental blocks and every shared nuisance in one residual vector."""
    def __init__(self,energy_order=2,pair_order=16,angle_order=16,bin_block_size=None):
        self.preliminary,self.babar=load_inputs()
        self.odd_data=json.loads((ROOT/'data/odd_parity/odd_low_energy_measurements.json').read_text())
        self.observations=[x for x in self.odd_data['observations'] if x.get('status')!='deferred']
        self.stage3=json.loads((ROOT/'data/vector_sector/stage3_results.json').read_text())
        self.assumptions=self.stage3['assumptions']
        self.pdg=json.loads((ROOT/'data/vector_sector/pdg2024.json').read_text())
        self.alda=json.loads((ROOT/'data/vector_sector/alda_2026nrd_arxiv_v1.json').read_text())
        self.gv_data=json.loads((ROOT/'data/vector_sector/gv_measurements.json').read_text())
        self.point=np.array([self.preliminary[key] for key in EVEN_NAMES])
        self.tree=calibration.J.module('alda_neutral_inputs').alda_2026nrd_neutral_inputs(mass_scheme='tree',F=F,fpi=F)
        self.fast=calibration.FastChannels(self.tree,self.alda['physical_masses_gev'],
            np.sqrt(4*np.pi*self.assumptions['alpha_em']))
        self.embedding=calibration.model_maps()['original_11']
        self.odd_original_map=original_odd_matrix()
        self.odd_target,self.odd_sigma=np.array([calibration.J.data_units(x) for x in self.observations]).T
        self.even_target=np.array([self.pdg['masses'][key]['value_mev']/1000 for key in
            ('rho_charged','kstar_charged','omega','phi')]+[
            self.pdg['leptonic_widths'][key]['value_kev']*1e-6 for key in ('rho','omega','phi')])
        self.even_sigma=np.array([self.pdg['masses'][key]['sigma_mev']/1000 for key in
            ('rho_charged','kstar_charged','omega','phi')]+[
            self.pdg['leptonic_widths'][key]['sigma_kev']*1e-6 for key in ('rho','omega','phi')])
        self.total_order=tuple(sorted(self.odd_data['total_widths']))
        self.total_target=np.array([self.odd_data['total_widths'][name]['value_gev'] for name in self.total_order])
        self.total_sigma=np.array([self.odd_data['total_widths'][name]['sigma_gev'] for name in self.total_order])
        self.phi_target=np.array([row['value'] for row in self.gv_data['rows'][2:]])
        self.phi_cholesky=np.linalg.cholesky(np.array(self.gv_data['covariance'])[2:,2:])
        self.kaon_flavors=gv_fit.flavor_factors(self.preliminary['omega_phi_mixing_angle_rad'])
        for row,total in zip(self.gv_data['rows'][:2],('Kstarplus','Kstar0')):
            width=self.odd_data['total_widths'][total]
            assert row['value']==width['value_gev'] and row['sigma']==width['sigma_gev']
        self.bin_block_size=bin_block_size
        self.phase=self.quadrature_phase((energy_order,pair_order,angle_order))
        self.cache=None
        assert len(self.observations)==12 and sum(len(self.babar[c]['observed']) for c in CHANNELS)==45
        self.measurement_count=12+7+7+2+45
        assert self.measurement_count==73
        midpoint=4*(self.point[1]**2-self.point[0]**2)/(3*(self.assumptions['xl_gev2']-self.assumptions['xs_gev2']))
        self.branch_sign=np.sign(self.point[3]-midpoint)

    def unpack(self,values):
        """Recompute the shared mass/source sector and all seven total widths."""
        point=self.point+EVEN_SCALES*values[:6]
        sector,mk=calibration.J.action_sector(point,'single_trace_negative',self.assumptions)
        midpoint=4*(point[1]**2-point[0]**2)/(3*(self.assumptions['xl_gev2']-self.assumptions['xs_gev2']))
        if np.sign(point[3]-midpoint)!=self.branch_sign:raise ValueError('inverse-mass branch crossed')
        light=point[4]+2*point[5]*self.assumptions['xl_gev2']
        strange=point[4]+2*point[5]*self.assumptions['xs_gev2']
        if light<=0 or strange>=0:raise ValueError('selected source branch crossed')
        gv=values[6]
        if gv<=0 or np.any(abs(values[26:30])>5):raise ValueError('selected GV/pole-residue domain crossed')
        totals={name:self.odd_data['total_widths'][name]['value_gev']
            +self.odd_data['total_widths'][name]['sigma_gev']*values[7+index]
            for index,name in enumerate(FREE_TOTALS)}
        physical=self.alda['physical_masses_gev']
        for row,name in zip(self.gv_data['rows'][:2],('Kstarplus','Kstar0')):
            modes=([('K_charged','pi0','kstar_charged_K+_pi0'),('K_neutral','pi_charged','kstar_charged_K0_pi+')]
                if name=='Kstarplus' else [('K_neutral','pi0','kstar_neutral_K0_pi0'),('K_charged','pi_charged','kstar_neutral_K+_pi-')])
            totals[name]=gv**2*sum(gv_fit.width_coefficient(row['physical_vector_mass_GeV'],mk,
                physical[first],physical[second],self.kaon_flavors[flavor]) for first,second,flavor in modes)
        if min(totals.values())<=0:raise ValueError('nonpositive total width')
        return point,sector,mk,gv,totals

    def old_blocks(self,values,state):
        """Return the 28 unique preliminary/GV summaries and odd derivatives."""
        point,sector,mk,gv,totals=state
        channels=self.fast(sector,mk)
        original=odd_representative(values[ODD_START:ODD_START+12])
        coefficients=np.array([original[key] for key in [*(f'c{i}' for i in range(1,8)),*(f'd{i}' for i in range(1,5))]])
        predictions=[];derivatives=[];adapter_defect=0.
        for row in self.observations:
            channel=channels[row['channel']]
            gradient=channel['row']@self.embedding
            amplitude=-channel['offset']+gradient@coefficients
            reference=channel['offset']+gradient@(-coefficients)
            adapter_defect=max(adapter_defect,abs(amplitude**2-reference**2))
            factor=channel['phase_factor']/(totals[row['total_width_key']] if row['observable']=='branching_fraction' else 1.)
            predictions.append(factor*amplitude**2)
            derivatives.append(2*factor*amplitude*(gradient@self.odd_original_map))
        odd=(np.array(predictions)-self.odd_target)/self.odd_sigma
        odd_jac=np.array(derivatives)/self.odd_sigma[:,None]
        even=(calibration.J.even_observables(sector,mk,self.assumptions)-self.even_target)/self.even_sigma
        total=(np.array([totals[name] for name in self.total_order])-self.total_target)/self.total_sigma
        phi=sector.generators[2]
        phi_traces=[float((phi[index,index]-phi[2,2]).real)**2 for index in (0,1)]
        phi_mass=sector.masses[2]
        phi_predictions=[]
        for row,daughter,phi_trace in zip(self.gv_data['rows'][2:],('K_charged','K_neutral'),phi_traces):
            phi_predictions.append(gv**2*gv_fit.width_coefficient(row['physical_vector_mass_GeV'],phi_mass,
                self.alda['physical_masses_gev'][daughter],self.alda['physical_masses_gev'][daughter],phi_trace)/totals['phi'])
        phi_residual=np.linalg.solve(self.phi_cholesky,np.array(phi_predictions)-self.phi_target)
        return np.r_[odd,even,total,phi_residual],odd_jac,dict(
            odd_predictions=list(map(float,predictions)),phi_branching_predictions=phi_predictions,
            canonical_reference_squared_amplitude_defect=adapter_defect,shared_totals=totals)

    def quadrature_phase(self,orders):
        """Choose full grids or bounded-memory energy-bin evaluation."""
        if self.bin_block_size is None:return phase_points(self.babar,*orders)
        return EnergyBinBlocks(self.babar,orders,self.bin_block_size)

    def babar_blocks(self,values,state,derivatives=False):
        """Integrate each bin in fixed order, optionally building one block at a time."""
        if not isinstance(self.phase,EnergyBinBlocks):
            return self._babar_blocks_full(values,state,derivatives)
        original_phase=self.phase;original_data=self.babar
        channel_predictions={channel:[] for channel in CHANNELS}
        try:
            for data in original_phase.blocks():
                self.babar=data
                self.phase=phase_points(data,*original_phase.orders)
                self.babar_cache_key=None
                _,predictions=self._babar_blocks_full(values,state)
                if derivatives:
                    rows=self._babar_blocks_full(values,state,True)
                    offset=0
                    for channel in CHANNELS:
                        count=len(data[channel]['edges'])
                        channel_predictions[channel].append(rows[offset:offset+count])
                        offset+=count
                else:
                    for channel in CHANNELS:channel_predictions[channel].append(predictions[channel])
            predictions={channel:np.concatenate(parts,axis=0)
                for channel,parts in channel_predictions.items()}
            if derivatives:return np.vstack([predictions[channel] for channel in CHANNELS])
            residual=np.concatenate([(predictions[channel]*(1+RELATIVE_SYSTEMATICS[channel]*values[30+index])
                -original_data[channel]['observed'])/original_data[channel]['errors']
                for index,channel in enumerate(CHANNELS)])
            return residual,predictions
        finally:
            self.phase=original_phase;self.babar=original_data
            self.babar_cache_key=None;self.babar_cached={};self.current_amplitudes={}

    def _babar_blocks_full(self,values,state,derivatives=False):
        """Integrate the coherent current with the fitted neutral response.

        The ground-state light/strange propagator is obtained by rotating both
        ends of the physical-pole propagator. Excited multiplets retain the
        ideal flavor assignment of Qin:2024ulb. Running rho/Kstar widths have
        the reference threshold shapes, normalized to the shared on-shell
        widths; the Kstar width is the charged/neutral mean (isospin limit).
        """
        point,sector,mk,gv,totals=state
        a,b,p,c3,c4,d1,d2,d3,g,g2,g4,bg=(np.eye(12).reshape(12,12,1,1,1,1) if derivatives else values[12:24])
        cm=p+8*c3;dm=d1+8*d2-d3
        light=point[4]+2*point[5]*self.assumptions['xl_gev2']
        strange=point[4]+2*point[5]*self.assumptions['xs_gev2']
        reference=dict(self.preliminary,m_rho_GeV=point[0],m_Kstar_GeV=mk,
            Gamma_total_omega_GeV=totals['omega'],Gamma_total_phi_GeV=totals['phi'])
        light_generator=np.diag([1,1,0])/np.sqrt(2)
        strange_generator=np.diag([0,0,1])
        rotation=np.array([[np.trace(generator@basis).real for generator in sector.generators[1:]]
            for basis in (light_generator,strange_generator)])
        def excited(name,invariant,beta):
            return reference_propagator(name,invariant,reference,beta)-reference_propagator(name,invariant,reference,0.)
        def neutral(invariant,beta):
            omega=1/(sector.masses[1]**2-invariant-1j*sector.masses[1]*totals['omega'])
            phi=1/(sector.masses[2]**2-invariant-1j*sector.masses[2]*totals['phi'])
            ll=rotation[0,0]**2*omega+rotation[0,1]**2*phi+excited('omega',invariant,beta)
            ss=rotation[1,0]**2*omega+rotation[1,1]**2*phi+excited('phi',invariant,beta)
            ls=rotation[0,0]*rotation[1,0]*omega+rotation[0,1]*rotation[1,1]*phi
            return ll,ss,ls
        def running(name,invariant,beta):
            mass=point[0] if name=='rho' else mk
            shared=totals['rho0'] if name=='rho' else (totals['Kstar0']+totals['Kstarplus'])/2
            return running_propagator(name,invariant,beta,mass,shared,reference)
        cache_key=tuple(values[:12])+tuple(values[26:30])
        if getattr(self,'babar_cache_key',None)!=cache_key:
            self.babar_cache_key=cache_key;self.babar_cached={}
        if not derivatives:self.current_amplitudes={}
        predictions={}
        for index,channel in enumerate(CHANNELS):
            phase=self.phase[channel];q,s,t=phase['q'],phase['s'],phase['t']
            u=q+2*MK**2+MPI**2-s-t
            beta=(values[26+index],values[28+index])
            if channel not in self.babar_cached:
                rq=running('rho',q,beta);rs=running('rho',s,beta)
                llq,ssq,lsq=neutral(q,beta);lls,_,lss=neutral(s,beta)
                lls=lls-np.sqrt(2)*lss
                kt=running('kstar',t,beta);ku=running('kstar',u,beta)
                source_rho=light*rq
                source_light=light*llq-np.sqrt(2)*strange*lsq
                source_strange=strange*ssq-light*lsq/np.sqrt(2)
                self.babar_cached[channel]=(rs,lls,kt,ku,source_rho,source_light,source_strange)
            rs,lls,kt,ku,source_rho,source_light,source_strange=self.babar_cached[channel]
            def cr(pair,mass2):return a*q-b*pair+cm*mass2
            def crk(pair):return cr(pair,MK**2)+24*c4*(MK**2-MPI**2)
            def dr(pair,mass2):return dm*mass2+d3*(q+pair)
            gr1=g*(q+s-4*MK**2-MPI**2)+4*bg*MK**2-4*g2*(q-2*MK**2-MPI**2)-4*g4*(MK**2-MPI**2)
            gr2=g*(q-s-MPI**2)+2*bg*MPI**2-2*g2*(q-2*MK**2-MPI**2)+4*g4*(MK**2-MPI**2)
            contact=2*np.sqrt(2)/(F**3*MV);direct=2*np.sqrt(2)*gv/(3*F**3*MV);exchange=2*gv/F**3
            if channel=='kkpi0':
                amplitude=(0. if derivatives else 3/(12*np.pi**2*F**3))+contact*((source_rho+source_light/3)*gr1-2*source_strange*gr2/3)
                amplitude=amplitude-direct*((rs+3*lls)*cr(s,MPI**2)+crk(t)*kt+crk(u)*ku)
                amplitude=amplitude+exchange*(2*dr(s,MPI**2)*(source_rho*lls+source_light*rs/3)+(source_rho+source_light/3-2*source_strange/3)*(dr(t,MK**2)*kt+dr(u,MK**2)*ku))
            else:
                amplitude=contact*(g*(t-u)*source_rho-source_light*gr1/3+2*source_strange*gr2/3)
                amplitude=amplitude+direct*(rs*cr(s,MPI**2)-2*cr(t,MK**2)*kt+crk(u)*ku)
                amplitude=amplitude-exchange*(2*source_light*rs*dr(s,MPI**2)/3-(source_rho-source_light/3+2*source_strange/3)*dr(t,MK**2)*kt+(source_rho+source_light/3-2*source_strange/3)*dr(u,MK**2)*ku)
                if not derivatives:
                    amplitude=amplitude+(t-u)/(1.3182**2-s-1j*1.3182*.107)*(-2*values[24]/F**3+np.sqrt(2)*values[25]*source_rho/F**3)
            if derivatives:
                predictions[channel]=np.sum(phase['weight'][None,...]*2*np.real(self.current_amplitudes[channel][None,...].conj()*amplitude),axis=(2,3,4)).T
            else:
                self.current_amplitudes[channel]=amplitude
                predictions[channel]=np.sum(phase['weight']*abs(amplitude)**2,axis=(1,2,3))
        if derivatives:
            return np.vstack([predictions[channel]*(1+RELATIVE_SYSTEMATICS[channel]*values[30+index])/self.babar[channel]['errors'][:,None] for index,channel in enumerate(CHANNELS)])
        residual=np.concatenate([(predictions[channel]*(1+RELATIVE_SYSTEMATICS[channel]*values[30+index])-self.babar[channel]['observed'])/self.babar[channel]['errors'] for index,channel in enumerate(CHANNELS)])
        return residual,predictions

    def residual(self,values):
        """Concatenate 73 observed summaries and two normalization constraints."""
        try:
            state=self.unpack(values)
            old,_,_=self.old_blocks(values,state)
            babar,_=self.babar_blocks(values,state)
            return np.r_[old,babar,values[30:32]]
        except (ValueError,FloatingPointError,np.linalg.LinAlgError):
            return np.full(75,1e8+np.linalg.norm(values))

    def jacobian(self,values,finite_difference_steps=None):
        """Exact odd derivatives; forward differences for remaining parameters."""
        baseline=self.residual(values)
        jacobian=np.zeros((75,32))
        state=self.unpack(values)
        _,old_odd,_=self.old_blocks(values,state)
        jacobian[:12,12:24]=old_odd
        jacobian[28:73,12:24]=self.babar_blocks(values,state,derivatives=True)
        for index in list(range(12))+list(range(24,32)):
            step=(np.sqrt(np.finfo(float).eps)*max(1.,abs(values[index]))
                if finite_difference_steps is None else finite_difference_steps[index])
            shifted=values.copy();shifted[index]+=step
            jacobian[:,index]=(self.residual(shifted)-baseline)/step
        return jacobian

    def initial(self):
        """Use previous coefficients only to initialize the nonlinear search."""
        values=np.zeros(32);values[6]=self.preliminary['G_V_GeV']
        for index,name in enumerate(FREE_TOTALS):
            values[7+index]=(self.preliminary['Gamma_total_'+name+'_GeV']-self.odd_data['total_widths'][name]['value_gev'])/self.odd_data['total_widths'][name]['sigma_gev']
        c=self.preliminary
        values[12:24]=[c['c1']-c['c2']+c['c5'],c['c1']-c['c2']-c['c5']+2*c['c6'],c['c1']+c['c2']-c['c5'],c['c3'],c['c4'],c['d1'],c['d2'],c['d3'],0,0,0,0]
        return values


QUADRATURE_TOLERANCES=dict(prediction_sigma=.01,residual=.01,
    jacobian_relative=1e-3,jacobian_absolute=1e-3,jacobian_negligible=1e-12,chi2=.01)


def compare_quadrature_outputs(previous,current,errors,parameter_scales,tolerances):
    """Compare nb predictions and dimensionless residuals/scaled derivatives."""
    prediction_changes={channel:float(np.max(np.abs(current['predictions'][channel]
        -previous['predictions'][channel])/errors[channel])) for channel in CHANNELS}
    residual_change=float(np.max(np.abs(current['residual']-previous['residual'])))
    previous_jacobian=previous['jacobian']*parameter_scales[None,:]
    current_jacobian=current['jacobian']*parameter_scales[None,:]
    column_changes=np.linalg.norm(current_jacobian-previous_jacobian,axis=0)
    column_norms=np.linalg.norm(current_jacobian,axis=0)
    previous_norms=np.linalg.norm(previous_jacobian,axis=0)
    columns=[]
    for index,(change,norm) in enumerate(zip(column_changes,column_norms)):
        comparison_norm=max(norm,previous_norms[index])
        negligible=bool(comparison_norm<=tolerances['jacobian_negligible'])
        relative=None if negligible else float(change/comparison_norm)
        passed=bool(change<=tolerances['jacobian_absolute'] if negligible
            else relative<=tolerances['jacobian_relative'])
        columns.append(dict(parameter=PARAMETER_NAMES[index],absolute_change=float(change),
            previous_norm=float(previous_norms[index]),refined_norm=float(norm),relative_change=relative,
            comparison='absolute' if negligible else 'relative',passed=passed))
    chi2_change=float(abs(current['chi2']-previous['chi2']))
    comparisons=dict(predictions=dict(maximum_change_in_sigma=max(prediction_changes.values()),
            by_channel=prediction_changes,passed=bool(max(prediction_changes.values())<=tolerances['prediction_sigma'])),
        residuals=dict(maximum_absolute_change=residual_change,passed=bool(residual_change<=tolerances['residual'])),
        jacobian=dict(columns=columns,passed=all(column['passed'] for column in columns)),
        chi2=dict(absolute_change=chi2_change,passed=bool(chi2_change<=tolerances['chi2'])))
    return dict(comparisons,passed=all(comparison['passed'] for comparison in comparisons.values()))


def audit_quadrature(model,values,orders,parameter_scales=None,tolerances=None):
    """Refine each dimension separately and then all dimensions together.

    Each order sequence includes its baseline and at least one larger order.
    Every reported successive comparison must pass. All derivatives use
    the same steps and scales.
    """
    dimensions=('energy','pair','angle')
    orders={name:tuple(orders[name]) for name in dimensions}
    for name,sequence in orders.items():
        if len(sequence)<2 or any(int(order)!=order or order<1 for order in sequence) or any(
                higher<=lower for lower,higher in zip(sequence,sequence[1:])):
            raise ValueError(f'{name} requires at least two strictly increasing positive integer orders')
    tolerances=dict(QUADRATURE_TOLERANCES,**(tolerances or {}))
    parameter_scales=(np.ones(len(values)) if parameter_scales is None else np.asarray(parameter_scales))
    steps=np.sqrt(np.finfo(float).eps)*np.maximum(1.,np.abs(values))
    errors={channel:model.babar[channel]['errors'] for channel in CHANNELS}
    baseline=tuple(orders[name][0] for name in dimensions)
    original_phase=model.phase
    evaluations={}
    refinements={}
    attempted_orders=None
    def evaluate(order_tuple):
        nonlocal attempted_orders
        if order_tuple not in evaluations:
            attempted_orders=dict(zip(dimensions,order_tuple))
            model.phase=(model.quadrature_phase(order_tuple) if hasattr(model,'quadrature_phase')
                else phase_points(model.babar,*order_tuple))
            model.babar_cache_key=None
            residual=model.residual(values)
            jacobian=model.jacobian(values,finite_difference_steps=steps)
            _,predictions=model.babar_blocks(values,model.unpack(values))
            evaluations[order_tuple]=dict(residual=residual.copy(),jacobian=jacobian.copy(),
                predictions={channel:prediction.copy() for channel,prediction in predictions.items()},
                chi2=float(residual@residual))
        return evaluations[order_tuple]
    def compare(lower,higher):
        result=compare_quadrature_outputs(evaluate(lower),evaluate(higher),errors,parameter_scales,tolerances)
        return dict(from_orders=dict(zip(dimensions,lower)),to_orders=dict(zip(dimensions,higher)),**result)
    def compare_sequence(name,sequence):
        refinements[name]=dict(comparisons=[],passed=False)
        for lower,higher in zip(sequence,sequence[1:]):
            refinements[name]['comparisons'].append(compare(lower,higher))
        refinements[name]['passed']=all(comparison['passed']
            for comparison in refinements[name]['comparisons'])
    try:
        for index,name in enumerate(dimensions):
            sequence=[]
            for order in orders[name]:
                current=list(baseline);current[index]=order;sequence.append(tuple(current))
            compare_sequence(name,sequence)
        combined=[tuple(orders[name][min(level,len(orders[name])-1)] for name in dimensions)
            for level in range(max(map(len,orders.values()))) ]
        compare_sequence('combined',combined)
        passed=all(result['passed'] for result in refinements.values())
        return dict(status='certified' if passed else 'uncertified',passed=passed,
            tolerances=tolerances,orders={name:list(sequence) for name,sequence in orders.items()},
            derivative_settings=dict(method='analytic odd columns; forward differences for remaining columns',
                finite_difference_steps=steps.tolist(),parameter_scales=parameter_scales.tolist()),
            pole_mapping='phi s-pole mapping and both Kstar t/u-pole mappings',refinements=refinements)
    except (ValueError,FloatingPointError,np.linalg.LinAlgError) as error:
        return dict(status='unevaluated',passed=False,reason=str(error),
            error_type=type(error).__name__,attempted_orders=attempted_orders,
            orders={name:list(sequence) for name,sequence in orders.items()},
            refinements=refinements)
    finally:
        model.phase=original_phase
        model.babar_cache_key=None


def projected_rank(jacobian):
    """Project all non-odd derivatives before measuring odd information."""
    odd=jacobian[:,12:24]
    nuisance=np.column_stack((jacobian[:,:12],jacobian[:,24:]))
    projected=odd-nuisance@np.linalg.lstsq(nuisance,odd,rcond=1e-12)[0]
    singular=np.linalg.svd(projected,compute_uv=False)
    return dict(singular_values=singular.tolist(),rank_by_relative_threshold={str(t):int(np.sum(singular>singular[0]*t)) for t in (1e-6,1e-8,1e-10)})


def even_gv_residual_diagnostic(model,values):
    """Diagnose shared even/GV tension without changing the global objective.

    The twelve selected raw rows are refitted on their own solely to locate
    the incompatible observations. Their minimum is never used as data or
    as a prior in the joint fit. Two independent one-parameter ratio fits
    also give exact minimum chi-squares for the charge pairs.
    """
    point,sector,mk,gv,totals=model.unpack(values)
    kaon_coefficients=np.array([totals[name]*mk**2/gv**2 for name in ('Kstarplus','Kstar0')])
    kaon_target=np.array([row['value'] for row in model.gv_data['rows'][:2]])
    kaon_error=np.array([row['sigma'] for row in model.gv_data['rows'][:2]])
    kaon_scale=np.sum(kaon_coefficients*kaon_target/kaon_error**2)/np.sum((kaon_coefficients/kaon_error)**2)
    kaon_residual=(kaon_scale*kaon_coefficients-kaon_target)/kaon_error
    physical=model.alda['physical_masses_gev']
    phi_coefficients=np.array([gv_fit.width_coefficient(row['physical_vector_mass_GeV'],1.,physical[name],physical[name],1.) for row,name in zip(model.gv_data['rows'][2:],('K_charged','K_neutral'))])
    phi_inverse=np.linalg.inv(model.phi_cholesky@model.phi_cholesky.T)
    phi_scale=phi_coefficients@phi_inverse@model.phi_target/(phi_coefficients@phi_inverse@phi_coefficients)
    phi_residual=np.linalg.solve(model.phi_cholesky,phi_scale*phi_coefficients-model.phi_target)
    initial=values.copy();initial[7:12]=0.
    free=np.r_[np.arange(7),10]
    selected=np.r_[np.arange(12,19),[19+model.total_order.index(name) for name in ('Kstar0','Kstarplus','phi')],26,27]
    def residual(free_values):
        current=initial.copy();current[free]=free_values
        try:return model.old_blocks(current,model.unpack(current))[0][selected]
        except ValueError:return np.full(12,1e8+np.linalg.norm(free_values))
    diagnostic=least_squares(residual,initial[free],x_scale='jac',max_nfev=1000,ftol=1e-10,xtol=1e-10,gtol=1e-8)
    return dict(scope='Diagnostic only; never added to the joint likelihood.',
        selected_residual_indices=selected.tolist(),chi2=float(diagnostic.fun@diagnostic.fun),residuals=diagnostic.fun.tolist(),
        success=bool(diagnostic.success),nfev=int(diagnostic.nfev),message=diagnostic.message,
        exact_Kstar_pair_minimum_chi2=float(kaon_residual@kaon_residual),
        predicted_Kstar_neutral_over_charged_ratio=float(kaon_coefficients[1]/kaon_coefficients[0]),
        observed_Kstar_neutral_over_charged_ratio=float(kaon_target[1]/kaon_target[0]),
        exact_phi_pair_minimum_chi2=float(phi_residual@phi_residual),
        predicted_phi_charged_over_neutral_ratio=float(phi_coefficients[0]/phi_coefficients[1]),
        observed_phi_charged_over_neutral_ratio=float(model.phi_target[0]/model.phi_target[1]))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',default='data/vector_sector/babar_kkpi_global_vector_sector_fit.json')
    parser.add_argument('--max-nfev',type=int,default=1500)
    parser.add_argument('--method',choices=('trf','lm'),default='trf')
    parser.add_argument('--svd-scale',action='store_true')
    parser.add_argument('--ftol',type=float,default=1e-8)
    parser.add_argument('--energy-order',type=int,default=2)
    parser.add_argument('--pair-order',type=int,default=16)
    parser.add_argument('--angle-order',type=int,default=16)
    parser.add_argument('--start')
    parser.add_argument('--initialization',choices=('preliminary','babar','negative_c4','negative_odd'),default='preliminary')
    parser.add_argument('--profiles',action='store_true')
    parser.add_argument('--g-profiles',action='store_true')
    parser.add_argument('--audit-quadrature',action='store_true')
    for dimension in ('energy','pair','angle'):
        parser.add_argument(f'--audit-{dimension}-orders',type=int,nargs='+',
            help='Increasing orders including the fit order; default: fit order, twice, four times')
    parser.add_argument('--profile-max-nfev',type=int,default=200)
    args=parser.parse_args()
    model=GlobalLikelihood(args.energy_order,args.pair_order,args.angle_order)
    initial=model.initial()
    if args.start:initial=np.array(json.loads(Path(args.start).read_text())['parameter_values'])
    if args.initialization=='babar':
        previous=json.loads((ROOT/'data/vector_sector/babar_kkpi_unified_vector_coefficients_fit.json').read_text())['fit']['parameters']
        c3=model.preliminary['c3'];d2=model.preliminary['d2']
        initial[12:24]=[previous[ODD_NAMES[0]],previous[ODD_NAMES[1]],previous['c1_plus_c2_plus_8c3_minus_c5']-8*c3,c3,previous['c4'],previous['d1_plus_8d2_minus_d3']-8*d2+previous['d3'],d2,previous['d3'],*(previous[name] for name in ODD_NAMES[8:])]
        initial[24:32]=[previous[name] for name in PARAMETER_NAMES[24:32]]
    elif args.initialization=='negative_c4':initial[16]*=-1
    elif args.initialization=='negative_odd':initial[12:24]*=-1
    lower=np.full(32,-np.inf);upper=np.full(32,np.inf)
    lower[:6]=[-100,-100,-20,-10,-15,-10];upper[:6]=[100,100,20,10,15,10]
    lower[6]=.005;upper[6]=.2
    for index,name in enumerate(FREE_TOTALS):lower[7+index]=-.99*model.odd_data['total_widths'][name]['value_gev']/model.odd_data['total_widths'][name]['sigma_gev']
    lower[26:30]=-5;upper[26:30]=5
    destination=ROOT/args.output
    destination.parent.mkdir(parents=True,exist_ok=True)
    def checkpoint(intermediate_result):
        destination.with_suffix('.progress.json').write_text(json.dumps(dict(parameter_values=intermediate_result.x.tolist(),chi2=float(intermediate_result.fun@intermediate_result.fun),nfev=int(intermediate_result.nfev))))
    best_cost=[np.inf]
    def fit_residual(values):
        if np.any(values<lower) or np.any(values>upper):return np.full(75,1e8+np.linalg.norm(values))
        residual=model.residual(values)
        cost=float(residual@residual)
        if cost<best_cost[0]:
            best_cost[0]=cost
            destination.with_suffix('.progress.json').write_text(json.dumps(dict(parameter_values=values.tolist(),chi2=cost)))
        return residual
    options=dict(jac=model.jacobian,method=args.method,x_scale='jac',max_nfev=args.max_nfev,ftol=args.ftol,xtol=1e-9,gtol=1e-7)
    if args.method=='trf':options.update(callback=checkpoint,bounds=(lower,upper))
    if args.svd_scale:
        initial_jacobian=model.jacobian(initial)
        initial_norms=np.linalg.norm(initial_jacobian,axis=0)
        _,initial_singular,initial_right=np.linalg.svd(initial_jacobian/initial_norms,full_matrices=False)
        transformation=(initial_right.T/np.maximum(initial_singular,1e-3))/initial_norms[:,None]
        def transformed_residual(coordinates):return fit_residual(initial+transformation@coordinates)
        def transformed_jacobian(coordinates):return model.jacobian(initial+transformation@coordinates)@transformation
        options.pop('bounds',None);options.pop('callback',None)
        options.update(jac=transformed_jacobian,x_scale=1.)
        fit=least_squares(transformed_residual,np.zeros(32),**options)
        scaled_optimality=float(fit.optimality)
        fit.x=initial+transformation@fit.x
        fit.jac=model.jacobian(fit.x)
        fit.grad=fit.jac.T@fit.fun
        fit.optimality=float(np.linalg.norm(fit.grad,np.inf))
    else:
        fit=least_squares(fit_residual,initial,**options)
        scaled_optimality=float(fit.optimality)
    state=model.unpack(fit.x);old,_,details=model.old_blocks(fit.x,state);babar,predictions=model.babar_blocks(fit.x,state)
    column_norms=np.linalg.norm(fit.jac,axis=0)
    normalized_jacobian=fit.jac/column_norms
    left,singular,right=np.linalg.svd(normalized_jacobian,full_matrices=False)
    covariance=((right.T/np.maximum(singular,1e-14)**2)@right)/column_norms[:,None]/column_norms[None,:]
    boundaries=(('odd',0,12),('masses_and_leptonic',12,19),('shared_widths',19,26),('phi_KK',26,28),('BaBar',28,73),('normalization_constraints',73,75))
    output=dict(parameter_names=PARAMETER_NAMES,parameter_values=fit.x.tolist(),local_errors=np.sqrt(np.diag(covariance)).tolist(),
        original_odd_representative=odd_representative(fit.x[12:24]),physical_even_parameters=dict(zip(EVEN_NAMES,map(float,state[0]))),
        bibtex_keys=['Alda:2026nrd','Ruiz-Femenia:2003jdx','ParticleDataGroup:2024cfk','BaBar:2007ceh','Qin:2024ulb'],
        input_sha256={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in ['data/odd_parity/odd_low_energy_measurements.json','data/vector_sector/pdg2024.json','data/vector_sector/gv_measurements.json','data/vector_sector/alda_2026nrd_arxiv_v1.json',*(f'data/vector_sector/babar_2007_{channel}_inclusive_below_1p78.csv' for channel in CHANNELS)]},
        optimizer_method=args.method,svd_optimizer_coordinates=args.svd_scale,scaled_optimality=scaled_optimality,optimizer_ftol=args.ftol,initialization=args.initialization,measurement_count=73,constraint_count=2,parameter_count=32,chi2=float(fit.fun@fit.fun),dof=43,
        p_value=float(chi2.sf(fit.fun@fit.fun,43)),model_adequacy='inadequate' if chi2.sf(fit.fun@fit.fun,43)<.01 else 'not_rejected_by_nominal_chi2',global_minimum_proven=False,coefficient_confidence_intervals_certified=False,success=bool(fit.success),message=fit.message,nfev=int(fit.nfev),njev=int(fit.njev),optimality=float(fit.optimality),gradient_infinity_norm=float(np.linalg.norm(fit.grad,np.inf)),column_normalized_gradient_infinity_norm=float(np.max(abs(fit.grad/column_norms))),gauss_newton_predicted_chi2_reduction=float(np.linalg.norm(left.T@fit.fun)**2),gradient_by_parameter=dict(zip(PARAMETER_NAMES,map(float,fit.grad))),active_bounds=fit.active_mask.tolist(),
        residuals=fit.fun.tolist(),blocks={name:dict(chi2=float(fit.fun[start:stop]@fit.fun[start:stop]),count=stop-start) for name,start,stop in boundaries},
        coefficient_map=exact_odd_map(),odd_projected_information=projected_rank(fit.jac),
        sector_information={label:int(np.linalg.matrix_rank(fit.jac[:,columns]-fit.jac[:,[j for j in range(32) if j not in columns]]@np.linalg.lstsq(fit.jac[:,[j for j in range(32) if j not in columns]],fit.jac[:,columns],rcond=1e-12)[0],tol=1e-6)) for label,columns in [('VJP',list(range(12,17))),('VVP',list(range(17,20))),('VPPP',list(range(20,24)))]},details=details,
        babar_predictions={key:value.tolist() for key,value in predictions.items()},quadrature=dict(energy=args.energy_order,pair_per_segment=args.pair_order,pair_segments=3,angle_per_segment=args.angle_order,angle_segments=2,method='phi s-pole mapping and both Kstar t/u-pole mappings'),
        profiles_status='Only local quadratic profiles supplied unless --profiles or --g-profiles is requested; no confidence intervals certified.',
        assumptions=['Negative single-trace source; near-ideal inverse-mass branch.', 'No coefficient or previous-fit chi-square priors.',
        'Kstar Kpi saturation; each of the two duplicate experimental widths counted once.',
        'PDG summaries diagonal apart from phi KK BR correlation -0.72; missing cross-covariances set to zero.',
        'Full rotated ground neutral propagator; excited multiplets ideal; averaged vertex masses and physical phase-space masses.',
        'Ground rho/Kstar reference running shapes rescaled to shared on-shell widths; Kstar propagator uses charge-mean width.',
        'Fixed source excited poles, regulator 2.5 GeV; no K2star contribution.'],local_quadratic_profiles={name:dict(variance=float(covariance[12+index,12+index]),zero_delta_chi2=float(fit.x[12+index]**2/covariance[12+index,12+index])) for index,name in enumerate(ODD_NAMES)},profiles={})
    output['shared_even_gv_diagnostic']=even_gv_residual_diagnostic(model,fit.x)
    if args.audit_quadrature:
        audit_orders={dimension:getattr(args,f'audit_{dimension}_orders') or [order,2*order,4*order]
            for dimension,order in zip(('energy','pair','angle'),(args.energy_order,args.pair_order,args.angle_order))}
        if any(sequence[0]!=order for sequence,order in zip(audit_orders.values(),
                (args.energy_order,args.pair_order,args.angle_order))):
            parser.error('audit order sequences must begin at the corresponding fit order')
        output['quadrature_refinement']=audit_quadrature(model,fit.x,audit_orders)
    destination=ROOT/args.output;destination.parent.mkdir(parents=True,exist_ok=True)
    destination.write_text(json.dumps(output,indent=2)+'\n')
    if args.profiles or args.g_profiles:
        for index,name in enumerate(ODD_NAMES):
            if args.g_profiles and not args.profiles and index<8:continue
            fixed_index=12+index;profile=[]
            fixed_values=[0.]
            if args.profiles and index>=8:fixed_values.extend(fit.x[fixed_index]+displacement*max(np.sqrt(covariance[fixed_index,fixed_index]),1e-3) for displacement in (-1.,1.))
            for fixed_value in fixed_values:
                free=np.delete(np.arange(32),fixed_index)
                def profiled_residual(free_values):
                    values=fit.x.copy();values[free]=free_values;values[fixed_index]=fixed_value
                    return model.residual(values)
                def profiled_jacobian(free_values):
                    values=fit.x.copy();values[free]=free_values;values[fixed_index]=fixed_value
                    return model.jacobian(values)[:,free]
                result=least_squares(profiled_residual,fit.x[free],jac=profiled_jacobian,bounds=(lower[free],upper[free]),x_scale='jac',max_nfev=args.profile_max_nfev,ftol=1e-7)
                profile.append(dict(value=float(fixed_value),delta_chi2=float(result.fun@result.fun-fit.fun@fit.fun),success=bool(result.success),nfev=int(result.nfev)))
            output['profiles'][name]=profile;destination.write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps(dict(path=str(destination.relative_to(ROOT)),chi2=output['chi2'],success=output['success'],rank=output['odd_projected_information'])))


if __name__=='__main__':main()
