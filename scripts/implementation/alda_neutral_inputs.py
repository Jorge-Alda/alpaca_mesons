"""ArXiv-v1 Alda neutral generators from its full 3x3 mass diagonalization.

Canonical fields are converted by TP=F lambda_P/(sqrt(2) fpi).
Select physical or tree masses explicitly; physical eta is a phenomenological
kinematic replacement of the LO eigenmass, not a new diagonalization fit.
"""
from pathlib import Path
import json
import numpy as np

SNAPSHOT = Path(__file__).resolve().parents[2] / 'data/vector_sector/alda_arxiv_v1.json'


def alda_neutral_inputs(*, mass_scheme, F=.093, fpi=.093):
    """Return pi0/eta/etap dictionaries with TP, chi, F, mp and tree_mass.

    mass_scheme is 'tree' or 'physical'; masses in GeV, chi in GeV squared.
    Real diagonal matrices are immutable. No independent state rescaling.
    """
    if mass_scheme not in ('tree', 'physical'):
        raise ValueError("mass_scheme must be 'tree' or 'physical'")
    for label,value in (('F',F),('fpi',fpi)):
        if isinstance(value,bool) or not np.isrealobj(value) or not np.isscalar(value):
            raise ValueError(label+' must be a positive finite scalar')
        if not np.isfinite(value) or value <= 0:
            raise ValueError(label+' must be positive and finite')
    data=json.loads(SNAPSHOT.read_text())
    u,d,s=data['B0mq_gev2']; m0=data['m0_squared_gev2']
    mass=np.array([[u+d,(u-d)/np.sqrt(3),np.sqrt(2/3)*(u-d)],
                   [(u-d)/np.sqrt(3),(u+d+4*s)/3,np.sqrt(2)*(u+d-2*s)/3],
                   [np.sqrt(2/3)*(u-d),np.sqrt(2)*(u+d-2*s)/3,2*(u+d+s)/3+m0]])
    squared,rotation=np.linalg.eigh(mass)
    if np.any(squared <= 0):
        raise ValueError('neutral mass matrix must be positive definite')
    original=np.array([[1.,-1.,0.],[1.,1.,-2.]/np.sqrt(3),np.sqrt(2/3)*np.ones(3)])
    chi=np.diag(2*np.array([u,d,s]));chi.setflags(write=False)
    result={}
    for i,name in enumerate(('pi0','eta','etap')):
        diagonal=rotation[:,i]@original
        if diagonal@np.array(data['published_lambda_diagonals'][name])<0:
            rotation[:,i]*=-1; diagonal*=-1
        TP=np.diag(float(F)*diagonal/(np.sqrt(2)*float(fpi)));TP.setflags(write=False)
        result[name]=dict(TP=TP,chi=chi,F=float(F),mp=float(np.sqrt(squared[i]))
            if mass_scheme=='tree' else data['physical_masses_gev'][name],
            tree_mass=float(np.sqrt(squared[i])),mass_scheme=mass_scheme)
    return result
