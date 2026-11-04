"""Certify saved KKpi fit points without refitting; GeV and nb conventions.

Run: .venv/bin/python scripts/certify_kkpi_quadrature.py --output data/vector_sector/kkpi_quadrature_certification.json
Repeat: .venv/bin/python scripts/certify_kkpi_quadrature.py --replay data/vector_sector/kkpi_quadrature_certification.json --output /tmp/kkpi_repeated.json
"""
import argparse
import hashlib
import json
from pathlib import Path
import shlex
import subprocess
import sys

import numpy as np
import scipy

sys.path.insert(0,str(Path(__file__).resolve().parent))
import fit_babar_kkpi_global_vector_sector as fit

ROOT=fit.ROOT
REFERENCE='data/vector_sector/babar_kkpi_global_vector_sector_fit.json'


def file_hash(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_fit(path):
    record=json.loads(Path(path).read_text())
    if record['parameter_names']!=list(fit.PARAMETER_NAMES):
        raise ValueError(f'{path}: different fitted coordinates')
    if len(record['parameter_values'])!=32:
        raise ValueError(f'{path}: incomplete fitted point')
    return record


def valid_point(model,values):
    try:
        model.unpack(values)
        residual=model.residual(values)
        return bool(np.all(np.isfinite(residual)) and np.max(np.abs(residual))<1e8)
    except (ValueError,FloatingPointError,np.linalg.LinAlgError):
        return False


def select_points(model,reference,reference_path,competing=()):
    """Use signed strongest/weakest normalized-Jacobian directions, no refit.

    Unit-norm columns remove coordinate units. A displacement of radius 0.25
    is reduced by deterministic halves until it remains on the physical branch.
    These are numerical coverage points, not confidence intervals.
    """
    values=np.asarray(reference['parameter_values'],dtype=float)
    points=[dict(label='reference_candidate',parameter_values=values.tolist(),
        source=str(reference_path),rationale='saved reference candidate; no refit',
        source_success=reference.get('success'),source_chi2=reference.get('chi2'))]
    for path,record in sorted(competing,key=lambda item:str(item[0])):
        if record.get('input_sha256')!=reference.get('input_sha256'):
            raise ValueError(f'{path}: competing fit has different likelihood inputs')
        points.append(dict(label=f'competing_candidate_{len(points)}',
            parameter_values=record['parameter_values'],source=str(path),
            rationale='additional saved fit with identical coordinates and likelihood inputs',
            source_success=record.get('success'),source_chi2=record.get('chi2')))
    if not valid_point(model,values):
        return points,dict(basis='unavailable',reason='reference point is invalid',
            competing_sources=[str(path) for path,_ in competing])
    jacobian=model.jacobian(values)
    norms=np.linalg.norm(jacobian,axis=0)
    scales=1./np.maximum(norms,1e-12)
    _,singular,right=np.linalg.svd(jacobian*scales[None,:],full_matrices=False)
    for label,index in (('strong',0),('weak',len(singular)-1)):
        direction=right[index].copy()
        if direction[np.argmax(abs(direction))]<0:direction*=-1
        for sign in (-1,1):
            radius=.25
            for reduction in range(41):
                displaced=values+sign*radius*scales*direction
                if valid_point(model,displaced):break
                radius/=2
            else:
                raise ValueError(f'could not construct a valid {label} displacement')
            points.append(dict(label=f'{label}_{"minus" if sign<0 else "plus"}',
                parameter_values=displaced.tolist(),source=str(reference_path),
                rationale=f'{label} singular direction of column-normalized reference residual Jacobian; numerical coverage only',
                singular_value=float(singular[index]),signed_radius=sign*radius,
                validity_halvings=reduction))
    return points,dict(basis='column-normalized reference residual Jacobian',
        coordinate_scales=scales.tolist(),singular_values=singular.tolist(),
        initial_radius=.25,maximum_validity_halvings=40,
        covariance_used=False,confidence_interpretation=False,
        competing_sources=[str(path) for path,_ in competing],
        competing_search='same global-fit coordinate names and input hashes; progress snapshots excluded')


def source_metadata(reference_path,reference,competing_paths):
    paths=set(reference.get('input_sha256',{}))|{str(reference_path),*map(str,competing_paths)}
    hashes={path:file_hash(ROOT/path) for path in sorted(paths)}
    scripts=sorted({'scripts/certify_kkpi_quadrature.py'}|{str(Path(module.__file__).resolve().relative_to(ROOT))
        for module in tuple(sys.modules.values()) if getattr(module,'__file__',None)
        and Path(module.__file__).resolve().is_relative_to(ROOT)
        and not Path(module.__file__).resolve().is_relative_to(ROOT/'.venv')
        and Path(module.__file__).suffix=='.py'})
    revision=subprocess.run(['git','rev-parse','HEAD'],cwd=ROOT,text=True,capture_output=True,check=True).stdout.strip()
    return dict(input_sha256=hashes,reference_recorded_input_sha256=reference.get('input_sha256',{}),
        software=dict(git_revision=revision,working_tree_changes=subprocess.run(
            ['git','status','--porcelain'],cwd=ROOT,text=True,capture_output=True,check=True).stdout.splitlines(),
            source_sha256={path:file_hash(ROOT/path) for path in scripts},
            python=sys.version.split()[0],numpy=np.__version__,scipy=scipy.__version__))


def certify_points(model,points,orders):
    results=[]
    for point in points:
        values=np.asarray(point['parameter_values'],dtype=float)
        outcome=fit.audit_quadrature(model,values,orders)
        outcome.setdefault('derivative_settings',dict(
            method='analytic odd columns; forward differences for remaining columns',
            finite_difference_steps=(np.sqrt(np.finfo(float).eps)*np.maximum(1.,abs(values))).tolist(),
            parameter_scales=np.ones(32).tolist()))
        outcome.setdefault('tolerances',dict(fit.QUADRATURE_TOLERANCES))
        results.append(dict(point,**outcome))
    return results


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference',default=REFERENCE)
    parser.add_argument('--competing-fit',nargs='*',default=None)
    parser.add_argument('--replay',help='Repeat exact saved point vectors and integration orders')
    for dimension in ('energy','pair','angle'):
        parser.add_argument(f'--{dimension}-orders',nargs='+',type=int)
    parser.add_argument('--output',required=True)
    args=parser.parse_args(argv)
    command=shlex.join([sys.executable,'scripts/certify_kkpi_quadrature.py',*(sys.argv[1:] if argv is None else argv)])
    if args.replay:
        saved=json.loads((ROOT/args.replay).read_text())
        orders=saved['orders']
        if any(getattr(args,f'{name}_orders') is not None for name in orders):
            parser.error('replay uses recorded orders; do not override them')
        for path,expected in saved['provenance']['input_sha256'].items():
            if file_hash(ROOT/path)!=expected:raise ValueError(f'input hash changed: {path}')
        for path,expected in saved['provenance']['software']['source_sha256'].items():
            if file_hash(ROOT/path)!=expected:raise ValueError(f'software hash changed: {path}')
        points=saved['selected_points'];selection=saved['point_selection']
        provenance=saved['provenance']
    else:
        reference=read_fit(ROOT/args.reference)
        recorded=reference['quadrature']
        base=dict(energy=recorded['energy'],pair=recorded['pair_per_segment'],angle=recorded['angle_per_segment'])
        orders={name:getattr(args,f'{name}_orders') or [order,2*order,4*order] for name,order in base.items()}
        if any(orders[name][0]!=base[name] for name in base):
            parser.error('order sequences must start at the saved fit orders')
        model=fit.GlobalLikelihood(*(base[name] for name in ('energy','pair','angle')))
        competing_paths=(args.competing_fit if args.competing_fit is not None else
            [str(path.relative_to(ROOT)) for path in sorted((ROOT/'data/vector_sector').glob('babar_kkpi_global_vector_sector_fit*.json'))
                if path.resolve()!=(ROOT/args.reference).resolve() and not path.name.endswith('.progress.json')])
        competing=[(path,read_fit(ROOT/path)) for path in competing_paths]
        points,selection=select_points(model,reference,args.reference,competing)
        provenance=source_metadata(args.reference,reference,competing_paths)
    model=fit.GlobalLikelihood(*(orders[name][0] for name in ('energy','pair','angle')))
    results=certify_points(model,points,orders)
    counts={name:sum(point['status']==name for point in results) for name in ('certified','uncertified','unevaluated')}
    output=dict(schema_version=1,command=command,working_directory=str(ROOT),
        parameter_names=list(fit.PARAMETER_NAMES),selected_points=points,point_selection=selection,
        orders=orders,tolerances=fit.QUADRATURE_TOLERANCES,
        pole_mapping='phi s-pole mapping and both Kstar t/u-pole mappings',
        derivative_settings=dict(method='analytic odd columns; forward differences for remaining columns',
            step_rule='sqrt(machine epsilon) * max(1, abs(fitted coordinate)); fixed per point',
            parameter_scales=np.ones(32).tolist()),
        provenance=provenance,points=results,counts=counts,
        status='certified' if counts['certified']==len(points) else 'uncertified',
        scope='numerical integration at the listed points; no optimizer or fit-quality certification')
    if args.replay:output['replayed_from']=args.replay
    destination=ROOT/args.output;destination.parent.mkdir(parents=True,exist_ok=True)
    destination.write_text(json.dumps(output,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':main()
