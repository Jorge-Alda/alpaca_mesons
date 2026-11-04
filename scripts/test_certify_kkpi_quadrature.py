"""Point selection and saved-invocation reproduction tests.

Run: .venv/bin/python -m unittest scripts.test_certify_kkpi_quadrature
"""
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
from scripts import certify_kkpi_quadrature as certification


class Model:
    def __init__(self,*orders,**options):pass
    def unpack(self,values):
        if np.max(abs(values))>1:raise ValueError('outside physical branch')
    def residual(self,values):return np.zeros(75)
    def jacobian(self,values):return np.vstack((np.diag(np.arange(1.,33.)),np.zeros((43,32))))


def reference():
    return dict(parameter_names=list(certification.fit.PARAMETER_NAMES),parameter_values=[0.]*32,
        success=False,chi2=10.,input_sha256={},quadrature=dict(energy=1,pair_per_segment=2,angle_per_segment=2))


class PointSelectionTests(unittest.TestCase):
    def test_repeated_selection_has_identical_vectors_and_provenance(self):
        first=certification.select_points(Model(),reference(),'reference.json')
        second=certification.select_points(Model(),reference(),'reference.json')
        self.assertEqual(first,second)
        points,selection=first
        self.assertEqual(len(points),5)
        self.assertEqual([point['label'] for point in points],
            ['reference_candidate','strong_minus','strong_plus','weak_minus','weak_plus'])
        self.assertTrue(all(len(point['parameter_values'])==32 and point['rationale'] for point in points))
        self.assertTrue(all(certification.valid_point(Model(),np.array(point['parameter_values'])) for point in points))
        self.assertFalse(selection['confidence_interpretation'])

    def test_competing_source_included_and_different_likelihood_rejected(self):
        competing=reference();competing['success']=True
        points,_=certification.select_points(Model(),reference(),'reference.json',[('other.json',competing)])
        self.assertEqual(points[1]['source'],'other.json')
        self.assertTrue(points[1]['source_success'])
        competing['input_sha256']={'different':'hash'}
        with self.assertRaisesRegex(ValueError,'different likelihood inputs'):
            certification.select_points(Model(),reference(),'reference.json',[('other.json',competing)])

    def test_invalid_reference_remains_in_selected_point_list(self):
        record=reference();record['parameter_values'][0]=2.
        points,selection=certification.select_points(Model(),record,'invalid.json')
        self.assertEqual(len(points),1)
        self.assertEqual(points[0]['parameter_values'][0],2.)
        self.assertEqual(selection['basis'],'unavailable')

    def test_replay_rejects_changed_input_or_software_hash(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            (root/'input.json').write_text('changed')
            (root/'calculation.py').write_text('changed')
            for category in ('input','software'):
                saved=dict(orders=dict(energy=[1,2],pair=[2,3],angle=[2,3]),
                    provenance=dict(input_sha256={'input.json':'outdated'} if category=='input' else {},
                        software=dict(source_sha256={'calculation.py':'outdated'} if category=='software' else {})))
                (root/'saved.json').write_text(json.dumps(saved))
                with patch.object(certification,'ROOT',root),self.assertRaisesRegex(ValueError,f'{category} hash changed'):
                    certification.main(['--replay','saved.json','--output','repeated.json'])

    def test_high_order_run_reuses_historical_vectors_and_selection(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            (root/'reference.json').write_text(json.dumps(reference()))
            points,selection=certification.select_points(Model(),reference(),'reference.json')
            historical=dict(parameter_names=list(certification.fit.PARAMETER_NAMES),
                selected_points=points,point_selection=selection,
                provenance=dict(input_sha256={}))
            (root/'historical.json').write_text(json.dumps(historical))
            provenance=dict(input_sha256={},software=dict(source_sha256={}))
            with patch.object(certification,'ROOT',root),patch.object(certification.fit,'GlobalLikelihood',Model),\
                    patch.object(certification,'source_metadata',return_value=provenance),\
                    patch.object(certification,'select_points',side_effect=AssertionError('must reuse saved points')),\
                    patch.object(certification.fit,'audit_quadrature',return_value=dict(status='certified',passed=True)):
                certification.main(['--reference','reference.json','--competing-fit',
                    '--points-from','historical.json','--energy-orders','2','4',
                    '--pair-orders','64','128','--angle-orders','64','128','--output','higher.json'])
            result=json.loads((root/'higher.json').read_text())
            self.assertEqual(result['selected_points'],historical['selected_points'])
            self.assertEqual(result['point_selection'],historical['point_selection'])
            self.assertEqual(result['orders'],dict(energy=[2,4],pair=[64,128],angle=[64,128]))
            expected=certification.file_hash(root/'historical.json')
            self.assertEqual(result['selected_points_source'],dict(path='historical.json',sha256=expected))
            self.assertEqual(result['provenance']['input_sha256']['historical.json'],expected)
            self.assertEqual(result['quadrature_evaluation']['bin_block_size'],1)

    def test_saved_command_replays_every_selected_point_and_reports_invalid_points(self):
        def audit(model,values,orders):
            return dict(status='unevaluated' if values[0]>1 else 'certified',passed=bool(values[0]<=1),
                reason='outside physical branch' if values[0]>1 else '',orders=orders,refinements={})
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            (root/'reference.json').write_text(json.dumps(reference()))
            provenance=dict(input_sha256={},software=dict(source_sha256={}))
            with patch.object(certification,'ROOT',root),patch.object(certification.fit,'GlobalLikelihood',Model),\
                    patch.object(certification,'source_metadata',return_value=provenance),\
                    patch.object(certification.fit,'audit_quadrature',side_effect=audit):
                certification.main(['--reference','reference.json','--competing-fit','--energy-orders','1','2',
                    '--pair-orders','2','3','--angle-orders','2','3','--output','first.json'])
                first=json.loads((root/'first.json').read_text())
                self.assertEqual(len(first['points']),5)
                # Replay the recorded invocation, changing only its output destination.
                import shlex
                arguments=shlex.split(first['command'])[2:]
                arguments[arguments.index('--output')+1]='repeated.json'
                certification.main(arguments)
                repeated=json.loads((root/'repeated.json').read_text())
                self.assertEqual(first['selected_points'],repeated['selected_points'])
                self.assertEqual(first['points'],repeated['points'])
                first['selected_points'].append(dict(label='invalid',parameter_values=[2.]+[0.]*31,source='explicit',rationale='invalid test'))
                (root/'saved.json').write_text(json.dumps(first))
                certification.main(['--replay','saved.json','--output','replay.json'])
                replay=json.loads((root/'replay.json').read_text())
                self.assertEqual(len(replay['points']),6)
                self.assertEqual(replay['counts']['unevaluated'],1)
                self.assertEqual(replay['status'],'uncertified')
                self.assertEqual(replay['points'][-1]['reason'],'outside physical branch')
                self.assertEqual(len(replay['points'][-1]['derivative_settings']['finite_difference_steps']),32)


if __name__=='__main__':unittest.main()
