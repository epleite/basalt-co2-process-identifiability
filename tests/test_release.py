from pathlib import Path
import unittest
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'data' / 'frozen'

class ReleaseTests(unittest.TestCase):
    def test_required_core_files(self):
        required = [
            ROOT / 'README.md',
            ROOT / 'CITATION.cff',
            ROOT / '.zenodo.json',
            ROOT / 'code' / 'run_all.py',
            ROOT / 'code' / 'reproduce_v1_3_fem.py',
            ROOT / 'code' / 'reproduce_v1_4.py',
            DATA / 'basalt_co2_v0_6_reflectivity_energy.csv',
            DATA / 'basalt_co2_v0_7_basis_similarity.csv',
            DATA / 'basalt_co2_v1_0_nominal_resistivity.csv',
            DATA / 'basalt_co2_v1_1_nominal_SIP_parameters.csv',
        ]
        for path in required:
            self.assertTrue(path.exists(), str(path))

    def test_v06_cancellation(self):
        d = pd.read_csv(DATA / 'basalt_co2_v0_6_reflectivity_energy.csv')
        r = d[(d.model_form == 'DEM-tangent') & (d.wave_solver == 'exact-multiples') &
              (d.time_days == 1825) & (d.acquisition == 'Geophone VSP')].iloc[0]
        self.assertAlmostEqual(r.cancellation_energy_ratio, 0.022448, places=5)

    def test_v07_elastic_antiparallel(self):
        d = pd.read_csv(DATA / 'basalt_co2_v0_7_basis_similarity.csv')
        r = d[(d.model_form == 'DEM-tangent') & (d.time_days == 1825) &
              (d.channel == 'PP+PS normalized') &
              (d.process_A == 'dissolution') & (d.process_B == 'cementation')].iloc[0]
        self.assertGreater(r.absolute_cosine_similarity, 0.99)
        self.assertLess(r.signed_cosine_similarity, -0.99)

    def test_v13_target_optimized(self):
        d = pd.read_csv(DATA / 'basalt_co2_v1_3_FEM_target_optimized.csv')
        r = d[(d.time_days == 1825) & (d.phase_noise_mrad == 5)].iloc[0]
        self.assertAlmostEqual(r.phase_RMS_D_mrad, 23.5366, places=3)
        self.assertAlmostEqual(r.phase_RMS_C_mrad, 6.1280, places=3)
        self.assertAlmostEqual(r.joint_abs_cos, 0.7070, places=3)

    def test_v14_null_direction(self):
        d = pd.read_csv(DATA / 'basalt_co2_v1_4_null_direction_variance.csv').set_index('data_case')
        ratio = d.loc['seismic + DC', 'variance_DplusC'] / d.loc['seismic + DC + SIP', 'variance_DplusC']
        self.assertGreater(ratio, 4.9)
        self.assertLess(ratio, 5.1)

if __name__ == '__main__':
    unittest.main()
