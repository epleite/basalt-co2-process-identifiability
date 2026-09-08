"""Consistency checks for the initial published tables, not full model reproduction."""
from pathlib import Path
import unittest
import pandas as pd

DATA = Path(__file__).resolve().parents[1] / 'data' / 'frozen'

class PublishedTableTests(unittest.TestCase):
    def test_null_direction_variance_ratio(self):
        table = pd.read_csv(DATA / 'basalt_co2_v1_4_null_direction_variance.csv').set_index('data_case')
        ratio = table.loc['seismic + DC', 'variance_DplusC'] / table.loc['seismic + DC + SIP', 'variance_DplusC']
        self.assertGreater(ratio, 4.9)
        self.assertLess(ratio, 5.1)

    def test_fem_target_optimized(self):
        table = pd.read_csv(DATA / 'basalt_co2_v1_3_FEM_target_optimized.csv')
        row = table[(table.time_days == 1825) & (table.phase_noise_mrad == 5)].iloc[0]
        self.assertAlmostEqual(row.phase_RMS_D_mrad, 23.5366, places=3)
        self.assertAlmostEqual(row.phase_RMS_C_mrad, 6.1280, places=3)
        self.assertAlmostEqual(row.joint_abs_cos, 0.7070, places=3)

    def test_posterior_variance_consistency(self):
        posterior = pd.read_csv(DATA / 'basalt_co2_v1_4_headline_posterior.csv').set_index('data_case')
        variances = pd.read_csv(DATA / 'basalt_co2_v1_4_null_direction_variance.csv').set_index('data_case')
        for name, row in posterior.iterrows():
            covariance = row.corr_D_C * row.sd_aD * row.sd_aC
            predicted = (row.sd_aD**2 + row.sd_aC**2 + 2 * covariance) / 2
            self.assertAlmostEqual(predicted, variances.loc[name, 'variance_DplusC'], places=12)

if __name__ == '__main__':
    unittest.main()
