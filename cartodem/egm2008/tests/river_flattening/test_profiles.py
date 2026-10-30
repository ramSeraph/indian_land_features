import unittest

import numpy as np

from river_flattening.profiles import fit_network


class ProfileTests(unittest.TestCase):
    def test_shared_confluence_and_braided_reaches(self):
        # Two tributaries merge at 2, split below 3, then join again at 6.
        reaches = [[0, 1, 2, 3], [4, 5, 2], [3, 7, 6], [3, 8, 6]]
        distances = [np.arange(len(r))*120.0 for r in reaches]
        edges = np.array([pair for r in reaches for pair in zip(r[:-1], r[1:])])
        truth = np.array([12, 10, 8, 6, 12, 10, 2, 4, 4], dtype=float)
        ids = np.repeat(np.arange(len(truth)), 5)
        observed = np.repeat(truth, 5) + np.tile([0, .1, .2, .4, 25], len(truth))
        fitted, report = fit_network(len(truth), edges, ids, observed,
                                      reaches=reaches, distances=distances,
                                      quantile=.2, stiffness=.2)
        self.assertTrue(np.all(fitted[edges[:, 0]] >= fitted[edges[:, 1]]))
        np.testing.assert_allclose(fitted, truth, atol=.25)
        self.assertEqual(fitted.dtype, np.float32)
        self.assertEqual(report['solver_status'], 'solved')

    def test_noisy_profile_does_not_propagate_one_low_outlier(self):
        n = 31
        x = np.arange(n)*120.0
        truth = 40 - np.arange(n)*.4
        ids = np.repeat(np.arange(n), 8)
        samples = np.repeat(truth, 8) + np.tile([-.1, 0, .1, .2, .3, .4, 10, 30], n)
        samples[11*8] = -80
        edges = np.column_stack([np.arange(n-1), np.arange(1, n)])
        fitted, _ = fit_network(n, edges, ids, samples, reaches=[list(range(n))],
                                distances=[x], stiffness=1)
        self.assertLess(float(np.max(np.abs(fitted-truth))), .5)
        self.assertTrue(np.all(np.diff(fitted) <= 0))

    def test_missing_observations_are_allowed_only_in_supported_components(self):
        fitted, _ = fit_network(3, [[0, 1], [1, 2]], [0, 2], [3, 1],
                                reaches=[[0, 1, 2]], distances=[[0, 100, 200]])
        np.testing.assert_allclose(fitted, [3, 2, 1], atol=.01)
        with self.assertRaisesRegex(ValueError, 'no valid bank'):
            fit_network(4, [[0, 1], [2, 3]], [0, 1], [3, 1])

    def test_cycle_and_reversed_reach_are_rejected(self):
        with self.assertRaisesRegex(ValueError, 'cycle'):
            fit_network(3, [[0, 1], [1, 2], [2, 0]], [0, 1, 2], [3, 2, 1])
        with self.assertRaisesRegex(ValueError, 'disagrees'):
            fit_network(3, [[0, 1], [1, 2]], [0, 1, 2], [3, 2, 1],
                        reaches=[[2, 1, 0]], distances=[[0, 100, 200]])

    def test_bank_bounds_are_respected_without_clamping_to_zero(self):
        fitted, _ = fit_network(3, [[0, 1], [1, 2]], [0, 1, 2], [2, 1, 0],
                                reaches=[[0, 1, 2]], distances=[[0, 100, 200]],
                                upper_bounds=[-.5, -1, -1.5])
        self.assertTrue(np.all(fitted <= np.array([-.5, -1, -1.5])+.0001))
        self.assertTrue(np.all(fitted < 0))


if __name__ == '__main__':
    unittest.main()
