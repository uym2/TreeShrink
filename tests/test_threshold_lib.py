import math
import os
import sys
import unittest

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

from treeshrink.threshold_lib import threshold_l_kernel, threshold_loglnorm


KERNEL_VALUES = np.array([
    1.02, 1.05, 1.08, 1.11, 1.25, 1.32,
    1.41, 1.55, 1.76, 2.05, 2.40, 2.90,
])

LOGLNORM_VALUES = np.array([
    1.02, 1.04, 1.07, 1.12, 1.18, 1.31,
    1.52, 1.91, 2.65, 3.80, 5.40, 7.90,
])


class ThresholdLibTests(unittest.TestCase):
    def test_threshold_l_kernel_matches_validated_values(self):
        cases = [
            (0.10, 2.590107),
            (0.05, 2.975282),
            (0.01, 3.664773),
        ]
        for e, expected in cases:
            with self.subTest(e=e):
                self.assertEqual(threshold_l_kernel(KERNEL_VALUES, e=e), expected)

    def test_threshold_loglnorm_matches_validated_values(self):
        cases = [
            (0.10, 7.935737),
            (0.05, 37.467218),
            (0.01, 31077.847904),
        ]
        for e, expected in cases:
            with self.subTest(e=e):
                self.assertEqual(threshold_loglnorm(LOGLNORM_VALUES, e=e), expected)

    def test_threshold_l_kernel_ignores_non_positive_and_non_finite_values(self):
        noisy_values = np.concatenate((
            KERNEL_VALUES,
            np.array([-10.0, -1.0, 0.0, math.inf, -math.inf, math.nan]),
        ))

        self.assertEqual(
            threshold_l_kernel(noisy_values, e=0.05),
            threshold_l_kernel(KERNEL_VALUES, e=0.05),
        )

    def test_threshold_loglnorm_ignores_values_with_non_positive_or_non_finite_logs(self):
        noisy_values = np.concatenate((
            LOGLNORM_VALUES,
            np.array([0.25, 0.5, 1.0, math.inf, math.nan]),
        ))

        self.assertEqual(
            threshold_loglnorm(noisy_values, e=0.05),
            threshold_loglnorm(LOGLNORM_VALUES, e=0.05),
        )

    def test_threshold_l_kernel_requires_two_positive_finite_values(self):
        cases = [
            [],
            [0.0],
            [1.0],
            [math.nan, math.inf],
        ]
        for values in cases:
            with self.subTest(values=values):
                with self.assertRaisesRegex(ValueError, "Need at least 2 positive finite values"):
                    threshold_l_kernel(values)

    def test_threshold_loglnorm_requires_two_values_with_positive_finite_logs(self):
        cases = [
            [],
            [0.5, 1.0],
            [1.0, math.inf, math.nan],
        ]
        for values in cases:
            with self.subTest(values=values):
                with self.assertRaisesRegex(ValueError, "Need at least 2 values with positive finite logs"):
                    threshold_loglnorm(values)

    def test_threshold_l_kernel_rejects_constant_positive_values(self):
        with self.assertRaisesRegex(ValueError, "Data are constant"):
            threshold_l_kernel([2.0, 2.0, 2.0])


if __name__ == "__main__":
    unittest.main()
