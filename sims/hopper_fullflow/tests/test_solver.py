"""Ensure work limits preserve FullFlow's residual tolerance and step retry."""
import unittest
import numpy as np
from fullflow import Network
from fullflow.Solvers.steady_state.settings import LeastSquaresSettings
from Hopper_Solver import transient_solver


class SolverTests(unittest.TestCase):
    def test_work_budget_is_passed_to_scipy_without_relaxing_tolerances(self):
        class Cache:
            def lower_bound_array(self): return np.array([-np.inf])
            def upper_bound_array(self): return np.array([np.inf])
            def keep_feasible_array(self): return np.array([False])
        solver = transient_solver(Network('test'), max_evaluations=12)
        settings = LeastSquaresSettings(rtol=1e-6, ftol=1e-12, xtol=1e-12)
        kwargs = solver.step_solver._least_squares_kwargs(Cache(), np.array([1.]), settings)
        self.assertEqual(kwargs['max_nfev'], 12)
        self.assertEqual(kwargs['ftol'], 1e-12)
        self.assertEqual(kwargs['xtol'], 1e-12)
        self.assertEqual(settings.rtol, 1e-6)
        self.assertEqual(kwargs['method'], 'trf')

    def test_invalid_budget_rejected(self):
        for value in (0, -1, True, 1.5):
            with self.assertRaises(ValueError):
                transient_solver(Network('test'), max_evaluations=value)
