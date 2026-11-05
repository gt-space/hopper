"""Bound nonlinear work per attempt so hard timesteps retry promptly."""
from fullflow import Transient
from fullflow.Solvers.transient.operations import TransientStepSolve


class BoundedStepSolve(TransientStepSolve):
    def __init__(self, *args, max_evaluations=20, **kwargs):
        if not isinstance(max_evaluations, int) or isinstance(max_evaluations, bool) or max_evaluations < 1:
            raise ValueError('max_evaluations must be a positive integer.')
        super().__init__(*args, **kwargs)
        self.max_evaluations = max_evaluations

    def _least_squares_kwargs(self, cache, x0, settings):
        kwargs = super()._least_squares_kwargs(cache, x0, settings)
        kwargs['max_nfev'] = self.max_evaluations
        return kwargs


def transient_solver(network, max_evaluations=20):
    solver = Transient(network)
    solver.step_solver = BoundedStepSolve(network, solver._cache, solver.evaluator,
                                         max_evaluations=max_evaluations)
    return solver
