"""
pyIPOPT - A Python pyOpt interface to IPOPT.

Copyright (c) 2011 by Sylvain Arreckx
All rights reserved. Not to be used for commercial purposes.
Revision: 0.1   $Date: 23/03/2011 12:00$


Developers:
-----------
- Mr. Sylvain Arreckx (SA)
- Mr. Stefan Bethge (SB)


History
-------
        v. 0.2   - Python 3 compatibility (SB, 2017)
        v. 0.1  - Initial Class Creation (SA, 2011)
"""

__version__ = '$Revision: $'

# =============================================================================
# IPOPT Library
# =============================================================================
try:
    import pyipopt
except ImportError:
    print('Error: IPOPT shared library failed to import')

import copy
# =============================================================================
# Standard Python modules
# =============================================================================
import os
import time

# =============================================================================
# External Python modules
# =============================================================================
import numpy

# =============================================================================
# Extension modules
# =============================================================================
from pyOpt import Gradient, History, Optimizer

# =============================================================================
# Misc Definitions
# =============================================================================
inf = 10.E+20  # define a value for infinity
# =============================================================================
eps = 1.0  # define a value for machine precision
while (eps / 2.0 + 1.0) > 1.0:
    eps = eps / 2.0

eps = 2.0 * eps

# =============================================================================
# IPOPT Optimizer Class
# =============================================================================


class IPOPT(Optimizer):
    """IPOPT Optimizer Class - Inherited from Optimizer Abstract Class."""

    def __init__(self, pll_type=None, *args, **kwargs):
        """IPOPT Optimizer Class Initialization.

        Documentation last updated:  Feb. 16, 2010 - Peter W. Jansen
        """

        name = 'IPOPT'
        category = 'Local Optimizer'
        def_opts = {
            # IPOPT Printing Options
            # Print Control (0 - None, 1 - Final,2,3,4,5 - Debug)
            'IPRINT': [int, 2],
            'IOUT': [int, 6],  # Output Unit Number
            'output_file': [str, 'IPOPT.out'],  # Output File Name

            # Output options
            'print_level': [int, 5],  # Output verbosity level
            # Print all options set by the user
            'print_user_options': [str, 'no'],
            # Switch to print all algorithmic options
            'print_options_documentation': [str, 'no'],
            'file_print_level': [int, 5],  # Verbosity level for output file
            'option_file_name': [str, ''],  # File name of options file

            # Termination options
            'tol': [float, 1e-8],  # relative convergence tolerance
            'max_iter': [int, 3000],  # Maximum number of iterations
            'max_cpu_time': [float, 1e+6],  # Maximum number of CPU seconds.
            # Desired threshold for the dual infeasibility
            'dual_inf_tol': [float, 1],
            # Desired threshold for the constraint violation
            'constr_viol_tol': [float, 1e-4],
            # Desired threshold for the complementarity conditions
            'compl_inf_tol': [float, 1e-4],
            # "Acceptable" convergence tolerance (relative)
            'acceptable_tol': [float, 1e-6],
            # Number of "acceptable" iterates before triggering termination
            'acceptable_iter': [int, 15],
            # "Acceptance" threshold for the constraint violation
            'acceptable_constr_viol_tol': [float, 1e-2],
            # "Acceptance" threshold for the dual infeasibility.
            'acceptable_dual_inf_tol': [float, 1e+10],
            # "Acceptance" threshold for the complementarity conditions
            'acceptable_compl_inf_tol': [float, 1e-2],
            # "Acceptance" stopping criterion based on objective function change
            'acceptable_obj_change_tol': [float, 1e+20],
            # Threshold for maximal value of primal iterates
            'diverging_iterates_tol': [float, 1e+20],

            # NLP scaling options
            # Scaling factor for the objective function
            'obj_scaling_factor': [float, 1],
            # Select the technique use for scaling the NLP ('none',
            # 'user-scaling', 'gradient-based', 'equilibration-based')
            'nlp_scaling_method': [str, 'gradient-based'],
            # Maximum gradient after NLP scaling
            'nlp_scaling_max_gradient': [float, 100],
            # Minimum value of gradient-based scaling values
            'nlp_scaling_min_value': [float, 1e-8],

            # NLP options
            # Factor for initial relaxation of the bounds
            'bound_relax_factor': [float, 1e-8],
            # Indicates whether final points should be projected into original
            # bounds
            'honor_original_bounds': [str, 'yes'],
            # Indicates whether it is desired to check for Nan/Inf in derivative
            # matrices
            'check_derivatives_for_naninf': [str, 'no'],
            # any bound less or equal this value will be considered -inf (i.e.
            # not lower bounded)
            'nlp_lower_bound_inf': [float, -1e+19],
            # any bound greater or this value will be considered +inf (i.e. not
            # upper bounded)
            'nlp_upper_bound_inf': [float, 1e+19],
            # Determines how fixed variables should be handled
            # ('make_parameter', 'make_constraint', 'relax_bounds')
            'fixed_variable_treatment': [str, 'make_parameter'],
            # Indicates whether all equality constraints are linear
            'jac_c_constant': [str, 'no'],
            # Indicates whether all inequality constraints are linear
            'jac_d_constant': [str, 'no'],
            # Indicates whether the problem is a quadratic problem
            'hessian_constant': [str, 'no'],
            # Initialization options
            # Desired minimum relative distance from the initial point to bound
            'bound_frac': [float, 0.01],
            # Desired minimum absolute distance from the initial point to bound
            'bound_push': [float, 0.01],
            # Desired minimum relative distance from the initial slack to bound
            'slack_bound_frac': [float, 0.01],
            # Desired minimum absolute distance from the initial slack to bound
            'slack_bound_push': [float, 0.01],
            # Initial value for the bound multipliers
            'bound_mult_init_val': [float, 1],
            # Maximum allowed least-square guess of constraint multipliers
            'constr_mult_init_max': [float, 1000],
            # Initialization method for bound multipliers ('constant',
            # 'mu_based')
            'bound_mult_init_method': [str, 'constant'],

            # Barrier parameter options
            # Indicates if we want to do Mehrotra's algorithm
            'mehrotra_algorithm': [str, 'no'],
            # Update strategy for barrier parameter ('monotone', 'adaptive')
            'mu_strategy': [str, 'monotone'],
            # Oracle for a new barrier parameter in the adaptive strategy
            # ('probing', 'loqo', 'quality-function')
            'mu_oracle': [str, 'quality-function'],
            # Maximum number of search steps during direct search procedure
            # determining the optimal centering parameter
            'quality_function_max_section_steps': [int, 8],
            # Oracle for the barrier parameter when switching to fixed mode
            # ('probing', 'loqo', 'quality-function', 'average_compl')
            'fixed_mu_oracle': [str, 'average_compl'],
            'mu_init': [float, 0.1],  # Initial value for the barrier parameter
            # Factor for initialization of maximum value for barrier parameter
            'mu_max_fact': [float, 1000],
            'mu_max': [float, 1e+5],  # Maximum value for barrier parameter
            'mu_min': [float, 1e-11],  # Minimum value for barrier parameter
            'mu_target': [float, 0],  # Desired value of complementarity
            # Factor for mu in barrier stop test.
            'barrier_tol_factor': [float, 10],
            # Determines linear decrease rate of barrier parameter
            'mu_linear_decrease_factor': [float, 0.2],
            # Determines superlinear decrease rate of barrier parameter (between
            # 1 and 2)
            'mu_superlinear_decrease_power': [float, 1.5],

            # Multiplier updates
            # Method to determine the step size for constraint multiplier
            # ('primal', 'bound-mult', 'min', 'max', 'full', 'min-dual-infeas',
            # 'safer-min-dual-infeas', 'primal-and-full', 'dual-and-full',
            # 'acceptor')
            'alpha_for_y': [str, 'primal'],
            # Tolerance for switching to full equality multiplier steps
            'alpha_for_y_tol': [float, 10],
            # Tells the algorithm to recalculate the equality and inequality
            # multipliers as least square estimates
            'recalc_y': [str, 'no'],
            # Feasibility threshold for recomputation of multipliers
            'recalc_y_feas_tol': [float, 1e-6],

            # Line search options
            # Maximum number of second order correction trial steps at each
            # iteration
            'max_soc': [int, 4],
            # Number of shortened iterations that trigger the watchdog
            'watchdog_shortened_iter_trigger': [int, 10],
            # Maximum number of watchdog iterations
            'watchdog_trial_iter_max': [int, 3],
            # Always accept the first trial step
            'accept_every_trial_step': [str, 'no'],
            # The type of corrector steps that should be taken (unsupported!)
            # ('none', 'affine', 'primal-dual')
            'corrector_type': [str, 'none'],
            'dependency_detector': [str, 'none'],
            # Which linear solver should be used to detect linearly dependent equality constraints (experimental)
            # ('ma28, 'mumps', 'none')

            # Warm start options
            # Warm-start for initial point
            'warm_start_init_point': [str, 'no'],
            # same as bound_push for the regular initializer
            'warm_start_bound_push': [float, 0.001],
            # same as bound_frac for the regular initializer
            'warm_start_bound_frac': [float, 0.001],
            # same as slack_bound_frac for the regular initializer
            'warm_start_slack_bound_frac': [float, 0.001],
            # same as slack_bound_push for the regular initializer
            'warm_start_slack_bound_push': [float, 0.001],
            # same as mult_bound_push for the regular initializer.
            'warm_start_mult_bound_push': [float, 0.001],
            # Maximum initial value for the equality multipliers
            'warm_start_mult_init_max': [float, 1e+6],

            # Restoration Phase
            # Enable heuristics to quickly detect an infeasible problem
            'expect_infeasible_problem': [str, 'no'],
            # Threshold for disabling "expect_infeasible_problem" option
            'expect_infeasible_problem_ctol': [float, 0.001],
            # Multiplier threshold for activating "expect_infeasible_problem"
            # option
            'expect_infeasible_problem_ytol': [float, 1e+8],
            # Tells algorithm to switch to restoration phase in first iteration
            'start_with_resto': [str, 'no'],
            # Required reduction in primal-dual error in the soft restoration
            # phase
            'soft_resto_pderror_reduction_factor': [float, 0.9999],
            # Required reduction of infeasibility before leaving restoration
            # phase
            'required_infeasibility_reduction': [float, 0.9],
            # Threshold for resetting bound multipliers after the restoration
            # phase
            'bound_mult_reset_threshold': [float, 1000],
            # Threshold for resetting equality and inequality multipliers after
            # restoration phase
            'constr_mult_reset_threshold': [float, 0],
            # Determines if the original objective function should be evaluated
            # at restoration phase trial points
            'evaluate_orig_obj_at_resto_trial': [str, 'yes'],

            # Linear Solver
            'linear_solver': [str, 'ma57'],
            # Method for scaling the linear system ('none', 'mc19',
            # 'slack-based')
            'linear_system_scaling': [str, 'mc19'],
            # Flag indicating that linear scaling is only done if it seems
            # required
            'linear_scaling_on_demand': [str, 'yes'],
            # Maximum number of iterative refinement steps per linear system
            # solve
            'max_refinement_steps': [float, 10],
            # Minimum number of iterative refinement steps per linear system
            # solve
            'min_refinement_steps': [float, 1],

            # Hessian Perturbation
            # Maximum value of regularization parameter for handling negative
            # curvature
            'max_hessian_perturbation': [float, 1e+20],
            # Smallest perturbation of the Hessian block
            'min_hessian_perturbation': [float, 1e-20],
            # Size of first x-s perturbation tried
            'first_hessian_perturbation': [float, 0.0001],
            # Increase factor for x-s perturbation for very first perturbation
            'perturb_inc_fact_first': [float, 100],
            # Increase factor for x-s perturbation
            'perturb_inc_fact': [float, 8],
            # Decrease factor for x-s perturbation
            'perturb_dec_fact': [float, 0.333333],
            # Size of the regularization for rank-deficient constraint Jacobians
            'jacobian_regularization_value': [float, 1e-8],

            # Quasi-Newton
            # Indicates what Hessian information is to be used ('exact',
            # 'limited-memory')
            'hessian_approximation': [str, 'exact'],
            # Quasi-Newton update formula for the limited memory approximation
            # ('bfgs', 'sr1')
            'limited_memory_update_type': [str, 'bfgs'],
            # Maximum size of the history for the limited quasi-Newton Hessian
            # approximation
            'limited_memory_max_history': [int, 6],
            # Threshold for successive iterations where update is skipped
            'limited_memory_max_skipping': [int, 2],
            # Initialization strategy for the limited memory quasi-Newton
            # approximation ('scalar1', 'scalar2', 'scalar3', 'scalar4',
            # 'constant')
            'limited_memory_initialization': [str, 'scalar1'],
            # Value for B0 in low-rank update
            'limited_memory_init_val': [float, 1],
            # Upper bound on value for B0 in low-rank update
            'limited_memory_init_val_max': [float, 1e+8],
            # Lower bound on value for B0 in low-rank update
            'limited_memory_init_val_min': [float, 1e-8],
            # Determines if the quasi-Newton updates should be special during
            # the restoration phase
            'limited_memory_special_for_resto': [str, 'no'],

            # Derivative Test options
            # Enable derivative checker ('none', 'first-order', 'second-order',
            # 'only-second-order')
            'derivative_test': [str, 'none'],
            # Size of the finite difference perturbation in derivative test
            'derivative_test_perturbation': [float, 1e-8],
            # Threshold for indicating wrong derivative
            'derivative_test_tol': [float, 0.0001],
            # Indicates whether information for all estimated derivatives should
            # be printed
            'derivative_test_print_all': [str, 'no'],
            # Index of first quantity to be checked by derivative checker
            'derivative_test_first_index': [int, -2],
            # Maximal perturbation of an evaluation point
            'point_perturbation_radius': [float, 10],

            # MA57 Linear Solver options
            # Pivot tolerance for the linear solver MA57
            'ma57_pivtol': [float, 1e-8],
            # Maximum pivot tolerance for the linear solver MA57
            'ma57_pivtolmax': [float, 0.0001],
            # Safety factor for work space memory allocation for the linear
            # solver MA57
            'ma57_pre_alloc': [float, 1.05],
            'ma57_pivot_order': [int, 5],  # Controls pivot order in MA57
            # Controls MA57 automatic scaling
            'ma57_automatic_scaling': [str, 'yes'],
            # Controls block size used by Level 3 BLAS in MA57BD
            'ma57_block_size': [int, 16],
            'ma57_node_amalgamation': [int, 16],  # Node amalgamation parameter
            'ma57_small_pivot_flag': [int, 0],
        }
        informs = {
            0: 'Solve Succeeded',
            1: 'Solved To Acceptable Level',
            2: 'Infeasible Problem Detected',
            3: 'Search Direction Becomes Too Small',
            4: 'Diverging Iterates',
            5: 'User Requested Stop',
            6: 'Feasible Point Found',
            -1: 'Maximum Iterations Exceeded',
            -2: 'Restoration Failed',
            -3: 'Error In Step Computation',
            -4: 'Maximum CpuTime Exceeded',
            -10: 'Not Enough Degrees Of Freedom',
            -11: 'Invalid Problem Definition',
            -12: 'Invalid Option',
            -13: 'Invalid Number Detected',
            -100: 'Unrecoverable Exception',
            -101: 'NonIpopt Exception Thrown',
            -102: 'Insufficient Memory',
            -199: 'Internal Error'
        }
        self.set_options = []
        Optimizer.__init__(self, name, category, def_opts, informs, *args,
                           **kwargs)

    def __solve__(self,
                  opt_problem={},
                  sens_type='FD',
                  store_sol=True,
                  disp_opts=False,
                  store_hst=False,
                  hot_start=False,
                  sens_mode='',
                  sens_step={},
                  *args,
                  **kwargs):
        """Run Optimizer (Optimize Routine)

        **Keyword arguments:**

        - opt_problem -> INST: Optimization instance
        - sens_type -> STR/FUNC: Gradient type, *Default* = 'FD'
        - store_sol -> BOOL: Store solution in Optimization class flag,
          *Default* = True
        - disp_opts -> BOOL: Flag to display options in solution text, *Default*
          = False
        - store_hst -> BOOL/STR: Flag/filename to store optimization history,
          *Default* = False
        - hot_start -> BOOL/STR: Flag/filename to read optimization history,
          *Default* = False
        - sens_mode -> STR: Flag for parallel gradient calculation, *Default* =
          ''
        - sens_step -> FLOAT: Sensitivity setp size, *Default* = {} [corresponds
          to 1e-6 (FD), 1e-20(CS)]

        Documentation last updated:  Feb. 2, 2011 - Peter W. Jansen
        """

        self.pll = False
        self.myrank = 0

        myrank = self.myrank

        tmp_file = False
        def_fname = self.options['output_file'][1].split('.')[0]
        if isinstance(store_hst, str):
            if isinstance(hot_start, str):
                if (myrank == 0):
                    if (store_hst == hot_start):
                        hos_file = History(hot_start, 'r', self)
                        log_file = History(store_hst + '_tmp', 'w', self,
                                           opt_problem.name)
                        tmp_file = True
                    else:
                        hos_file = History(hot_start, 'r', self)
                        log_file = History(store_hst, 'w', self,
                                           opt_problem.name)

                self.sto_hst = True
                self.hot_start = True
            elif hot_start:
                if (myrank == 0):
                    hos_file = History(store_hst, 'r', self)
                    log_file = History(store_hst + '_tmp', 'w', self,
                                       opt_problem.name)
                    tmp_file = True

                self.sto_hst = True
                self.hot_start = True
            else:
                if (myrank == 0):
                    log_file = History(store_hst, 'w', self, opt_problem.name)

                self.sto_hst = True
                self.hot_start = False

        elif store_hst:
            if isinstance(hot_start, str):
                if (hot_start == def_fname):
                    if (myrank == 0):
                        hos_file = History(hot_start, 'r', self)
                        log_file = History(def_fname + '_tmp', 'w', self,
                                           opt_problem.name)
                        tmp_file = True

                else:
                    if (myrank == 0):
                        hos_file = History(hot_start, 'r', self)
                        log_file = History(def_fname, 'w', self,
                                           opt_problem.name)

                self.sto_hst = True
                self.hot_start = True
            elif hot_start:
                if (myrank == 0):
                    hos_file = History(def_fname, 'r', self)
                    log_file = History(def_fname + '_tmp', 'w', self,
                                       opt_problem.name)
                    tmp_file = True

                self.sto_hst = True
                self.hot_start = True
            else:
                if (myrank == 0):
                    log_file = History(def_fname, 'w', self, opt_problem.name)

                self.sto_hst = True
                self.hot_start = False

        else:
            self.sto_hst = False
            self.hot_start = False

        gradient = Gradient(opt_problem, sens_type, sens_mode, sens_step,
                            *args, **kwargs)

        def eval_f(x, user_data=None):
            """IPOPT - Objective Value Function."""
            # Variables Groups Handling
            if opt_problem.use_groups:
                xg = {}
                for group in group_ids.keys():
                    if (group_ids[group][1] - group_ids[group][0] == 1):
                        xg[group] = x[group_ids[group][0]]
                    else:
                        xg[group] = x[group_ids[group][0]:group_ids[group][1]]

                xn = xg
            else:
                xn = x

            # Flush Output Files
            self.flushFiles()

            # Evaluate User Function
            fail = 0
            # if (myrank == 0):
            #    if self.hot_start:
            #        [vals,hist_end] = hos_file.read(ident=['obj', 'con', 'fail'])
            #        if hist_end:
            #            self.hot_start = False
            #            hos_file.close()
            #        else:
            #            [ff,gg,fail] = [vals['obj'][0][0],vals['con'][0],int(vals['fail'][0][0])]
            #
            #

            # if self.pll:
            #    self.hot_start = Bcast(self.hot_start,root=0)

            # if self.hot_start and self.pll:
            #    [ff,gg,fail] = Bcast([ff,gg,fail],root=0)
            # else:
            [ff, gg, fail] = opt_problem.obj_fun(xn, *args, **kwargs)

            # Store History
            if (myrank == 0):
                if self.sto_hst:
                    log_file.write(x, 'x')
                    log_file.write(ff, 'obj')
                    log_file.write(gg, 'con')
                    log_file.write(fail, 'fail')

                # Objective Assigment
            if isinstance(ff, complex):
                f = ff.astype(float)
            else:
                f = ff

            # Constraints Assigment
            g = numpy.zeros(len(opt_problem._constraints.keys()))
            for i in range(len(opt_problem._constraints.keys())):
                if isinstance(gg[i], complex):
                    g[i] = gg[i].astype(float)
                else:
                    g[i] = gg[i]

            return f

        def eval_g(x, user_data=None):

            # Variables Groups Handling
            if opt_problem.use_groups:
                xg = {}
                for group in group_ids.keys():
                    if (group_ids[group][1] - group_ids[group][0] == 1):
                        xg[group] = x[group_ids[group][0]]
                    else:
                        xg[group] = x[group_ids[group][0]:group_ids[group][1]]

                xn = xg
            else:
                xn = x

            # Flush Output Files
            self.flushFiles()

            # Evaluate User Function
            fail = 0
            #            if (myrank == 0):
            #                if self.hot_start:
            #                    [vals,hist_end] = hos_file.read(ident=['obj', 'con', 'fail'])
            #                    if hist_end:
            #                        self.hot_start = False
            #                        hos_file.close()
            #                    else:
            #                        [ff,gg,fail] = [vals['obj'][0][0],vals['con'][0],int(vals['fail'][0][0])]

            # if self.pll:
            #   self.hot_start = Bcast(self.hot_start,root=0)

            # if self.hot_start and self.pll:
            #    [ff,gg,fail] = Bcast([ff,gg,fail],root=0)
            # else:
            [ff, gg, fail] = opt_problem.obj_fun(xn, *args, **kwargs)

            # Store History
            if (myrank == 0):
                if self.sto_hst:
                    log_file.write(x, 'x')
                    log_file.write(ff, 'obj')
                    log_file.write(gg, 'con')
                    log_file.write(fail, 'fail')

                # Objective Assigment
            if isinstance(ff, complex):
                f = ff.astype(float)
            else:
                f = ff

            # Constraints Assigment
            g = numpy.zeros(len(opt_problem._constraints.keys()))
            for i in range(len(opt_problem._constraints.keys())):
                if isinstance(gg[i], complex):
                    g[i] = gg[i].astype(float)
                else:
                    g[i] = gg[i]

            return g

        def eval_grad_f(x, user_data=None):
            """IPOPT - Objective/Constraint Gradients Function."""
            # if self.hot_start:
            #    if (myrank == 0):
            #        [vals,hist_end] = hos_file.read(ident=['grad_obj','grad_con'])
            #        if hist_end:
            #            self.hot_start = False
            #            hos_file.close()
            #        else:
            #            dff = vals['grad_obj'][0].reshape((len(opt_problem._objectives.keys()),len(opt_problem._variables.keys())))
            #            dgg = vals['grad_con'][0].reshape((len(opt_problem._constraints.keys()),len(opt_problem._variables.keys())))
            #
            #
            #    if self.pll:
            #        self.hot_start = Bcast(self.hot_start,root=0)
            #
            #    if self.hot_start and self.pll:
            #        [dff,dgg] = Bcast([dff,dgg],root=0)
            #

            # if not self.hot_start:

            [f, g, fail] = opt_problem.obj_fun(x, *args, **kwargs)
            dff, dgg = gradient.getGrad(x, group_ids, [f], g, *args, **kwargs)

            # Store History
            if self.sto_hst and (myrank == 0):
                log_file.write(dff, 'grad_obj')
                log_file.write(dgg, 'grad_con')

            # Gradient Assignment
            df = numpy.zeros(len(opt_problem._variables.keys()))

            for i in range(len(opt_problem._variables.keys())):
                df[i] = dff[0, i]

            return df

        def eval_grad_g(x, flag, user_data=None):

            # if self.hot_start:
            #    if (myrank == 0):
            #        [vals,hist_end] = hos_file.read(ident=['grad_obj','grad_con'])
            #        if hist_end:
            #            self.hot_start = False
            #            hos_file.close()
            #        else:
            #            dff = vals['grad_obj'][0].reshape((len(opt_problem._objectives.keys()),len(opt_problem._variables.keys())))
            #            dgg = vals['grad_con'][0].reshape((len(opt_problem._constraints.keys()),len(opt_problem._variables.keys())))
            #
            #
            #    if self.pll:
            #        self.hot_start = Bcast(self.hot_start,root=0)
            #
            #    if self.hot_start and self.pll:
            #        [dff,dgg] = Bcast([dff,dgg],root=0)
            #

            # if not self.hot_start:

            if flag:
                a = numpy.zeros(
                    len(opt_problem._variables.keys()) *
                    len(opt_problem._constraints.keys()), int)
                b = numpy.zeros(
                    len(opt_problem._variables.keys()) *
                    len(opt_problem._constraints.keys()), int)

                for i in range(len(opt_problem._constraints.keys())):
                    for j in range(len(opt_problem._variables.keys())):
                        a[i * len(opt_problem._variables.keys()) + j] = i
                        b[i * len(opt_problem._variables.keys()) + j] = j
                return (a, b)

            else:
                [f, g, fail] = opt_problem.obj_fun(x, *args, **kwargs)
                dff, dgg = gradient.getGrad(x, group_ids, [f], g, *args,
                                            **kwargs)

                # Store History
                if self.sto_hst and (myrank == 0):
                    log_file.write(dff, 'grad_obj')
                    log_file.write(dgg, 'grad_con')

                # Gradient Assignment
                a = numpy.zeros([
                    len(opt_problem._variables.keys()) *
                    len(opt_problem._constraints.keys())
                ])
                for i in range(len(opt_problem._constraints.keys())):
                    for j in range(len(opt_problem._variables.keys())):
                        a[i * len(opt_problem._variables.keys()) +
                          j] = dgg[i, j]

                return a

        # Variables Handling
        nvar = len(opt_problem._variables.keys())
        xl = []
        xu = []
        xx = []
        for key in opt_problem._variables.keys():
            if opt_problem._variables[key].type == 'c':
                xl.append(opt_problem._variables[key].lower)
                xu.append(opt_problem._variables[key].upper)
                xx.append(opt_problem._variables[key].value)
            elif opt_problem._variables[key].type == 'i':
                raise OSError('IPOPT cannot handle integer design variables')
            elif opt_problem._variables[key].type == 'd':
                raise OSError('IPOPT cannot handle discrete design variables')

        xl = numpy.array(xl)
        xu = numpy.array(xu)
        xx = numpy.array(xx)

        # Variables Groups Handling
        group_ids = {}
        if opt_problem.use_groups:
            k = 0
            for key in opt_problem._vargroups.keys():
                group_len = len(opt_problem._vargroups[key]['ids'])
                group_ids[opt_problem._vargroups[key][
                    'name']] = [k, k + group_len]
                k += group_len

            # Constraints Handling
        ncon = len(opt_problem._constraints.keys())
        blc = []
        buc = []
        if ncon > 0:
            for key in opt_problem._constraints.keys():
                if opt_problem._constraints[key].type == 'e':
                    blc.append(opt_problem._constraints[key].equal)
                    buc.append(opt_problem._constraints[key].equal)
                elif opt_problem._constraints[key].type == 'i':
                    blc.append(opt_problem._constraints[key].lower)
                    buc.append(opt_problem._constraints[key].upper)

        else:
            if ((store_sol) and (myrank == 0)):
                print("Optimization Problem Does Not Have Constraints\n")
                print("Unconstrained Optimization Initiated\n")

            ncon = 1
            blc.append(-inf)
            buc.append(inf)

        blc = numpy.array(blc)
        buc = numpy.array(buc)

        # Objective Handling
        objfunc = opt_problem.obj_fun
        nobj = len(opt_problem._objectives.keys())
        ff = []
        for key in opt_problem._objectives.keys():
            ff.append(opt_problem._objectives[key].value)

        ff = numpy.array(ff)

        # Create an IPOPT instance problem
        nnzj = nvar * ncon
        nnzh = nvar * nvar
        ipopt = pyipopt.create(nvar, xl, xu, ncon, blc, buc, nnzj, nnzh,
                               eval_f, eval_grad_f, eval_g, eval_grad_g)

        # Setup Options
        optionss = self.options.copy()
        del optionss['defaults']

        for i in optionss:
            if not self.options['defaults'][i][1] == optionss[i][1]:
                if self.options[i][0].__name__ == 'int':
                    ipopt.int_option(i, self.options[i][1])

                if self.options[i][0].__name__ == 'float':
                    ipopt.num_option(i, self.options[i][1])

                if self.options[i][0].__name__ == 'str':
                    ipopt.str_option(i, self.options[i][1])

        # Run IPOPT

        t0 = time.time()
        r = ipopt.solve(xx)
        sol_time = time.time() - t0

        if (myrank == 0):
            if self.sto_hst:
                log_file.close()
                if tmp_file:
                    hos_file.close()
                    name = hos_file.filename
                    os.remove(name + '.cue')
                    os.remove(name + '.bin')
                    os.rename(name + '_tmp.cue', name + '.cue')
                    os.rename(name + '_tmp.bin', name + '.bin')

        ipopt.close()

        # Store Results
        sol_inform = {}
        print(r)
        sol_inform['value'] = r[-1]  # ifail[0]
        sol_inform['text'] = self.getInform(r[-1])  # self.getInform(ifail[0])

        if store_sol:
            sol_name = 'IPOPT Solution to ' + opt_problem.name

            sol_options = copy.copy(self.options)
            if 'default' in sol_options:
                del sol_options['defaults']

            sol_evals = 0

            sol_vars = copy.deepcopy(opt_problem._variables)
            i = 0
            x = r[0]
            for key in sol_vars.keys():
                sol_vars[key].value = x[i]
                i += 1

            sol_objs = copy.deepcopy(opt_problem._objectives)
            sol_objs[0].value = r[4]

            sol_cons = {}

            if ncon > 0:
                sol_lambda = r[3]
            else:
                sol_lambda = {}

            opt_problem.addSol(
                self.__class__.__name__,
                sol_name,
                objfunc,
                sol_time,
                sol_evals,
                sol_inform,
                sol_vars,
                sol_objs,
                sol_cons,
                sol_options,
                display_opts=disp_opts,
                Lambda=sol_lambda,
                Sensitivities=sens_type,
                myrank=myrank,
                arguments=args,
                **kwargs)

        return ff, xx, sol_inform  # ifail[0]

    def _on_setOption(self, name, value):
        """Set Optimizer Option Value (Optimizer Specific Routine).

        Documentation last updated:  May. 07, 2008 - Ruben E. Perez
        """

        pass

    def _on_getOption(self, name):
        """Get Optimizer Option Value (Optimizer Specific Routine).

        Documentation last updated:  May. 07, 2008 - Ruben E. Perez
        """

        pass

    def _on_getInform(self, infocode):
        """Get Optimizer Result Information (Optimizer Specific Routine).

        Keyword arguments:
        -----------------
        id -> STRING: Option Name

        Documentation last updated:  May. 07, 2008 - Ruben E. Perez
        """

        try:
            inform_text = self.informs[infocode]
        except:
            inform_text = 'Unknown Exit Status'

        return inform_text

        pass

    def _on_flushFiles(self):
        """Flush the Output Files (Optimizer Specific Routine).

        Documentation last updated:  August. 09, 2009 - Ruben E. Perez
        """
        pass
