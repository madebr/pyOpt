def test_rosenbrock():
    # !/usr/bin/env python
    """Solves Rosenbrock's Unconstrained Problem.

    min     100*(x2-x1^2)**2 + (1-x1)^2
    s.t.:   -10 <= xi <= 10,  i = 1,2

    f* = 0 , x* = [1, 1]
    """

    from pyOpt import Optimization

    def objfunc(x):
        f = 100 * (x[1] - x[0] ** 2) ** 2 + (1 - x[0]) ** 2
        g = []

        fail = 0
        return f, g, fail

    def getlastsolution(prob: Optimization):
        new_index = prob.firstavailableindex(prob.getSolSet())
        return prob.getSol(new_index - 1)

    # Instantiate Optimization Problem
    opt_prob = Optimization('Rosenbrock Unconstraint Problem', objfunc)
    opt_prob.addVar('x1', 'c', lower=-10.0, upper=10.0, value=-3.0)
    opt_prob.addVar('x2', 'c', lower=-10.0, upper=10.0, value=-4.0)
    opt_prob.addObj('f')

    # Instantiate Optimizer (PSQP) & Solve Problem
    from pyOpt import PSQP

    psqp = PSQP()
    psqp.setOption('IPRINT', 0)
    psqp(opt_prob, sens_type='FD')

    # Instantiate Optimizer (SLSQP) & Solve Problem
    from pyOpt import SLSQP

    slsqp = SLSQP()
    slsqp.setOption('IPRINT', -1)
    slsqp(opt_prob, sens_type='FD')

    # Instantiate Optimizer (CONMIN) & Solve Problem
    from pyOpt import CONMIN

    conmin = CONMIN()
    conmin.setOption('IPRINT', 0)
    conmin(opt_prob, sens_type='CS')

    # Instantiate Optimizer (COBYLA) & Solve Problem
    from pyOpt import COBYLA

    cobyla = COBYLA()
    cobyla.setOption('IPRINT', 0)
    cobyla(opt_prob)

    # Instantiate Optimizer (SOLVOPT) & Solve Problem
    from pyOpt import SOLVOPT

    solvopt = SOLVOPT()
    solvopt.setOption('iprint', -1)
    solvopt(opt_prob, sens_type='FD')

    # Instantiate Optimizer (KSOPT) & Solve Problem
    from pyOpt import KSOPT

    ksopt = KSOPT()
    ksopt.setOption('IPRINT', 0)
    ksopt(opt_prob, sens_type='FD')

    # Instantiate Optimizer (NSGA2) & Solve Problem
    # TODO: reactivate, currently NSGA2 fails (when building the wheels and then testing)
    # from pyOpt import NSGA2
    #
    # nsga2 = NSGA2()
    # nsga2.setOption('PrintOut', 0)
    # nsga2(opt_prob)

    # Instantiate Optimizer (SDPEN) & Solve Problem
    from pyOpt import SDPEN

    sdpen = SDPEN()
    sdpen.setOption('iprint', -1)
    sdpen(opt_prob)
