from pyOpt import NSGA2, Optimization


def objfunc(xdict):
    x = xdict['x']
    y = xdict['y']

    ff = [
        (x - 0.0)**2 + (y - 0.0)**2,
        (x - 1.0)**2 + (y - 1.0)**2,
    ]
    gg = []
    fail = False

    return ff, gg, fail


# Instantiate Optimization Problem
optProb = Optimization('Rosenbrock function', objfunc, use_groups=True)
optProb.addVar('x', 'c', value=0, lower=-600, upper=600)
optProb.addVar('y', 'c', value=0, lower=-600, upper=600)

optProb.addObj('obj1')
optProb.addObj('obj2')

# 300 generations will find x=(0,0), 200 or less will find x=(1,1)
options = {
    'maxGen': 200,
}
opt = NSGA2(options=options)
opt.setOption('PrintOut', 0)
opt(optProb)

print(optProb.getSol(0))
