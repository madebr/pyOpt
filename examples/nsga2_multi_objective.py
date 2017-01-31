from pyOpt import Optimization, NSGA2

def objfunc(xdict):
    x = xdict['x']
    y = xdict['y']

    ff = 2*[0.0]
    ff[0] = (x - 0.0)**2 + (y - 0.0)**2
    ff[1] = (x - 1.0)**2 + (y - 1.0)**2
    gg = []
    fail = False

    return ff, gg, fail

# Instantiate Optimization Problem
optProb = Optimization('Rosenbrock function', objfunc)
optProb.addVar('x', 'c', value=0, lower=-600, upper=600)
optProb.addVar('y', 'c', value=0, lower=-600, upper=600)

optProb.addObj('obj1')
optProb.addObj('obj2')

#300 generations will find x=(0,0), 200 or less will find x=(1,1)
options = {
    'maxGen':200
}
opt = NSGA2(options=options)
sol = opt(optProb)

print(sol)
