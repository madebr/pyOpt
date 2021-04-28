# pyOpt Quick Reference Guide
Copyright (c) 2008-2014, pyOpt Developers

This is a quick guide to begin solving optimization problems with pyOpt.


## Optimization Problem Definition

pyOpt is design to solve general constrained nonlinear optimization problems:

```
    min  f(x)
     x

    s.t. g_j(x)  = 0, j = 1, ..., m_e
         g_j(x) <= 0, j = m_e + 1, ..., m

        x_i_L <= x_i <= x_i_U, i = 1, ..., n
```

  where:

  * `x` is the vector of design variables
  * `f(x)` is a nonlinear function
  * `g(x)` is a linear or nonlinear function
  * `n` is the number of design variables
  * `m_e` is the number of equality constraints
  * `m` is the total number of constraints (number of equality constraints: `m_i = m - m_e`)


## Optimization Class

Instantiating an Optimization Problem:

```python
opt_prob = Optimization(name, obj_fun, var_set={}, obj_set={}, con_set={})
```

where:
-   `name`: name of the problem (e.g. `'name'`)
-   `obj_fun`: objective function
-   `var_set`: dict containing the variables
-   `obj_set`: dict containing the objectives
-   `con_set`: dict containing the constraints

### Objective Function Template

```python
def obj_fun(x, *args, **kwargs):
    f = <your_function>(x, *args, **kwargs)
    g = <your_function>(x, *args, **kwargs)
    fail = 0

    return f, g, fail
```

where:

*   `f`: objective value
*   `g`: list of constraint values
    -   If the optimization problem is unconstrained,
        `g` must be an empty list: `g = []`
    -   Inequality constraints are handled as `<=`.
*   `fail`:
    -   `0`: successful function evaluation
    -   `1`: unsuccessful function evaluation (test must be provided by user)

### Adding an Objective

```python
opt_prob.addObj('name', value=0.0, optimum=0.0)
```

### Adding Design Variables

#### Single Design variable

```python
opt_prob.addVar('name', type='c', value=0.0, lower=-inf, upper=inf, choices=listofchoices)
```

#### Group of Design Variables

```python
opt_prob.addVarGroup('name', numerinGroup, type='c', value=value, lower=lb, upper=up,choices=listochoices)
```

where:

*   `value`, `lb`, `ub`: float, int or list
*   Supported types:
    -   `'c'`: continous design variable
    -   `'i'`: integer design variable
    -   `'d'`: discrete design variable (based on choices, e.g.: list/dict of materials)


### Adding Constraints

#### Single Constraint

```python
opt_prob.addCon('name', type='i', lower=-inf, upper=inf, equal=0.0)
```

#### A Group of Constraints

```python
opt_prob.addConGroup('name', numberinGroup, type='i', lower=lb, upper=up, equal=eq)
```

where:

*   `lb`, `ub`, `eq`: float, int or list
*   Supported Types:
	-   `'i'`: inequality constraint.
	-   `'e'`: equality constraint.


## Optimizer Class

Instanciating an Optimizer (e.g.: Snopt):

```python
opt = pySNOPT.SNOPT()
```

### Setting Optimizer Options

During instantiation:

```python
opt = pySNOPT.SNOPT(options={'name':value,...})
```

or one by one:
```python
opt.setOption('name',value)
```


### Getting Optimizer Options/Attributes

```python
opt.getOption('name')

opt.ListAttributes()
```

## Optimizing

### Solving the Optimization Problem

```python
opt(opt_prob, sens_type='FD', disp_opts=False, sens_mode='',*args, **kwargs)
```
where:

*   `sens_type`: sensitivity type
    - `'FD'`: finite differences
    - `'CS'`: complex step
*   `disp_opts`: flag for displaying the options in the solution output
*   `opt_prob`: user provided function
    -   format: `grad_function = lambda x, f, g : g_obj, g_con, fail` (See Objective FUnction Template section)
*   `sens_mode`: parallel sensitivity flag (`'serial'`,`'pgc'`, `'parallel'`)
*   Additional arguments and keyword arguments (e.g.: parameters) can be passed to the objective function

## Output

### Console output

*   Print Optimization problem with the initial values:
```python
print(opt_prob)
```

* Print specific solution of the Optimization problem:
```python
print(opt_prob._solutions[key])
```
where:
*   `key`: index in order of optimizer call.

### File output

```python
opt_prob.write2file(outfile='', disp_sols=False, solutions=[])
```

where:

*   `outfile`: filename or file instance (default name=`opt_prob name[0].txt`)
*   `disp_sols`: True will display all the stored solutions
*   `solutions`: list of indices of stored solutions to display

### Output as Input

The solution can be used directly as a optimization problem for
refinement by the same or a new optimizer:

```python
optimizer(opt_prob._solutions[key])
```
where:
*   `key`: index in order of optimizer call

The new solution will be stored as a sub-solution of the previous solution:

```python
print(opt_prob._solutions[key]._solutions[nkey])
```

## History and Hot Start

The history flag stores all function evaluations from an optimizer in
binary format in a .bin and .cue file:

```python
optimizer(opt_prob, store_hst=True)
```
where:
*   `store_hst`:
    -   `True`: use default file name for the history
    -   string type: custom filename

The binary history file can be used to hot start the optimizer if the
optimization was interrupted. The flag needs the filename of the
history (True will use the default name)

```python
optimizer(opt_prob, store_hst=True, hot_start=True)
```

If the store history flag is set as the same as the hot start flag a
temporary file will be created during the run and teh original file
will be overwritten at the end.

For hot start to work properly all options must be the same as when
the history was created.
