/* Routine for evaluating population members  */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "nsga2.h"
# include "rand.h"

/* Routine to evaluate objective function values and constraints for a population */
int evaluate_pop (population *pop, Global global)
{
    int i, res;
    for (i=0; i<global.popsize; i++)
    {
        res = evaluate_ind (&(pop->ind[i]),global);
        if (res)
        {
            return res;
        }
    }
    return 0;
}

/* Routine to evaluate objective function values and constraints for an individual */
int evaluate_ind (individual *ind, Global global)
{
    int j, res;
    res = nsga2func (global.nreal, global.nbin, global.nobj, global.ncon, ind->xreal, ind->xbin, ind->gene, ind->obj, ind->constr);
    if (res)
    {
        return res;
    }
    if (global.ncon==0)
    {
        ind->constr_violation = 0.0;
    }
    else
    {
        ind->constr_violation = 0.0;
        for (j=0; j<global.ncon; j++)
        {
            if (ind->constr[j]<0.0)
            {
                ind->constr_violation += ind->constr[j];
            }
        }
    }
    return res;
}
