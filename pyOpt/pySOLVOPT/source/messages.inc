c This is an include file to compile with solvopt.f(or)
       character*(*) errmes,wrnmes,error2,error32,error42,error43
       character*(*) error52,error62,error63,error5,error6
       character*(*) warn1,warn20,warn21,warn4,warn31,warn32
       character*(*) warn09,warn08,termwarn0,termwarn1,appwarn
       character*(*) endwarn1, endwarn2, endwarn3, endwarn4
       parameter (errmes='SolvOpt error:') 
       parameter (wrnmes='SolvOpt warning:')
       parameter (error2='Improper space dimension.')
       parameter (error32='Function equals infinity at the point.')
       parameter (
     * error42='Gradient equals infinity at the starting point.',
     * error43='Gradient equals zero at the starting point.')
       parameter (
     * error52='<FUNC> returns infinite value at the point.',
     * error62='<GRADC> returns infinite vector at the point.',
     * error63='<GRADC> returns zero vector at an infeasible point.')
       parameter (
     * error5='Function is unbounded.',
     * error6='Choose another starting point.')
       parameter (warn1=
     * 'Gradient is zero, but stopping criteria are not fulfilled.',
     * warn20='Normal re-setting of a transformation matrix.',
     * warn21=
     * 'Re-setting due to the use of a new penalty coefficient.')
       parameter (
     * warn4='Iterations limit exceeded.',
     * warn31='The function is flat in certain directions.',
     * warn32='Trying to recover by shifting insensitive variables.',
     * warn09='Re-run from recorded point.',
     * warn08='Ravine with a flat bottom is detected.')
       parameter (
     * termwarn0='SolvOpt: Normal termination.',
     * termwarn1='SolvOpt: Termination warning:',
     * appwarn='The above warning may be reasoned'//
     *         ' by inaccurate gradient approximation',
     * endwarn1='Premature stop is possible.'//
     *    ' Try to re-run the routine from the obtained point.',
     * endwarn2='Result may not provide the optimum.'//
     *    ' The function apparently has many extremum points.',
     * endwarn3='Result may be inaccurate in the coordinates.'//
     *    ' The function is flat at the solution.',
     * endwarn4='Stopping criteria are not fulfilled.'//
     *    ' The function is very steep at the solution.')
