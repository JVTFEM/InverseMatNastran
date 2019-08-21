# -------------------------------------------------------------------
# This script will use DOT to minimize the difference between a
# measured (typically from DIC) 2D displacement field and that from
# a corresponding FE model.  The the RMSE is minimized.
#
# Gerhard Venter
# 29 July 2019
# -------------------------------------------------------------------
import dot as dot
import numpy as nm
import analysis as anal     # Library of utility functions for use
                            # with GENESIS linear static analysis
                            # using the pyNastran interface

# -------------------------------------------------------------------
# BDF, OP2 and experimental data files
# -------------------------------------------------------------------
BDF_FILE = 'platewithhole_FEM.dat'
EXP_FILE = 'exp_data.dat'   #(ID, x_loc, y_loc, z_loc, Tx, Ty, Tz)


# -------------------------------------------------------------------
# Evaluation function that is coupled to DOT
# Input:  x     - The design variable values from dot
#         obj   - The objective function value to change in place
#         g     - The constraint array to change in place
#         param - Possible pass through parameters - in this case the
#                 the transformation matrix
# -------------------------------------------------------------------
def myEvaluate(x, obj, g, param):

    # We unscale the design variables here.  We work with E and G 
    # since they have roughly the same order of magnitude (better for
    # the optimizer)
    E = x[0] * 1.e10
    G = x[1] * 1.e10

    # Change the GENESIS input file
    anal.changeMAT1Card( BDF_FILE, 1, E, G )
    
    # Run GENESIS
    anal.runGENESIS( BDF_FILE, unZipOP2=True )
    
    # Get the objective function value.  We multiply it with a large 
    # value to get optimum values in the range of 1 or so (better for
    # the optimizer)
    obj.value = 1.e8*anal.getObjectiveFn( BDF_FILE, EXP_FILE, 1 , 
                                         transMat=param).sum()
    
    return

# -------------------------------------------------------------------
# The main part of the script - setup and call DOT to do the 
# optimization
# -------------------------------------------------------------------

# Start by setting up a transformation matrix to transform points
# from the DIC to the FEM coordinate systems.  Do only once for each
# sample
#femPnts = np.array( [[0., 0.],[1., 0.],[0., 1.],[1., 1.]] )
#dicPnts = np.array( [[1., 1.],[3., 1.],[1., 4.],[3., 4.]] )
#transMat  = setupTransform(femPnts, dicPnts)
transMat = None

# Now start setting up the optimization environment
nDvar = 2  # We have two design variables
nCons = 0  # We have not constraints

# Setup the design variable initial, lower and upper bounds
x  = nm.empty(nDvar, float)
xl = nm.empty(nDvar, float)
xu = nm.empty(nDvar, float)

for i in range(nDvar):
    xl[i] = 2.
    xu[i] = 20.
    x[i]  = 10.

# Instantiate the DOT object
aDot = dot.dot()

# Set the DOT parameters
aDot.nPrint = 3
aDot.nMethod = 1
aDot.evaluate = myEvaluate
aDot.nmParam  = transMat    # Simply pass through the transformatoin matrix

# Call DOT to perform the optimization
aDot.dotcall(x, xl, xu, nCons)
