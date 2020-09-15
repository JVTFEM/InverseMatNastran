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
def exp_fname(fname):
    import pandas as pd
    import pathlib
    
    # Get the path of the current directory
    Path1 = str(pathlib.Path.cwd())
    
    # Get the name of all ".csv" files in current directory
    csv_s = list(pathlib.Path(Path1).glob('*.csv'))  
    
    if csv_s == []:
    # If no ".csv" files can be found use the default data file
        print("NO csv FILE FOUND USING DEFAULT FILE")
        fname = "exp_data.dat"
        
    else:
    # If ".csv" file(s) are found, use the first avalaible file and 
    # generate the experimental data file according to requested file name
        
        csv_file_name = str(csv_s[0]).replace(str(Path1)+'\\','')
        file_data = pd.read_csv(csv_file_name)
        new_file  = file_data.to_csv(header=None,index=False)
        f = open(fname, "w",newline='')
        f.write(str(new_file))
        f.close()
        
    return fname
# -------------------------------------------------------------------
# BDF, OP2 and experimental data files
# -------------------------------------------------------------------
BDF_FILE = 'platewithhole_FEM.dat'
EXP_FILE = exp_fname('exp_data.dat')   #(ID, x_loc, y_loc, z_loc, Tx, Ty, Tz)


# -------------------------------------------------------------------
# Evaluation function that is coupled to DOT
# Input:  x     - The design variable values from dot
#         obj   - The objective function value to change in place
#         g     - The constraint array to change in place
#         param - Possible pass through parameters - in this case the
#                 the transformation matrix
# -------------------------------------------------------------------
def get_units(fname):
    cur_file  = open(fname,"r")
    
    #Determine the total number of characters in file
    cur_file.seek(0,2)
    no_chars = cur_file.tell()
    
    #Set cursor back to starting position
    cur_file.seek(0,0)
    cur_pos = 0

    while cur_pos <= no_chars:

        line = cur_file.readline()
        cur_pos = cur_file.tell()

        if "$*                UNITS:" in line:
            
            #Look for presense of "mm" in units line
            if "mm" in line:
                units = "mm"
                cur_file.close()
                break
            
            #Look for presense of "Meter" in units line
            if "Meter" in line:
                units = "M"
                cur_file.close()
                break
                
    #if the units could not be found the default unit is assumed
    #to be in mm
    units = "mm"            
    cur_file.close()
    
    return units

units = get_units(BDF_FILE)

def myEvaluate(x, obj, g, param):

    # We unscale the design variables here.  We work with E and G 
    # since they have roughly the same order of magnitude (better for
    # the optimizer)
    if units == 'mm':
        # mN/mm^2 
        scalar = 1.e7   
    else:
        # N/m^2
        scalar = 1.e10   

    E = x[0] * scalar
    G = x[1] * scalar

    # Change the GENESIS input file
    anal.changeMAT1Card( BDF_FILE, 1, E, G )
    
    # Run GENESIS
    anal.runGENESIS( BDF_FILE, unZipOP2=False )
    
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
#
#-- Setup a transformation matrix if necessary
#femPnts = np.array( [[0., 0.],[1., 0.],[0., 1.],[1., 1.]] )
#dicPnts = np.array( [[1., 1.],[3., 1.],[1., 4.],[3., 4.]] )
#transMat  = setupTransform(femPnts, dicPnts)
#
#-- For the current example no transformation matrix is required
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
# Poision's ratio set to 0.5 as default value, therefore G = E/3
x[1] = x[0]/3 

# Instantiate the DOT object
aDot = dot.dot()

# Set the DOT parameters
aDot.nPrint   = 3
aDot.nMethod  = 1
aDot.evaluate = myEvaluate
aDot.nmParam  = transMat    # Simply pass through the transformationn matrix as a parameter

# Call DOT to perform the optimization
aDot.dotcall(x, xl, xu, nCons)
