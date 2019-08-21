# -------------------------------------------------------------------
# A set of functions that can be used to change a MAT1 material card
# in a GENESIS input data deck, then run the GENESIS analysis and
# extract and compare the new FE displacement values (from the op2
# file) with a set of reference (experimental) data points.  The FE
# and expermental data points are not necessarily located at the
# same position so interpolation is used.
#
# These functions assume that a linear static case is considered with
# MAT1 material properties and that the X,Y displacements are only
# compared in the X,Y plane.
#
# Gerhard Venter
# 29 July 2019
# -------------------------------------------------------------------
from pyNastran.bdf.bdf import BDF
from pyNastran.op2.op2 import OP2
from scipy.interpolate import Rbf
import math as m
import numpy as np
import pandas as pd
import subprocess as sp
import statsmodels.api as sm
import matplotlib.pylab as plt


# -------------------------------------------------------------------
# GENESIS and gzip installations
# These values are machine specific and may need to be updated
# -------------------------------------------------------------------
GENESIS = '/opt/vrand/bin/genesis'
ZIP     = '/usr/bin/gzip'


# -------------------------------------------------------------------
# Function to create the GENESIS OP2 filename from the BDF file name
# using the standard GENESIS conventio
# Input: bdfFile    - The name of the BDF file
#        iterID     - The GENESIS design iteration to use (default=0)
#        compressed - Did GENESIS compress the file?
# Output: The corresponding GENESIS OP2 filename
# -------------------------------------------------------------------
def getGenesisOP2Filename(bdfFile, iterID=0, compressed=True):
    
    # Find the right most occurance of '.' - we want to split off the 
    # file extension
    idx = bdfFile.rindex('.')
    
    # Keep everything up to the last occurance of '.'
    baseName = bdfFile[:idx]
    
    # Add the Genesis formatting for the OP2 file
    baseName = baseName + '%02d.op2'%(iterID)
    
    # Add the compressed extension if needed
    if (compressed):
        baseName = baseName + '.gz'
    
    return baseName


# -------------------------------------------------------------------
# Function to change the BDF input file with new isotropic material 
# data.  We work with E and G because the are in the same order of
# magnitude which works much better for the optimizer.  The Poisson's
# ratio will be calculated on the fly from the E and G values.
#
# Input: bdfFile  - The name of the BDF file to open and modify (in 
#                   place)
#         matID   - The Nastran ID of the MAT1 card to modify
#         newE    - The new Youngs Modulus value
#         newG    - The new Shear Modulus value
# -------------------------------------------------------------------
def changeMAT1Card(bdfFile, matID, newE, newG):
    
    # Read the BDF file
    modelBDF = BDF(debug=False)
    modelBDF.read_bdf(bdfFile)    

    # Get the material property to change
    mat = modelBDF.materials.get(matID)
    
    # Change the material property data in place
    mat.e  = newE
    mat.g  = newG
    mat.nu = (newE/(2.*newG)) - 1.  # Calculate the new Poisson ratio
    
    # Write the file back out - this will overwrite the old file
    modelBDF.write_bdf(bdfFile)
    
    return


# -------------------------------------------------------------------
# Function to run GENESIS and wait for it to finish.  Can also unzip 
# the op2 file produced by GENESIS if wanted.
#
# Input: bdfFile  - The name of the BDF file to open and modify (in 
#                   place)
#        unZipOP2 - Flag do unzip the OP2 file produced by GENESIS
#                   (default is not to unzip)
# -------------------------------------------------------------------
def runGENESIS(bdfFile, unZipOP2=False):
    
    # Launc the sub-process to run GENESIS.  This will wait for the 
    # process to finish before continueing (thus running in the
    # foreground)
    sp.run( [GENESIS, bdfFile] )
    
    # If requested, unzip the OP2 file
    if ( unZipOP2 ):
        sp.run([ZIP, '-f', '-d', getGenesisOP2Filename(bdfFile,compressed=True)])
        
    return


# -------------------------------------------------------------------
# Function that does the heavy lifting of calculating the objective 
# function value by comparing X,Y displacements from experimental 
# data with X, Y displacements from FE data.  This function assumes
# a static analysis (thus time = 0)
#
# Input: bdfFile  - The name of the BDF file to open and modify (in 
#                   place)
#        expFile  - File containing experimental data as ID, X-loc,
#                   Y-loc, Z-loc, X-disp, Y-disp, Z-disp
#        loadCase - Load case to extract
#        transMat - Transformation matrix obtained from setupTransform
#                   method.  If None - no transformation is performed
#                   (Default = None)
# Output: (x_err, y_err) - The RMSE for the x and y components of 
#                          the displacement field
# -------------------------------------------------------------------
def getObjectiveFn(bdfFile, expFile, loadCase, transMat=None):
    
    # Static data, so only the first timestep will be read
    TIME = 0  
    
    # Read the BDF file
    modelBDF = BDF(debug=False)
    modelBDF.read_bdf(bdfFile)

    # Read the OP2 file
    op2File = getGenesisOP2Filename(bdfFile, compressed=False)
    modelOP2 = OP2(debug=False)
    modelOP2.read_op2(op2File)

    # Get all the nodal displacements (X,Y,Z) from the OP2 file for
    # the specified loadCase
    disp = modelOP2.displacements[loadCase]
    txyz = disp.data[TIME, :, :3]

    # Get all the nodal data from the BDF file - loaded node by 
    # node to ensure the same order as the disp data
    iCnt = 0;
    gxyz = np.zeros( (txyz.shape[0], 3) )
    for (nid, ntype) in disp.node_gridtype:
        gxyz[iCnt] = modelBDF.nodes.get(nid).xyz
        iCnt = iCnt + 1

    # Create RBF interpolation for X, Y location of data and disps
    rbf_x = Rbf(gxyz[:,0], gxyz[:,1], txyz[:,0], function='linear')
    rbf_y = Rbf(gxyz[:,0], gxyz[:,1], txyz[:,1], function='linear')

    # Read the experimental data as Pandas dataframe
    expData = pd.read_csv( expFile )

    # Scale the experimental data points if required
    if (transMat != None) :
        
        # First scale the location of the grid points
        dic_gxy = transformDICPnts(transMat, expData.iloc[:,[1,2]])
    
        # Then scale and recalculate the displacements
        # -- One needs to first add the displacements to the grid locations, 
        # -- then scale, then subtract the scale grid locations to get the
        # -- scaled displacements
        dic_txy = transformDICPnts(transMat, expData.iloc[:,[1,2]] + expData.iloc[:,[4,5]])
        dic_txy = dic_txy - dic_gxy

    # Else just use the points as is
    else:
        dic_gxy = expData.iloc[:,[1,2]]
        dic_txy = expData.iloc[:,[4,5]]
    
    # Predict at each of the experimetnal data points using the 
    # RBF functions
    x_pred = rbf_x( dic_gxy.iloc[:,0], dic_gxy.iloc[:,1] )
    y_pred = rbf_y( dic_gxy.iloc[:,0], dic_gxy.iloc[:,1] )

    # Calculate the RMSE difference between the predicted and 
    # actual data points
    x_err = x_pred - dic_txy.iloc[:,0]
    y_err = y_pred - dic_txy.iloc[:,1]
    x_err = m.sqrt(np.dot( x_err, x_err)/x_err.shape[0])
    y_err = m.sqrt(np.dot( y_err, y_err)/y_err.shape[0])
    
    # Some debug plotting that was used to check the data and 
    # RBF fits
    #plt.plot( x_pred, expData.iloc[:,4] )
    #plt.show()
    #plt.plot( y_pred, expData.iloc[:,5] )
    #plt.show()
    
    return np.array( (x_err, y_err) )


# -------------------------------------------------------------------
# Function to setup Affine transformations from the DIC coordinate
# system to the FEM coordinate system.  The user can supply any
# number of points, but there must be at least 3.  The X, Y coordinates
# of the points in the FEM and DIC coordinate systems are provided by
# the user.  A re-weighted least squares method is used to obtain the
# Affine transformation.
#
# Input:  femPnts  - A (nPoints, 2) numpy array with the X,Y coordinates
#                    of the points in the FEM coordinate system
#         dicPnts  - A (nPoints, 2) numpy array with the X,Y coordinates
#                    of the same points in the DIC coordinate system
# Output: transMat - A (2,3) transformation matrix to be used with 
#                    the trasformDICPnts function
# -------------------------------------------------------------------
def setupTransform(femPnts, dicPnts):
    
    # Check that the dimensions of the femPnts and dicPnts are the same
    if (femPnts.shape != dicPnts.shape):
        raise ValueError('FEM Pnts and DIC Pnts must have the same shape')
    
    # Detect how many points are supplied - must be AT LEAST 3 - more is better
    nPnts = femPnts.shape[0]
    if (nPnts < 3):
        raise ValueError('Must specify at least 3 Pnts to construct the transformation matrix')
    
    # Allocate space for the X matrix in the least squares
    X = np.ones( (nPnts, 3) )

    # Now setup the X matrix
    X[:,1:] = dicPnts

    # Now perform a re-weighted least squares - X coordinates
    modelRLM = sm.RLM( femPnts[:,0], X )
    resRLM   = modelRLM.fit()
    param_x  = resRLM.params
    
    # Now perform a re-weighted least squares - Y coordinates
    modelRLM = sm.RLM( femPnts[:,1], X )
    resRLM   = modelRLM.fit()
    param_y  = resRLM.params
    
    # Return the transformation matrices
    return np.array( [param_x, param_y] )


# -------------------------------------------------------------------
# Function to transfrom points from the DIC coordinate system to the
# FEM coordinate system.  This function makes use of a transformation
# matrix created by the setupTransform method.
#
# Input:  transMat  - A (2,3) transformation matrix that was obtained
#                     from the setupTransform method
#         dicPnts   - A (nPoints, 2) numpy array with the X,Y 
#                     coordinates of the DIC points to transform to 
#                     the FEM coordinate system
# Output: transDICPnts - A (nPoints,2) numpy array with the X, Y
#                        coordinates of the transformed dicPnts in the 
#                        FEM coordinate system
# -------------------------------------------------------------------
def transformDICPnts( transMat, dicPnts ):
    
    # Setup the memory for the X matrix used in the transformation
    # and setup the matrix
    nPnts = dicPnts.shape[0]
    X = np.ones( (nPnts, 3) )
    X[:,1:] = dicPnts
    
    # Transform the x-coordinates, then the y-coordinates
    newX = np.dot(X, transMat[0,:])
    newY = np.dot(X, transMat[1,:])
    
    # Return a matrix with the transformed coordinates
    return np.array( [newX, newY] ).transpose()
            
# Do some basic testing
print( getGenesisOP2Filename('platewithhole_FEM.dat', 0, True) )
changeMAT1Card('platewithhole_FEM.dat', 1, 6.898e10, 2.5e10)
runGENESIS( 'platewithhole_FEM.dat', unZipOP2=True )
print( getObjectiveFn( 'platewithhole_FEM.dat', 'exp_data.dat', 1 ) )

#femPnts = np.array( [[0., 0.],[1., 0.],[0., 1.],[1., 1.]] )
#dicPnts = np.array( [[1., 1.],[3., 1.],[1., 4.],[3., 4.]] )
#transMat  = setupTransform(femPnts, dicPnts)
#transPnts = transformDICPnts( transMat, dicPnts )
#print( dicPnts )
#print( transPnts )
