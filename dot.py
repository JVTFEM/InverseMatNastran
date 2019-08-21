# Wrapper to call DOT from Python
#
#import os
import numpy as nm
import ctypes as ct
from ctypes import byref as B

class dot:

	#Set some local constants
	nInfo 	= 0
	nMethod = 0
	nPrint	= 0
	nMinMax	= 0
	nMaxInt	= 20000000
	nmParam = nm.empty(1, float)
	nmRPRM  = nm.zeros(20, float)
	nmIPRM  = nm.zeros(20, int)

	def __init__(self):
		#Load the shared DOT library
		self.dotlib = ct.cdll.LoadLibrary("libDOT2.so")

	# The DOT wrapper
	def dotcall(self, x, xl, xu, nCons):

		# Reset nInit
		nInit = 0

		#Initailize all array types
		nDvar = x.shape[0]
		ctDVAR	= ct.c_double * nDvar
		ctCONS	= ct.c_double * nCons
		ctRPRM	= ct.c_double * 20
		ctIPRM	= ct.c_int64 * 20

		#Initialize all arrays
		RPRM = ctRPRM(*(self.nmRPRM))		#Tells dot to use defaults
		IPRM = ctIPRM(*(self.nmIPRM))		#Tells dot to use defaults
		X    = ctDVAR(*(x))				#Initial values
		XL   = ctDVAR(*(xl))			#Lower bounds
		XU   = ctDVAR(*(xu))			#Upper bounds
		G    = ctCONS(*([0.0]*nCons))	#Constraints

		#Initialize constants
		METHOD  = ct.c_int64( self.nMethod )
		NDV     = ct.c_int64( nDvar )
		NCON    = ct.c_int64( nCons )
		IPRINT  = ct.c_int64( self.nPrint )
		MINMAX  = ct.c_int64( self.nMinMax )
		INFO    = ct.c_int64( self.nInfo )
		OBJ     = ct.c_double( 0.0 )
		MAXINT  = ct.c_int64( self.nMaxInt )

		# Call DOT510
		NRWK    = ct.c_int64()
		NRWKMN  = ct.c_int64()
		NRIWD   = ct.c_int64()
		NRWKMX  = ct.c_int64()
		NRIWK   = ct.c_int64()
		NSTORE  = ct.c_int64()
		NGMAX   = ct.c_int64()
		IERR    = ct.c_int64()

		self.dotlib.dot510_(B(NDV), B(NCON), B(METHOD), B(NRWK), B(NRWKMN), B(NRIWD), B(NRWKMX), B(NRIWK), B(NSTORE), B(NGMAX), B(XL), B(XU), B(MAXINT), B(IERR))

		ctRWK	= ct.c_double * NRWKMX.value
		ctIWK	= ct.c_int64 * NRIWK.value
		IWK	= ctIWK( *([0]*NRIWK.value) )
		WK	= ctRWK( *([0.0]*NRWKMX.value) )

		# Call DOT
		while (True):
			self.dotlib.dot_(B(INFO),B(METHOD),B(IPRINT), B(NDV),  B(NCON), B(X), B(XL), B(XU), B(OBJ), B(MINMAX), B(G), B(RPRM), B(IPRM), B(WK), B(NRWKMX), B(IWK), B(NRIWK))
			if ( INFO.value == 0 ) :
				break
			else:
				self.evaluate(X, OBJ, G, self.nmParam)

		rslt = nm.empty( 2+nDvar, float)
		rslt[0] = OBJ.value
		rslt[1] = 0.0
		if len(G) > 0 :
			rslt[1] = max(G)
		for i in range( nDvar ):
			rslt[2+i] = X[i]
		return rslt

	def evaluate(self, x, obj, g, param):
		obj.value = 2.0*(x[0]*x[1] + x[0]*x[2] + 2.0*x[1]*x[2])
		g[0] = 1.0 - 0.5*x[0]*x[1]*x[2]
		return
