import os
import numpy as np
import math
from scipy import interpolate
from decimal import *
import string
import pickle
import time
import pylab as P
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import matplotlib

from py21cmmc._21cmfast import initial_conditions, perturb_field, UserParams, CosmoParams, ionize_box, FlagOptions, AstroParams, spin_temperature, brightness_temperature, run_coeval, run_lightcone, compute_tau, compute_luminosity_function

np.seterr(invalid='ignore', divide='ignore')

# Setup ctypes
import ctypes
from ctypes import cdll
from ctypes import c_float,c_int,c_uint,c_double,c_char_p

# Location of the .so for the power spectrum module
SharedLibraryLocation = "../Compute_21cmPS.so"

Lib_21cmPS = cdll.LoadLibrary(SharedLibraryLocation)
Function_21cmPS = Lib_21cmPS.Compute_21cmPS

# Data structure for the returned data from the C code
class Point(ctypes.Structure):
    _fields_ = [('PSbins', c_int),
    	('PS_k', ctypes.POINTER(c_float)),
        ('PS', ctypes.POINTER(c_float)),
        ('PS_error', ctypes.POINTER(c_float))        
        ];

Function_21cmPS.restype = Point
Function_21cmPS.argtypes = [c_int, c_int, c_int, c_float, ctypes.POINTER(c_float)]

if __name__ == '__main__':


	# Setup and run 21cmFAST from Py21CMMC
	FFTW_WISDOM = False

	UParams = UserParams(USE_FFTW_WISDOM=FFTW_WISDOM,HMF=1)
	CParams = CosmoParams()
	AParams = AstroParams()
	FOptions = FlagOptions(USE_MASS_DEPENDENT_ZETA=True,USE_TS_FLUCT=True,INHOMO_RECO=True,SUBCELL_RSD=True,OUTPUT_AVE=True)

	RANDOM_SEED = 1

	# Final redshift for the light-cone box
	z_final = 6.0
	redshift = z_final*1.0001

	# Setup and return the initial conditions
	init_boxes = initial_conditions(user_params=UParams,cosmo_params=CParams, random_seed=RANDOM_SEED,regenerate=True)

	# Construct and return the perturbed field for the final redshift
	perturb_field_finalz = perturb_field(redshift=redshift,init_boxes=init_boxes)
	
	# Construct and return the 21cm light-cone
	LC = run_lightcone(redshift=redshift, user_params=UParams, cosmo_params=CParams, astro_params=AParams,  flag_options=FOptions, do_spin_temp=True, regenerate=True, init_box = init_boxes, perturb=perturb_field_finalz, use_interp_perturb_field=True)

	# Convert the 21cm light-cone into a cubic box
	brightness_temp = LC.brightness_temp
	
	# Chunch the lightcone into cubic boxes
	chunk_indices = list(range(0, LC.n_slices, UParams.HII_DIM))

	start = chunk_indices[0]
	end = chunk_indices[1]
	chunklen = (end-start) * LC.cell_size
	
	# Define the size of the data cube to be passed to C for constructing the 21cm power spectrum
	size = UParams.HII_DIM * UParams.HII_DIM * UParams.HII_DIM
	
	# Construct the array to be passed to C
	array_of_size_floats = c_float*size
	box_for_PS = array_of_size_floats()

	# Construct the 21cm power spectra from the resultant boxes
	LightConeFlag = 1
	DIMENSIONAL_T_POWER_SPEC = 1

	# Loop over the chunked cubes
	for ii in range(len(chunk_indices)-1):

		start = chunk_indices[ii]

		for i in range(UParams.HII_DIM):
			for j in range(UParams.HII_DIM):
				for k in range(UParams.HII_DIM):
					box_for_PS[k + (UParams.HII_DIM)*( j + UParams.HII_DIM*i )] = brightness_temp[i,j,k+start]

		output = Function_21cmPS(LightConeFlag,DIMENSIONAL_T_POWER_SPEC,UParams.HII_DIM,UParams.BOX_LEN,box_for_PS)

		# Output the 21cm power spectrum data
		for i in range(output.PSbins):
			print((output.PS_k[i]),(output.PS[i]),(output.PS_error[i]))	
