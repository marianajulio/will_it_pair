import os
import sys
import json
import warnings
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import gaussian_kde
from scipy.interpolate import UnivariateSpline
from astropy.coordinates import SkyCoord, Galactocentric
import astropy.coordinates as coord
import astropy.units as u
import matplotlib.transforms as transforms
warnings.filterwarnings("ignore")

plt.rc('text', usetex=True)
plt.rcParams.update({'font.size': 18, 'font.family': 'monospace'})

# Constants
V_SUN = np.array([9, 232, 7]) # velocity of the sun relative to the galactic centre (km/s)
G = 4.30091e-6   # kpc (km/s)^2 / Msun
M31_GALACTOCENTRIC = np.array([-378.9, 612.7, -283.1])  # M31's position in Galactocentric coordinates (McConnachie+2012) in kpc
VEL_M31 = -301 #km/s (heliocentric, McConnachie+2012)

FRAC_DIST_ERROR = 0.03 # Fractional distance error for satellites without a distance uncertainty in the catalogue

N_MC_SAMPLES = 1000 # Number of Monte Carlo samples

markers = ['o', 's', 'D', '^', 'v', '<', '>', 'p', '*', '+', 'h', 'H', 'X', 'd', 'P', '8', '1']

# Observation repository (updated to match "Local Volume Database" v1.0.6)
OBS_REPO = 'https://raw.githubusercontent.com/apace7/local_volume_database/main/data/'

# Key names for the columns in the dataset
ID_NAME, OBJECT_NAME = 'key', 'name'
LON_GAL, LAT_GAL, DIST_GAL = 'll', 'bb', 'distance_gc' #degrees, degrees, kpc
RA, DEC = 'ra', 'dec'
DIST_SUN, E_DIST_P, E_DIST_M = 'distance', 'distance_ep', 'distance_em' #kpc, kpc, kpc
DIST_HOST, E_DIST_HOST_P, E_DIST_HOST_M = 'dist_host', 'dist_host_err_p', 'dist_host_err_m' #kpc, kpc, kpc (estimated in parse_data.py, since the host in the catalogue is sometimes considered another dwarf galaxy)
VLOS, E_VLOS_P, E_VLOS_M , VLOS_HOST = 'vlos_systemic', 'vlos_systemic_ep', 'vlos_systemic_em', 'vlos_rel_host' #km/s
PM_ALPHA, E_PM_ALPHA_P, E_PM_ALPHA_M = 'pmra', 'pmra_ep', 'pmra_em' #mas/yr
PM_DEC, E_PM_DEC_P, E_PM_DEC_M = 'pmdec', 'pmdec_ep', 'pmdec_em' #mas/yr
STELLAR_MASS = 'mass_stellar' #log; solar masses 
X, Y, Z = 'X', 'Y', 'Z'
VX, VY, VZ = 'VX', 'VY', 'VZ'
MAG_V = 'M_V'

# Output folder
RESULTS_DIR = 'results'
if not os.path.exists(RESULTS_DIR):
	   os.makedirs(RESULTS_DIR)

# Utility functions
def skyPlotAxes(ax):
	for deg in np.linspace(-60, 60, 5):
		ax.plot(np.linspace(-np.pi, np.pi, 360), np.repeat(np.deg2rad(deg), 360), 'k:', lw=0.75)
	for deg in np.linspace(-135, 135, 7):
		ax.plot(np.repeat(np.deg2rad(deg), 180), np.linspace(-np.pi / 2.0, np.pi / 2.0, 180), 'k:', lw=0.75)
		
	ax.set_xticks(np.linspace(-np.pi, np.pi, 9)[1:])
	ax.set_xticklabels(np.linspace(-180, 180, 9).astype(int)[1:])
	ax.set_yticks(np.deg2rad(np.linspace(-60, 60, 5)))
	ax.set_yticklabels(np.linspace(-60, 60, 5).astype(int))

# Checks that the folder given in the command line was created by initialise.py. If not, prints the usage and the available folders, and exits
def check_results_folder(folder, script_name):
	required = ['simulated_populations.opts', 'simulated_populations.data']
	if os.path.isdir(folder) and all(os.path.exists(os.path.join(folder, f)) for f in required):
		return
	message = f"Usage: python {script_name} {RESULTS_DIR}/<folder_name>\n'{folder}' is not a folder created by initialise.py."
	available = sorted(d for d in os.listdir(RESULTS_DIR) if os.path.exists(os.path.join(RESULTS_DIR, d, required[0])))
	if available:
		message += "\nAvailable folders:\n  " + "\n  ".join(f"{RESULTS_DIR}/{d}" for d in available)
	else:
		message += f"\nNo folders were found in {RESULTS_DIR}/. Please run initialise.py first."
	sys.exit(message)