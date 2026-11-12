from constants import *

# Transforms sky coordinates (RA, Dec, Distance) into MW- or M31-centric cartesian coordinates (X, Y, Z, distance to host)
def radec_to_xyz(host, n_satellites, ra, dec, dist_sun): 
	x, y, z, d = [], [], [], []
	for i in range(n_satellites):
		coords = SkyCoord(ra = ra[i]*u.deg, dec = dec[i]*u.deg, distance = dist_sun[i]*u.kpc,frame='icrs').transform_to('galactocentric')
		if host == 'mw':
			x.append(coords.x.value) 
			y.append(coords.y.value)
			z.append(coords.z.value)
			d.append(np.sqrt(coords.x.value**2 + coords.y.value**2 + coords.z.value**2))

		elif host == 'm31':
			coords_m31centric = np.array([coords.x.value, coords.y.value, coords.z.value]) - M31_GALACTOCENTRIC
			x.append(coords_m31centric[0])
			y.append(coords_m31centric[1])
			z.append(coords_m31centric[2])	
			d.append(np.linalg.norm(coords_m31centric))

	return x, y, z, d

# Transforms sky positions, proper motions, and line-of-sight velocities into Galactocentric space velocity components (VX, VY, VZ)
def pmrapmdec_to_vxvyvz(ra, dec, distance, pmra, pmdec, vlos):
	coords = SkyCoord(ra=ra.values*u.degree, dec=dec.values*u.degree, distance=distance.values*u.kpc, pm_ra_cosdec=pmra.values*u.mas/u.yr, pm_dec=pmdec.values*u.mas/u.yr, radial_velocity=vlos.values*u.km/u.s, frame='icrs')
	galcen = coords.transform_to(Galactocentric())
	
	return galcen.v_x.value, galcen.v_y.value, galcen.v_z.value

# Loads observation data from local_volume_database, computes 3D coordinates and space velocities, and applies the observational cuts defined in initialise.py
# Returns a pandas dataframe with the final dataset
def parse_data_observations(dataset, host, filters):
	df = pd.read_csv(dataset)

	# Satellites without a distance uncertainty in the catalogue get a fractional error (proportional to the heliocentric distance)
	assumed_error = FRAC_DIST_ERROR * df[DIST_SUN]
	df['dist_err_assumed'] = df[[E_DIST_M, E_DIST_P]].isna().any(axis=1)
	df[E_DIST_M] = df[E_DIST_M].fillna(assumed_error)
	df[E_DIST_P] = df[E_DIST_P].fillna(assumed_error)

	# Coordinate transformations
	N_SATELLITES = len(df)
	x, y, z, d = radec_to_xyz(host, N_SATELLITES, df[RA], df[DEC], df[DIST_SUN])
	_, _, _, dist_low_raw = radec_to_xyz(host, N_SATELLITES, df[RA], df[DEC], df[DIST_SUN]-df[E_DIST_M])
	_, _, _, dist_high_raw = radec_to_xyz(host, N_SATELLITES, df[RA], df[DEC], df[DIST_SUN]+df[E_DIST_P])
	dist_sorted = np.sort(np.vstack([d, dist_low_raw, dist_high_raw]).T, axis=1)

	# Velocity transformations
	vx, vy, vz = pmrapmdec_to_vxvyvz(df[RA], df[DEC], df[DIST_SUN], df[PM_ALPHA], df[PM_DEC], df[VLOS])
	df = df.assign(X=x, Y=y, Z=z, dist_host=d, dist_host_err_m=dist_sorted[:, 1]-dist_sorted[:, 0], dist_host_err_p=dist_sorted[:, 2]-dist_sorted[:, 1])
	df = df.assign(VX=vx, VY=vy, VZ=vz)

	# Keep only the confirmed dwarfs
	df = df.loc[df['confirmed_galaxy'] == 1].reset_index(drop=True)
	# Remove the dwarfs that don't have distance measurements
	df = df.dropna(subset=DIST_HOST).reset_index(drop=True)
	# Remove the dwarfs that are possibly not bound to their host (we set a distance threshold of 350 kpc, which is the approximate virial radius of the MW and M31)
	df = df.loc[df[DIST_HOST] <= 350.0].reset_index(drop=True)

	# Apply the optional cuts
	if 'stellar_mass_min' in filters:
		df = df.loc[df[STELLAR_MASS] >= float(filters['stellar_mass_min'])].reset_index(drop=True)
	if 'magnitude_threshold' in filters:
		df = df.loc[df[MAG_V] <= float(filters['magnitude_threshold'])].reset_index(drop=True)

	# Correction of the line-of-sight velocity of the satellites for the motion of the Sun (for the MW) or for the motion of M31 (for M31) to get the velocity of the satellites relative to their host
	if host.lower() == 'mw':
		l_rad, b_rad = np.radians(df[LON_GAL]), np.radians(df[LAT_GAL])
		rhat = np.column_stack([np.cos(b_rad) * np.cos(l_rad), np.cos(b_rad) * np.sin(l_rad), np.sin(b_rad)])
		correction = np.sum(rhat * V_SUN, axis=1)
		vlos_rel_host = df[VLOS]+ correction
	else:
		vlos_rel_host = df[VLOS] - VEL_M31
	df = df.assign(vlos_rel_host=vlos_rel_host)

	return df

def parse(options, obs_message = False):
	dataset = f'dwarf_{options["host"]}'
	filename = OBS_REPO + dataset + '.csv'

	return parse_data_observations(filename, options['host'], options['filters'])

# Calculate the cumulative number of satellites from the radial distribution of the host and fit the observed CDF
def fit_cdf(distance, method):
	observed_distances_sorted = np.sort(distance)
	observed_satellite_counts = np.arange(1, len(distance) + 1) / len(observed_distances_sorted)
	
	if method == 'kde':
		fit = gaussian_kde(observed_distances_sorted)
	elif method == 'spline':
		fit = UnivariateSpline(observed_distances_sorted, observed_satellite_counts, s=0.001*len(distance))
	
	return fit

# Plot the CDFs and their fits
def plot_cdf(distance):
	observed_distances_sorted = np.sort(distance)
	n_satellites = np.arange(1, len(distance) + 1)
	x_cdf = np.linspace(np.min(observed_distances_sorted), np.max(observed_distances_sorted), 1000)
	
	fit_cdf_spline = fit_cdf(distance, method = 'spline')
	y_spline = fit_cdf_spline(x_cdf)*len(distance)
	fit_cdf_kde = fit_cdf(distance, method = 'kde')
	y_kde = fit_cdf_kde(x_cdf)
	
	plt.figure()
	plt.plot(observed_distances_sorted, n_satellites, lw=2, color = '#a7d0e4', label ='Observed CDF')
	plt.plot(x_cdf, y_spline, '--', lw=2, color='#225da8', label = 'Fit')
	plt.xlabel('Radial distribution (kpc)')
	plt.ylabel('Number of satellites')
	plt.legend()
	plt.title("Fitting method: spline")
	plt.tight_layout()
	plt.show(block=False)
	
	plt.figure()
	plt.hist(observed_distances_sorted, density = True, linewidth=2, facecolor='#a7d0e4', edgecolor='#a7d0e4', fill=True, alpha= 0.4, label='PDF')
	plt.plot(x_cdf, y_kde, '--', lw=2, color='#225da8', label = 'Fit')
	plt.xlabel('Radial distribution (kpc)')
	plt.ylabel('Probability density function')
	plt.legend()
	plt.title("Fitting method: kde")
	plt.tight_layout()
	plt.show(block=False)
	
	return
	
# Building the name of the output files
def build_sim_name(options):
	magnitude = f"mt{options['filters']['magnitude_threshold']}_" if 'magnitude_threshold' in options['filters'] else ""
	stellar_mass = f"smgt{options['filters']['stellar_mass_min']}_" if 'stellar_mass_min' in options['filters'] else ""
	sim_name = f"{options['host'].lower()}_{stellar_mass}{magnitude}{options['populations']}"
	return sim_name


