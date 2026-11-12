from constants import *
import parse_data

markers = markers * 10

############# DYNAMICAL ANALYSIS ##############
# Estimates the minimum total pair mass (in 10^9 M_sun) required for the pair to be gravitationally bound based on line-of-sight velocity and 3D distance
def to_be_bound_NN(sat_1, sat_2, x, y, z, vlos):
	delta_r = np.sqrt((sat_1[x] - sat_2[x])**2 + (sat_1[y] - sat_2[y])**2 + (sat_1[z] - sat_2[z])**2)
	delta_v = abs(sat_1[vlos] - sat_2[vlos])
	M_pair = delta_r*delta_v**2/(2*G)
	
	return M_pair/1e9

# To get errors in our pair metrics
def get_stats(samples):
	p16, p50, p84 = np.nanpercentile(samples, [16, 50, 84])
	return p50, p50 - p16, p84 - p50   # median, lower error, upper error

# Draws N Monte Carlo realisations of one satellite
def mc_metrics(satellite, N):
	ra, dec = satellite[parse_data.RA], satellite[parse_data.DEC]
	pm_ra, e_pm_ra_p, e_pm_ra_m  = satellite[parse_data.PM_ALPHA], satellite[parse_data.E_PM_ALPHA_P], satellite[parse_data.E_PM_ALPHA_M]
	pm_dec, e_pm_dec_p, e_pm_dec_m = satellite[parse_data.PM_DEC], satellite[parse_data.E_PM_DEC_P], satellite[parse_data.E_PM_DEC_M]
	dist, e_dist_p, e_dist_m = satellite[parse_data.DIST_SUN], satellite[parse_data.E_DIST_P], satellite[parse_data.E_DIST_M]
	vlos, e_vlos_p, e_vlos_m = satellite[parse_data.VLOS], satellite[parse_data.E_VLOS_P], satellite[parse_data.E_VLOS_M]
	dist_host, e_dist_host_p, e_dist_host_m = satellite[parse_data.DIST_HOST], satellite[parse_data.E_DIST_HOST_P], satellite[parse_data.E_DIST_HOST_M]

	# Galactic coordinates for the MW host correction
	if host.lower() == 'mw':
		l_rad = np.deg2rad(satellite[parse_data.LON_GAL])
		b_rad = np.deg2rad(satellite[parse_data.LAT_GAL])
		rhat = np.array([np.cos(b_rad) * np.cos(l_rad), np.cos(b_rad) * np.sin(l_rad), np.sin(b_rad)])
		# Solar reflex correction
		correction = np.dot(rhat, V_SUN)

	# Distance, proper motions and line-of-sight velocity are sampled from Gaussians with sigma = (err_p + err_m)/2
	pm_ra_mc = pm_ra + 0.5 * (e_pm_ra_p + e_pm_ra_m) * np.random.normal(size=N)
	pm_dec_mc = pm_dec + 0.5 * (e_pm_dec_p + e_pm_dec_m) * np.random.normal(size=N)
	d_sun_mc = np.clip(dist + 0.5 * (e_dist_p + e_dist_m) * np.random.normal(size=N), 1e-3, None)
	vlos_mc = vlos + 0.5 * (e_vlos_p + e_vlos_m) * np.random.normal(size=N)
	d_host_mc = dist_host + 0.5 * (e_dist_host_p + e_dist_host_m) * np.random.normal(size=N)

	# Conversion to Galactocentric positions and velocities
	if host.lower() == 'mw':
		vlos_rel_host = vlos_mc + correction
	else:
		vlos_rel_host = vlos_mc - VEL_M31
	coords = coord.SkyCoord(ra=np.full(N, ra) * u.deg, dec=np.full(N, dec) * u.deg, distance=d_sun_mc * u.kpc, pm_ra_cosdec=pm_ra_mc * u.mas / u.yr, pm_dec=pm_dec_mc * u.mas / u.yr, radial_velocity=vlos_mc * u.km / u.s)
	galcen = coords.transform_to(Galactocentric())

	return {'X': galcen.x.to_value(u.kpc), 'Y': galcen.y.to_value(u.kpc), 'Z': galcen.z.to_value(u.kpc), 'VX': galcen.v_x.to_value(u.km / u.s), 'VY': galcen.v_y.to_value(u.km / u.s), 'VZ': galcen.v_z.to_value(u.km / u.s), 'DIST_HOST': d_host_mc, 'VLOS_HOST': vlos_rel_host, 'VLOS': vlos_mc, 'DIST_SUN': d_sun_mc}

# Combines the Monte Carlo samples of the two satellites of a pair (sample by sample) into the pair metrics
def estimate_mc_metrics(sat1_mc, sat2_mc):
	# Distance metrics
	# 3D separation
	dx, dy, dz = sat1_mc['X'] - sat2_mc['X'], sat1_mc['Y'] - sat2_mc['Y'], sat1_mc['Z'] - sat2_mc['Z']
	dr_3d_samples = np.sqrt(dx**2 + dy**2 + dz**2)
	# LOS distance difference
	dr_los_samples = np.abs(sat1_mc['DIST_SUN'] - sat2_mc['DIST_SUN'])

	# Velocity metrics
	# LOS velocity difference
	dv_los_samples = np.abs(sat1_mc['VLOS'] - sat2_mc['VLOS'])
	# 3D velocity difference
	dvx, dvy, dvz = sat1_mc['VX'] - sat2_mc['VX'], sat1_mc['VY'] - sat2_mc['VY'], sat1_mc['VZ'] - sat2_mc['VZ']
	dv_3d_samples = np.sqrt(dvx**2 + dvy**2 + dvz**2)

	# Host metics
	avg_d_host_samples = (sat1_mc['DIST_HOST'] + sat2_mc['DIST_HOST']) / 2.0
	avg_v_host_samples = ( sat1_mc['VLOS_HOST'] + sat2_mc['VLOS_HOST']) / 2.0
	
	return {'dr_3d': get_stats(dr_3d_samples), 'dr_los': get_stats(dr_los_samples), 'dv_3d': get_stats(dv_3d_samples), 'dv_los': get_stats(dv_los_samples), 'avg_d_host': get_stats(avg_d_host_samples), 'avg_v_host': get_stats(avg_v_host_samples)}


############# PLOTS ##############

# Pair diagnostics, coloured by the combined stellar mass of the pair
def plot_pair_diagnostics(n_pairs, dist, mstar, vlos, object_name):
	
	pair_masses = [np.log10(10**pairs_data.loc[pairs['sat_1'].iloc[i], mstar] + 10**pairs_data.loc[pairs['sat_2'].iloc[i], mstar]) for i in range(n_pairs)]
	norm = plt.Normalize(vmin=np.nanmin(pair_masses), vmax=np.nanmax(pair_masses))

	fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
	cmap = plt.cm.YlGnBu

	for i in range(n_pairs):
		sat_1, sat_2 = pairs_data.loc[pairs['sat_1'][i]], pairs_data.loc[pairs['sat_2'][i]]
		dr = np.linalg.norm(np.array([sat_1[parse_data.X], sat_1[parse_data.Y], sat_1[parse_data.Z]]) - np.array([sat_2[parse_data.X], sat_2[parse_data.Y], sat_2[parse_data.Z]]))
		dv_rel = abs(sat_1[vlos] - sat_2[vlos])
		stellar_mass_log = np.log10(10**sat_1[mstar] + 10**sat_2[mstar])

		pair_names = f"{sat_1[object_name]}-{sat_2[object_name]}"

		# 3D separation vs line-of-sight velocity difference
		ax1.scatter(dr, dv_rel, c=[stellar_mass_log], cmap=cmap, norm=norm, marker='o', edgecolor='k', s=130)
		ax1.annotate(pair_names, (dr+0.5, dv_rel), fontsize=8, alpha=0.7)
		# average distance to the host vs average line-of-sight velocity
		ax2.scatter((sat_1[dist]+sat_2[dist])/2, (sat_1[vlos]+sat_2[vlos])/2, c=[stellar_mass_log], cmap=cmap, norm=norm, marker='o', edgecolor='k', s=130)
		ax2.annotate(pair_names, ((sat_1[dist]+sat_2[dist])/2, (sat_1[vlos]+sat_2[vlos])/2), fontsize=8, alpha=0.7)

	ax1.set_xlabel(r'$\Delta R$ (kpc)', fontsize=16)
	ax1.set_ylabel(r'$\Delta v_{\mathrm{los}}$ (km/s)', fontsize=16)
	
	ax2.set_xlabel(r'$\overline{d}_{\mathrm{host}}$ (kpc)', fontsize=16)
	ax2.set_ylabel(r'$\overline{v}_{\mathrm{los}}$ (km/s)', fontsize=16)
	ax2.axhline(0, color='k', lw=1, ls=':')

	fig.subplots_adjust(right=0.88)
	cbar_ax = fig.add_axes([0.91, 0.15, 0.02, 0.7])
	sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
	cbar = fig.colorbar(sm, cax=cbar_ax)
	cbar.set_label(r'$\log_{10} M_{\star\mathrm{, pair}} \ (M_{\odot})$', fontsize=14)
	
	plt.savefig(os.path.join(simulated_population_folder, 'pair_diagnostics.pdf'), bbox_inches='tight')
	plt.show()

	return

# Distributions of pairs 
def plot_pairs_combined_distribution(n_pairs, ra, dec, vlos, mstar, object_name):
	fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 7), sharey=True)
	plt.subplots_adjust(wspace=0.1)
		
	v_limit = max(abs(satellites_data[vlos].min()), abs(satellites_data[vlos].max()))
	norm_v = plt.Normalize(vmin=-v_limit, vmax=v_limit)
	cmap_v = plt.get_cmap('RdBu_r')
		
	norm_m = plt.Normalize(vmin=satellites_data[mstar].min(), vmax=satellites_data[mstar].max())
	cmap_m = plt.get_cmap('YlGnBu')

	# Positions in RA/Dec of the satellites in pairs, in two panels coloured by line-of-sight velocity (left) and stellar mass (right)
	for i in range(n_pairs):
		# Each pair has its own marker; white markers mean the value is missing
		sat_1, sat_2 = pairs_data.loc[pairs['sat_1'][i]], pairs_data.loc[pairs['sat_2'][i]]
				
		for ax in [ax1, ax2]:
			#ax.plot([sat_1[ra], sat_2[ra]], [sat_1[dec], sat_2[dec]], '--', lw=1, color='grey', alpha=0.5, zorder=1)
			ax.set_xlabel(r'$\alpha$ (º)', fontsize=18)
				
			# vlos plot
			for sat in [sat_1, sat_2]:
				if np.isnan(sat[vlos]):
					ax1.scatter(sat[ra], sat[dec], facecolor='white', marker=markers[i], edgecolor='k', s=100)
					#ax1.annotate(sat[object_name], (sat[ra], sat[dec]), fontsize=10)
				else:
					ax1.scatter(sat[ra], sat[dec], c=[sat[vlos]], cmap=cmap_v, norm=norm_v, marker=markers[i], edgecolor='k', s=100)
					#ax1.annotate(sat[object_name], (sat[ra], sat[dec]), fontsize=10)

			# M_star plot
			for sat in [sat_1, sat_2]:
				if np.isnan(sat[mstar]):
					ax2.scatter(sat[ra], sat[dec], facecolor='white', marker=markers[i], edgecolor='k', s=100)
				else:
					ax2.scatter(sat[ra], sat[dec], c=[sat[mstar]], cmap=cmap_m, norm=norm_m, marker=markers[i], edgecolor='k', s=100)

	ax1.set_ylabel(r'$\delta$ (º)', fontsize=18)

	sm1 = plt.cm.ScalarMappable(cmap=cmap_v, norm=norm_v)
	cbar1 = fig.colorbar(sm1, ax=ax1, orientation='horizontal', pad=0.12, shrink=0.8)
	cbar1.set_label(r'$v_{\rm los}$ (km/s)', fontsize=14)

	sm2 = plt.cm.ScalarMappable(cmap=cmap_m, norm=norm_m)
	cbar2 = fig.colorbar(sm2, ax=ax2, orientation='horizontal', pad=0.12, shrink=0.8)
	label_m = r'$\log_{10}M_\star$ (M$_\odot$)'
	cbar2.set_label(label_m, fontsize=14)

	plt.savefig(simulated_population_folder+'/pairs_combined_vlos_mstar.pdf', bbox_inches='tight')
	plt.show()

# (Milky Way only) sky map in Galactic coordinates (l, b) showing the stellar mass of the satellite in pairs
def pairs_sky_distribution_mstar(n_pairs, l, b, mstar, object_name):
		
	norm = plt.Normalize(vmin=satellites_data[mstar].min(), vmax=satellites_data[mstar].max())
	cmap = plt.get_cmap('YlGnBu')

	plt.figure(figsize=(13, 10))

	axs = plt.subplot(projection='mollweide')
	skyPlotAxes(axs)
	axs.set_xticks(np.linspace(-np.pi, np.pi, 9)[1:] -.05)
	axs.set_xticklabels(np.linspace(0, 360, 9).astype(int)[1:])
	axs.set_xlabel(r'$l$ (º)',fontsize=18)
	axs.set_ylabel(r'$b$ (º)',fontsize=18)

	# satellites not in pairs are represented in grey
	plt.plot(np.radians(not_pairs_data[l]) - np.pi, np.radians(not_pairs_data[b]), c='grey', linestyle='None', marker = 'o', ms=7, alpha = 0.3)

	# pairs are joined by dashed lines and coloured by stellar mass
	for i in range(n_pairs):
		sat_1, sat_2 = pairs_data.loc[pairs['sat_1'][i]], pairs_data.loc[pairs['sat_2'][i]]
		axs.plot([np.radians(sat_1[l])-np.pi, np.radians(sat_2[l])-np.pi], [np.radians(sat_1[b]), np.radians(sat_2[b])], '--', lw=1, color='black', alpha=0.4, zorder=2)

		for sat in [sat_1, sat_2]:
			if np.isnan(sat[mstar]):
				axs.scatter(np.radians(sat[l])-np.pi, np.radians(sat[b]), facecolor='white', marker=markers[i], edgecolor='k', s=100)
			else:
				axs.scatter(np.radians(sat[l])-np.pi, np.radians(sat[b]), c=[sat[mstar]], cmap=cmap, norm=norm, marker=markers[i], edgecolor='k', s=100, zorder=10)
				
			axs.annotate(sat[object_name], (np.radians(sat[l])-np.pi + 0.05, np.radians(sat[b] - 0.05)), fontsize=11)


	sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
	label_text = r"$\log_{10}(M_\star/M_\odot)$"
	cbar = plt.colorbar(sm, ax=axs, orientation='horizontal', pad=0.08, shrink=0.7)
	cbar.set_label(label_text, fontsize=14)
	cbar.ax.xaxis.set_tick_params(labelsize=12)
				
	plt.savefig(simulated_population_folder+'/pairs_sky_distribution_mstar.pdf', bbox_inches='tight')
	plt.show()

	return

# Same as above, but for the line-of-sight velocity of the satellites in pairs
def pairs_sky_distribution_vlos(n_pairs, l, b, vlos, object_name):		

	v_limit = max(abs(satellites_data[vlos].min()), abs(satellites_data[vlos].max()))
	norm = plt.Normalize(vmin=-v_limit, vmax=v_limit)
	cmap = plt.get_cmap('RdBu_r')

	plt.figure(figsize=(13, 10))

	axs = plt.subplot(projection='mollweide')
	skyPlotAxes(axs)
	axs.set_xticks(np.linspace(-np.pi, np.pi, 9)[1:] -.05)
	axs.set_xticklabels(np.linspace(0, 360, 9).astype(int)[1:])
	axs.set_xlabel(r'$l$ (º)',fontsize=18)
	axs.set_ylabel(r'$b$ (º)',fontsize=18)
		
	plt.plot(np.radians(not_pairs_data[l]) - np.pi, np.radians(not_pairs_data[b]), c='grey', linestyle='None', marker = 'o', ms=7, alpha = 0.3)

	for i in range(n_pairs):
		sat_1, sat_2 = pairs_data.loc[pairs['sat_1'][i]], pairs_data.loc[pairs['sat_2'][i]]
		axs.plot([np.radians(sat_1[l])-np.pi, np.radians(sat_2[l])-np.pi], [np.radians(sat_1[b]), np.radians(sat_2[b])], '--', lw=1, color='black', alpha=0.4, zorder=2)

		for sat in [sat_1, sat_2]:
			if np.isnan(sat[vlos]):
				axs.scatter(np.radians(sat[l])-np.pi, np.radians(sat[b]), facecolor='white', marker=markers[i], edgecolor='k', s=100)
			else:
				axs.scatter(np.radians(sat[l])-np.pi, np.radians(sat[b]), c=[sat[vlos]], cmap=cmap, norm=norm, marker=markers[i], edgecolor='k', s=100, zorder=10)
				
			axs.annotate(sat[object_name], (np.radians(sat[l])-np.pi + 0.05, np.radians(sat[b] - 0.05)), fontsize=11)

	sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
	label_text = r"$v_{\mathrm{los}} (km/s)$"
	cbar = plt.colorbar(sm, ax=axs, orientation='horizontal', pad=0.08, shrink=0.7)
	cbar.set_label(label_text, fontsize=14)
	cbar.ax.xaxis.set_tick_params(labelsize=12)
				
	plt.savefig(simulated_population_folder+'/pairs_sky_distribution_vlos.pdf', bbox_inches='tight')
	plt.show()

	return

# 3D view of the satellites around the host
def plot_3d_satellite_system(satellites_data, pairs, x, y, z, vlos):
	fig = plt.figure(figsize=(14, 11))
	ax = fig.add_subplot(111, projection='3d')

	# to have a faint grid + background
	ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
	ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
	ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
	ax.xaxis._axinfo["grid"]['color'] =  (0.1, 0.1, 0.1, 0.1) 
	ax.yaxis._axinfo["grid"]['color'] =  (0.1, 0.1, 0.1, 0.1)
	ax.zaxis._axinfo["grid"]['color'] =  (0.1, 0.1, 0.1, 0.1)

	ax.plot(0, 0, 0, color='black', marker='x', markersize=10)

	# plot non-paired satellites in grey
	paired_ids = set(pairs['sat_1']).union(set(pairs['sat_2']))
	non_paired_data = satellites_data.loc[~satellites_data.index.isin(paired_ids)]
	ax.scatter(non_paired_data[x], non_paired_data[y], non_paired_data[z], color='lightgrey', alpha=0.5, s=50, label='Non-paired Satellites', zorder=1)

	# colormap for individual velocities of paired satellites
	v_limit = max(abs(satellites_data[vlos].min()), abs(satellites_data[vlos].max()))
	norm = plt.Normalize(vmin=-v_limit, vmax=v_limit)
	cmap = plt.get_cmap('RdBu_r')

	# plot the pairs
	for i in range(len(pairs)):
		sat_1, sat_2 = satellites_data.loc[pairs['sat_1'].iloc[i]], satellites_data.loc[pairs['sat_2'].iloc[i]]
		ax.plot([sat_1[x], sat_2[x]], [sat_1[y], sat_2[y]], [sat_1[z], sat_2[z]], color='black', alpha=0.5, ls='--', lw=1, zorder=5)

		for sat in [sat_1, sat_2]:
			if np.isnan(sat[vlos]):
				ax.scatter(sat[x], sat[y], sat[z], c='white', s=80, marker='o', edgecolor='black', linewidth=1, zorder=10)
			else:
				ax.scatter(sat[x], sat[y], sat[z], c=[sat[vlos]], cmap=cmap, norm=norm, s=80, marker='o', edgecolor='black', linewidth=1, zorder=10)

	sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
	cbar = plt.colorbar(sm, ax=ax, shrink=0.5, pad=0.1, aspect=15)
	cbar.set_label(r'$v_{\rm los}$ (km/s)', fontsize=14)

	ax.set_xlabel('X (kpc)')
	ax.set_ylabel('Y (kpc)')
	ax.set_zlabel('Z (kpc)')
	ax.set_box_aspect([1,1,1])
	plt.savefig(simulated_population_folder+'/3d_satellite_distribution.pdf', bbox_inches = 'tight', facecolor = 'w', edgecolor = 'w', transparent = False, format='pdf')	
	plt.show()
	
	return 

# To run the code
if __name__ == '__main__':
	if len(sys.argv) < 2:
		sys.exit(f"Usage: python check_pairs.py {RESULTS_DIR}/<folder_name>")
	simulated_population_folder = sys.argv[1]
	check_results_folder(simulated_population_folder, 'check_pairs.py')
	
	# Get the pairs previously discovered in nearest_neighbours.py
	pairs_path = os.path.join(simulated_population_folder, 'pairs.txt')
	# Check if the file exists and is not empty
	if not os.path.exists(pairs_path):
		sys.exit(f"The file {pairs_path} was not found. Please run nearest_neighbours.py on this folder first.")
	if os.path.getsize(pairs_path) == 0:
		sys.exit("nearest_neighbours.py was run, but no pairs were found, so there is nothing to analyse.")
	
	pairs = pd.read_csv(pairs_path, sep=r'\s+', names=['sat_1', 'sat_2'])
	# Sort each pair in alphabetic order
	pairs[['sat_1', 'sat_2']] = np.sort(pairs[['sat_1', 'sat_2']], axis=1)
	# Drop duplicate rows
	pairs = pairs.drop_duplicates().reset_index(drop=True)
	# Get the names of all the satellites that are in pairs
	satellites_in_pairs = pd.concat([pairs['sat_1'], pairs['sat_2']]).drop_duplicates().reset_index(drop=True).sort_values()
	
	# Get the choices made in initialise.py
	options_filename = os.path.join(simulated_population_folder, 'simulated_populations.opts')
	f = open(options_filename, 'r')
	options = json.load(f)
	host = options['host']

	# Parse the whole dataset that we get from the choices with the ID as the new index
	satellites_data = parse_data.parse(options, obs_message=False).set_index(parse_data.ID_NAME, drop=False)
	
	# Keep the data only for the satellites that are in pairs
	pairs_data = satellites_data.loc[satellites_data.index.isin(satellites_in_pairs)]
	# Keep the data only for the satellites that are *not* in pairs
	not_pairs_data = satellites_data.loc[~satellites_data.index.isin(satellites_in_pairs)]
	
	host = options['host']
	n_simulations = options['populations']
	n_satellites_in_pairs = len(pairs_data)
	n_satellites_not_in_pairs = len(not_pairs_data)
	n_pairs = len(pairs)

	print("Number of pairs: ", n_pairs)

	# Generate plots showing the spatial distribution and properties of the identified satellite pairs
	plot_3d_satellite_system(satellites_data, pairs, parse_data.X, parse_data.Y, parse_data.Z, parse_data.VLOS)

	# Compare the distributions of the identified pairs in different observables
	plot_pairs_combined_distribution(n_pairs, parse_data.RA, parse_data.DEC, parse_data.VLOS, parse_data.STELLAR_MASS, parse_data.OBJECT_NAME)
	# For the Milky Way, additionally show the sky distribution of pairs, colour-coded by stellar mass and line-of-sight velocity
	if host.lower() == 'mw':
		pairs_sky_distribution_mstar(n_pairs, parse_data.LON_GAL, parse_data.LAT_GAL, parse_data.STELLAR_MASS, parse_data.OBJECT_NAME)
		pairs_sky_distribution_vlos(n_pairs, parse_data.LON_GAL, parse_data.LAT_GAL, parse_data.VLOS, parse_data.OBJECT_NAME)
	
	# Analyse the properties of individual satellite pairs
    # Open output files to store the minimum mass required for gravitational binding and the calculated pair properties, respectively
	f_to_be_bound_pairs = open(simulated_population_folder+'/to_be_bound_pairs.txt','w')
	
	# header of the pair-properties file: each row contains the names of the two satellites, 
	# their relative separations, velocities, their average host distances and velocities, angular separation,
    # stellar mass ratio, minimum binding mass and dynamical mass-to-light ratio
	metrics_list = open(simulated_population_folder+'/metrics.txt','w')
	metrics_list.write(
			"sat1_name sat2_name "
			"delta_R_3d delta_R_3d_em delta_R_3d_ep "
			"delta_R_los delta_R_los_em delta_R_los_ep "
			"delta_v_3d delta_v_3d_em delta_v_3d_ep "
			"delta_v_los delta_v_los_em delta_v_los_ep "
			"avg_d_host avg_d_host_em avg_d_host_ep "
			"avg_v_host avg_v_host_em avg_v_host_ep "
			"ang_sep_deg mass_ratio "
			"m_bound ml_dyn\n")

	# Loop over all identified pairs and get the corresponding satellite data
	for i in range(n_pairs):
		satellite_1, satellite_2 = pairs_data.loc[pairs['sat_1'][i]], pairs_data.loc[pairs['sat_2'][i]]

		# Angular separation between the two satellites
		ra1, dec1 = np.radians(satellite_1[parse_data.RA]), np.radians(satellite_1[parse_data.DEC])
		ra2, dec2 = np.radians(satellite_2[parse_data.RA]), np.radians(satellite_2[parse_data.DEC])
		cos_sep = (np.sin(dec1) * np.sin(dec2) + np.cos(dec1) * np.cos(dec2) * np.cos(ra1 - ra2))
		ang_sep = np.degrees(np.arccos(np.clip(cos_sep, -1.0, 1.0)))

		# Stellar mass ratio defined as the mass of the more massive satellite divided by that of the less massive satellite
		m1, m2 = 10**satellite_1[parse_data.STELLAR_MASS], 10**satellite_2[parse_data.STELLAR_MASS]
		if np.isfinite(m1) and np.isfinite(m2) and min(m1, m2) > 0:
			m_ratio = max(m1, m2) / min(m1, m2)
		else:
			m_ratio = np.nan # If either stellar mass is missing or non-positive, assign NaN

		# Propagate observational uncertainties through Monte Carlo sampling to estimate the metrics and their uncertainties
		sat1_mc, sat2_mc = mc_metrics(satellite_1, N_MC_SAMPLES), mc_metrics(satellite_2, N_MC_SAMPLES)
		mc = estimate_mc_metrics(sat1_mc, sat2_mc)

		# Estimate the minimum total mass required for the pair to be gravitationally bound
		m_bound_val_mc = (mc['dr_3d'][0]*mc['dv_los'][0]**2/(2*G)) / 1e9
			
		# Calculate the dynamical mass-to-light ratio required for the pair to be bound 
		# (we assume a stellar mass-to-light ratio of M*/L_V = 2, consistent with the assumption used to estimate the stellar masses in the Local Volume Database)
		L_total = (m1 + m2) / 2.0
		ml_dyn = (m_bound_val_mc * 1e9) / L_total

		metrics_list.write(f"{satellite_1[parse_data.ID_NAME]} " f"{satellite_2[parse_data.ID_NAME]} "
				f"{mc['dr_3d'][0]:.1f} " f"{mc['dr_3d'][1]:.1f} " f"{mc['dr_3d'][2]:.1f} "
				f"{mc['dr_los'][0]:.1f} " f"{mc['dr_los'][1]:.1f} " f"{mc['dr_los'][2]:.1f} "
				f"{mc['dv_3d'][0]:.1f} " f"{mc['dv_3d'][1]:.1f} " f"{mc['dv_3d'][2]:.1f} "
				f"{mc['dv_los'][0]:.1f} " f"{mc['dv_los'][1]:.1f} " f"{mc['dv_los'][2]:.1f} "
				f"{mc['avg_d_host'][0]:.1f} " f"{mc['avg_d_host'][1]:.1f} " f"{mc['avg_d_host'][2]:.1f} "
				f"{mc['avg_v_host'][0]:.1f} " f"{mc['avg_v_host'][1]:.1f} " f"{mc['avg_v_host'][2]:.1f} "
				f"{ang_sep:.1f} " f"{m_ratio:.1f} " f"{m_bound_val_mc:.1f} " f"{ml_dyn:.0f}\n")
		
		# Calculate the minimum mass required for the pair to be gravitationally bound using the original satellite measurements
		m_bound_val = to_be_bound_NN(satellite_1, satellite_2, parse_data.X, parse_data.Y, parse_data.Z, parse_data.VLOS)
		
		# Prints the binding mass if all required measurements are available
		if np.isfinite(m_bound_val):
			print(f"For {satellite_1[parse_data.OBJECT_NAME]} to be gravitationally bound to its nearest neighbour {satellite_2[parse_data.OBJECT_NAME]} the mass of the pair has to be {m_bound_val:.3g} 10^9 solar masses.\n")
		# Identify satellites with missing measurements and provide a reason for the inability to estimate the minimum mass
		else:
			missing = [s[parse_data.OBJECT_NAME] for s in (satellite_1, satellite_2) if np.isnan(s[parse_data.VLOS])]
			reason = f"no line-of-sight velocity is available for {' and '.join(missing)}" if missing else "some of the required measurements are missing"
			print(f"⚠︎ The minimum mass for {satellite_1[parse_data.OBJECT_NAME]} to be gravitationally bound to its nearest neighbour {satellite_2[parse_data.OBJECT_NAME]} cannot be estimated because {reason}.\n")

		# Save the pair IDs and the estimated minimum binding mass (NaN is written if the required measurements are unavailable)
		f_to_be_bound_pairs.write("%s %s %f\n" % (satellite_1[parse_data.ID_NAME], satellite_2[parse_data.ID_NAME], m_bound_val))
	
	f_to_be_bound_pairs.close()
	metrics_list.close()

	plot_pair_diagnostics(n_pairs, parse_data.DIST_HOST, parse_data.STELLAR_MASS, parse_data.VLOS_HOST, parse_data.OBJECT_NAME)