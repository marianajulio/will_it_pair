from constants import *
import parse_data

# Determine NN distances for observed satellites
def compute_nn_distances(observed_coords):		
	tree = cKDTree(observed_coords)
	nn_distances, nn_indices = tree.query(observed_coords, k=2)  # k=2 because k=1 is distance to self (0.0)
	return nn_distances[:, 1], nn_indices[:, 1]  # Return both distances and indices while removing the distance to self (indice [:, 0])
					  
# Function to use the errors on the distance of the satellite for which we want to obtain the distribution of NN
def include_error_distance(d, err_m, err_p, n_samples=10000):
	sigma = 0.5 * (err_m + err_p)
	return d + np.random.normal(0, sigma, n_samples)

# Estimate NN distances in relation to a specific satellite using the isotropic distribution around its host	
def compute_simulated_nn_distances(simulated_coords, satellite_index, observed_coords):
	target_pos = observed_coords[satellite_index]
	
	# Calculate distance to target for all generated satellites across all isotropic populations at once
	diff = simulated_coords - target_pos
	distances = np.linalg.norm(diff, axis=2) 
	
	# Take the minimum distance for each simulation
	min_distances = np.min(distances, axis=1)
	return min_distances

# Estimate the probability of finding a generated nearest neighbor closer than the one observed
def probability_nn(simulated_nn_distances, observed_nn_distance):
	return np.mean(simulated_nn_distances < observed_nn_distance)

# Estimate (and plot) the probability of finding a generated nearest neighbor closer than the one observed, considering the distance errors
def pvalues_NN(probability_nominal, probability_low, probability_high):

	probability_nominal_sorted = np.sort(probability_nominal)
	probability_nominal_counts = np.arange(1, len(probability_nominal_sorted) + 1) / len(probability_nominal_sorted)

	probability_low_sorted = np.sort(probability_low)
	probability_low_counts = np.arange(1, len(probability_low_sorted) + 1) / len(probability_low_sorted)

	probability_high_sorted = np.sort(probability_high)
	probability_high_counts = np.arange(1, len(probability_high_sorted) + 1) / len(probability_high_sorted)

	plt.figure()
	plt.plot(probability_nominal_sorted, probability_nominal_counts, linewidth=2, color="#c6e9b4") 
	plt.fill_betweenx(probability_nominal_counts, probability_low_sorted, probability_high_sorted, color='#c6e9b4', alpha=0.3, zorder=30) 
	plt.plot([0, 1], [0, 1], c="k", ls="--")
	plt.xlabel("p-value")
	plt.ylabel("CDF")
	plt.savefig(simulated_population_folder + "/CDF_pvalues_NN.pdf", bbox_inches="tight", facecolor="w", edgecolor="w", transparent=False, orientation="portrait", format="pdf")
	plt.show()

	plt.figure()
	plt.hist(probability_nominal, density=True, alpha=0.4, linewidth=2, color="#c6e9b4", edgecolor="#c6e9b4", label="Nominal")
	plt.hist(probability_low, density=True, histtype="step", linewidth=2, color="#a7d0e4", label="Low")
	plt.hist(probability_high, density=True, histtype="step", linewidth=2, color="#225da8", label="High")
	plt.axhline(1, color="k", linestyle="--")
	plt.xlabel("p-value")
	plt.ylabel("Density")
	plt.legend()
	plt.savefig(simulated_population_folder + "/hist_pvalues_NN.pdf", bbox_inches="tight", facecolor="w", edgecolor="w", transparent=False, orientation="portrait", format="pdf")
	plt.show()

	with open(simulated_population_folder + "/p_values.txt", "w") as f:
		f.write("# p_nominal cdf_nominal p_16 cdf_16 p_84 cdf_84\n")

		for i in range(len(probability_nominal_sorted)):
			f.write("%f %f %f %f %f %f\n" % (
					probability_nominal_sorted[i], probability_nominal_counts[i],
					probability_low_sorted[i], probability_low_counts[i],
					probability_high_sorted[i], probability_high_counts[i]))

	return


############# PLOTS ##############

# Plot the cumulative radial distribution of the observed satellites and compare it with the distributions obtained from the simulated populations
def plot_cdf_radial_distribution_filled(method, observed_distance_host, observed_distance_host_e_min, observed_distance_host_e_max, simulated_coords):
	
	# Calculate the lower and upper bounds of the observed radial distances by applying the corresponding observational uncertainties
	observed_lower_bound = observed_distance_host - observed_distance_host_e_min
	observed_upper_bound = observed_distance_host + observed_distance_host_e_max
	
	# Sort the observed distances and calculate the CDF (after applying uncertainties)
	sorted_observed_distances = np.sort(observed_distance_host)
	# Calculate the CDF rank
	n_sats = len(observed_distance_host)
	observed_rank = np.arange(1, n_sats + 1)

	# Extract and sort radial distances for all simulations
	sim_radii = np.linalg.norm(simulated_coords, axis=2)  # Shape: (n_simulations, n_satellites)
	all_sim_distances = np.sort(sim_radii, axis=1)

	# Calculate the median radial distance at each rank across all simulations and define the 1, 2 and 3 sigma ranges, respectively
	median_sim = np.median(all_sim_distances, axis=0)
	low_1s, up_1s = np.percentile(all_sim_distances, [16, 84], axis=0)
	low_2s, up_2s = np.percentile(all_sim_distances, [2.5, 97.5], axis=0)
	low_3s, up_3s = np.percentile(all_sim_distances, [0.15, 99.85], axis=0)

	plt.figure(figsize=(8, 6))

	plt.fill_betweenx(observed_rank, low_3s, up_3s, color='lightgrey', alpha=0.2)
	plt.fill_betweenx(observed_rank, low_2s, up_2s, color='lightgrey', alpha=0.4)
	plt.fill_betweenx(observed_rank, low_1s, up_1s, color='darkgrey', alpha=0.5)

	plt.plot(median_sim, observed_rank, color='dimgrey', lw=2, label='Median generated CDF')

	# Sort the lower and upper observational distance bounds independently to obtain the corresponding uncertainty bands for the observed distribution
	sorted_observed_lower_bound = np.sort(observed_lower_bound)
	sorted_observed_upper_bound = np.sort(observed_upper_bound)

	plt.fill_betweenx(observed_rank, sorted_observed_lower_bound, sorted_observed_upper_bound, color='#5e96b2', alpha=0.4, zorder=30)
	plt.plot(sorted_observed_distances, observed_rank, color="#5e96b2", label='Observed CDF', zorder=40)

	plt.ylabel('Number of satellites', fontsize=20)
	plt.xlabel('Radial distance (kpc)', fontsize=20)
	plt.legend(loc='lower right') 
	plt.tight_layout()
	
	plt.savefig(os.path.join(simulated_population_folder, f'simulated_CDFs_{method}_filled.pdf'), bbox_inches='tight', dpi=300, facecolor='w', edgecolor='w', transparent=False, orientation='portrait', format='pdf')   
	plt.show()

	return

# Plot satellite distribution around the host
def plot_satellite_distribution(n_satellites, sat_coordinates, object_name):
	fig = plt.figure(figsize=(10, 10))
	ax = fig.add_subplot(projection='3d')
	ax.set_aspect('equal')
	ax.set_xlabel("X (kpc)")
	ax.set_ylabel("Y (kpc)")
	ax.set_zlabel("Z (kpc)")

	ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0)), ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0)), ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
	ax.xaxis._axinfo["grid"]['color'], ax.yaxis._axinfo["grid"]['color'], ax.zaxis._axinfo["grid"]['color'] = (0.1, 0.1, 0.1, 0.1), (0.1, 0.1, 0.1, 0.1), (0.1, 0.1, 0.1, 0.1)
	
	ax.plot(0, 0, 0, 'x', c='k')

	for i in range(n_satellites):
		ax.scatter(sat_coordinates[i, 0], sat_coordinates[i, 1], sat_coordinates[i, 2], c='#40b5c4')
		ax.text(sat_coordinates[i, 0], sat_coordinates[i, 1], sat_coordinates[i, 2], f'{object_name[i]}', fontsize = 9)
	
	plt.savefig(os.path.join(simulated_population_folder, f'real_satellite_distribution.pdf'), bbox_inches = 'tight', facecolor = 'w', edgecolor = 'w', transparent = False, orientation='portrait', format='pdf')
	plt.show()
	
	return

# Plot the distribution of nearest-neighbour distances for a specific satellite, including the observed distance and the 16th and 84th percentiles	
def plot_histogram_nn(satellite_index, simulated_nn_distances, observed_nn_distance, sep_low, sep_high, observed_nn_index, object_name, probability_nominal, probability_low, probability_high, nbins):
	fig, ax = plt.subplots(figsize=(10, 7))
	trans = transforms.blended_transform_factory(ax.transData, ax.transAxes)

	ax.hist(simulated_nn_distances, bins=nbins, density=True, linewidth=2, facecolor="#40b5c4", edgecolor="#40b5c4", fill=True, alpha=0.3)

	# Observed NN separation
	ax.axvline(observed_nn_distance, color="#225da8", linestyle="-", linewidth=2)
	# 16th and 84th percentiles of observed pair separation
	ax.axvline(sep_low, color="#225da8", linestyle=":", linewidth=2)
	ax.axvline(sep_high, color="#225da8", linestyle=":", linewidth=2)

	ax.text(sep_low - 1, 0.98, f"{object_name[observed_nn_index]}", transform=trans, color="black", rotation=90, va="top", ha="right", fontsize=14)

	# Box with p-values 
	textstr = (r"$P_\mathrm{min}\,=\;%.2f$" "\n" r"$P\,=\;%.2f$" "\n" r"$P_\mathrm{max}=\;%.2f$") % (probability_low, probability_nominal, probability_high)
	props = dict(facecolor="white", edgecolor="lightgrey") 
	ax.text(0.97, 0.95, textstr, transform=ax.transAxes, fontsize=14, va="top",ha="right", bbox=props)

	ax.set_xlabel("Distance to the nearest-neighbour (kpc)", fontsize=16)
	ax.set_ylabel("Density", fontsize=16)
	ax.set_title(f"{object_name[satellite_index]}", fontsize=18)

	plt.savefig(os.path.join(simulated_population_folder, f"{DATASET[parse_data.ID_NAME][satellite_index]}_NN_distribution.pdf"), bbox_inches="tight", facecolor="w", edgecolor="w", transparent=False, format="pdf")
	plt.close()

	return

# To run the code
if __name__ == '__main__':
	if len(sys.argv) < 2:
		sys.exit(f"Usage: python nearest_neighbours.py {RESULTS_DIR}/<folder_name>")
	simulated_population_folder = sys.argv[1]
	check_results_folder(simulated_population_folder, 'nearest_neighbours.py')
	if not os.path.isdir(simulated_population_folder):
		sys.exit(f"The folder {simulated_population_folder} does not exist. Please check the folder name (it is printed at the end of initialise.py).")
	
	options_filename = os.path.join(simulated_population_folder, 'simulated_populations.opts')
	f = open(options_filename, 'r')
	options = json.load(f)

	data_filename = os.path.join(simulated_population_folder, 'simulated_populations.data')
	f = open(data_filename, 'rb')
	simulated_coords = np.load(f)
	
	DATASET = parse_data.parse(options, obs_message=False)

	host = options['host']
	n_simulations = options['populations']
	method = options['fitting_method']
	n_satellites = len(DATASET)
	coords = np.stack((DATASET[parse_data.X], DATASET[parse_data.Y], DATASET[parse_data.Z])).T
	coords_helio = np.stack((DATASET[parse_data.RA], DATASET[parse_data.DEC], DATASET[parse_data.DIST_SUN])).T
					
	nbins = 15 if n_simulations < 1000 else 30

	# Plot the real satellite distribution
	plot_satellite_distribution(n_satellites, coords, DATASET[parse_data.OBJECT_NAME])

	# Compute nearest-neighbour distances for observed satellites
	observed_nn_distances, observed_nn_indices = compute_nn_distances(coords)

	# Plot CDF of the simulated radial distribution against the observed one
	plot_cdf_radial_distribution_filled(method, DATASET[parse_data.DIST_HOST], DATASET[parse_data.E_DIST_HOST_M], DATASET[parse_data.E_DIST_HOST_P], simulated_coords)

	# Obtain the nearest-neighbour distributions for each satellite and evaluate pair significance
	prob, prob_low, prob_high = [], [], []
	pairs_file_path = os.path.join(simulated_population_folder, 'pairs.txt')
	nn_summary_path = os.path.join(simulated_population_folder, 'nn_pvalues_summary.txt')

	with open(pairs_file_path, 'w') as f_pairs, open(nn_summary_path, 'w') as f_summary:
		# header of summary file: (target satellite, nearest neighbour, nominal probability, lower and upper probabilities)
		f_summary.write("# ID_target ID_NN p_nominal p_low p_high\n")

		# loop over all observed satellites and identify their nearest neighbours
		for i in range(n_satellites):
			nn_idx = observed_nn_indices[i]
			print(f"{DATASET[parse_data.OBJECT_NAME][i]} has nearest neighbour {DATASET[parse_data.OBJECT_NAME][nn_idx]} at a distance of {observed_nn_distances[i]:.2f} kpc.\n")

			# Generate Monte Carlo samples of the heliocentric distances of both satellites accounting for their observational uncertainties
			n_samples = 10000
			d1_samples = include_error_distance(coords_helio[i, 2], DATASET[parse_data.E_DIST_M][i], DATASET[parse_data.E_DIST_P][i], n_samples)
			d2_samples = include_error_distance(coords_helio[nn_idx, 2], DATASET[parse_data.E_DIST_M][nn_idx], DATASET[parse_data.E_DIST_P][nn_idx], n_samples)

			# Construct sky coordinates for each Monte Carlo realisation using the sampled distances and the observed angular positions
			c1 = SkyCoord(ra=coords_helio[i, 0]*u.deg, dec=coords_helio[i, 1]*u.deg, distance=d1_samples*u.kpc, frame='icrs')
			c2 = SkyCoord(ra=coords_helio[nn_idx, 0]*u.deg, dec=coords_helio[nn_idx, 1]*u.deg, distance=d2_samples*u.kpc, frame='icrs')

			# Calculate the distribution of 3D separations between the two satellites and use its 16th and 84th percentiles to define the lower and upper bounds
			separations = c1.separation_3d(c2).value
			sep_low, sep_high = np.percentile(separations, [16, 84])

			# Calculate the nearest-neighbour distances in each simulated population
			simulated_nn_distances = compute_simulated_nn_distances(simulated_coords, i, coords)
			
			# Estimate how frequently the simulated nearest-neighbour distances are smaller than the observed separation
			probability_no_err = probability_nn(simulated_nn_distances, observed_nn_distances[i]) # uses the measured separation
			probability_low = probability_nn(simulated_nn_distances, sep_low) # uses the lower bound from the Monte Carlo distribution
			probability_high = probability_nn(simulated_nn_distances, sep_high) # uses the upper bound from the Monte Carlo distribution

			# Store results
			prob.append(probability_no_err)
			prob_low.append(probability_low)
			prob_high.append(probability_high)

			# Plot the simulated nearest-neighbour distance distribution, the observed separation and its uncertainty bounds
			plot_histogram_nn(i, simulated_nn_distances, observed_nn_distances[i], sep_low, sep_high, observed_nn_indices[i], DATASET[parse_data.OBJECT_NAME], probability_no_err, probability_low, probability_high, nbins)

			# Identify pairs that are significantly close (p <= 0.25)
			if probability_no_err <= 0.25:
				print(f"  -> {DATASET[parse_data.OBJECT_NAME][i]} is suspiciously close to {DATASET[parse_data.OBJECT_NAME][nn_idx]} (p = {probability_no_err:.4f})!\n")

				# Save the satellite IDs of the identified pair
				f_pairs.write(f"{DATASET[parse_data.ID_NAME][i]} {DATASET[parse_data.ID_NAME][nn_idx]}\n")
				
				# Save the nominal probability and the probabilities calculated using the lower and upper bounds of the observed separation
				f_summary.write("%s %s %.6f %.6f %.6f\n" % (DATASET[parse_data.ID_NAME][i], DATASET[parse_data.ID_NAME][nn_idx], probability_no_err, probability_low, probability_high))

	# Estimate and plot p-values
	pvalues_NN(prob, prob_low, prob_high)
	
	print("Nearest-neighbour analysis complete. You can now run 'check_pairs.py' to explore the properties of the identified pairs.")
	