from constants import *

# Draws a sample of n_satellites radial distances (matching the size of the observed satellite sample) using the fitted distribution (spline or kde)
def radial_distribution(n_satellites, distance, method):
	observed_distances_sorted = np.sort(distance)
	observed_satellite_counts = np.arange(1, len(distance) + 1) / len(observed_distances_sorted)
	
	if method == 'spline':
		fit = UnivariateSpline(observed_distances_sorted, observed_satellite_counts, s=0.001*len(distance))
		r_values = np.linspace(np.min(distance), np.max(distance), 1000)
		# Force the fitted CDF to be monotonic and to go from 0 to 1, so it can be inverted safely
		cdf_values = np.maximum.accumulate(np.clip(fit(r_values), 0, 1))
		cdf_values = (cdf_values - cdf_values[0]) / (cdf_values[-1] - cdf_values[0])
		generated_r_values = np.interp(np.random.uniform(0, 1, n_satellites), cdf_values, r_values)

	elif method == 'kde':
		fit = gaussian_kde(observed_distances_sorted)
		generated_r_values = fit.resample(n_satellites).flatten()
	
	return generated_r_values

# Generates isotropic mock distributions of satellites around MW or M31. Includes Galactic plane obscuration mask (|b| > 15 deg) when the host is the MW 
def generate_isotropic_population(host, n_satellites, n_simulations, distance, method):
	generated_coords = []
	obscuration_latitude_cut = 15
	if host.lower() == 'mw': # To include the obscuration of the MW
		print(f"\nNote: since your host is the Milky Way, an obscuration of ±{obscuration_latitude_cut}º will be included when generating the position of the satellites.")
		for _ in range(n_simulations):
			# Coordinates that pass the obscuration filter
			valid = np.empty((0, 3))
			
			# To generate the number of desired satellites while applying the obscuration filter
			while len(valid) < n_satellites:
				# Generate more than needed, since some will be rejected
				m = 2 * n_satellites 
				
				# Generate random spherical coordinates
				theta = np.arccos(1 - 2 * np.random.uniform(size=m)) # Random inclination
				phi = 2 * np.pi * np.random.uniform(size=m) # Random azimuthal angle
				r = radial_distribution(m, distance, method) # Radial distribution (in kpc)

				# Convert spherical to Cartesian coordinates (in kpc)
				x = r * np.sin(theta) * np.cos(phi) 
				y = r * np.sin(theta) * np.sin(phi)
				z = r * np.cos(theta)

				# Decide which points survive the obscuration
				keep = np.abs(np.arcsin(z / r)) > np.radians(obscuration_latitude_cut)
				valid = np.vstack([valid, np.column_stack([x, y, z])[keep]])
			
			generated_coords.append(valid[:n_satellites])
	else:
		for _ in range(n_simulations):
			# Generate random spherical coordinates
			theta = np.arccos(1 - 2 * np.random.rand(n_satellites))  # Random inclination
			phi = 2 * np.pi * np.random.rand(n_satellites)           # Random azimuthal angle
			r = radial_distribution(n_satellites, distance, method)       # Radial distribution (in kpc)
	
			# Convert spherical to Cartesian coordinates (in kpc)
			x = r * np.sin(theta) * np.cos(phi) 
			y = r * np.sin(theta) * np.sin(phi)
			z = r * np.cos(theta)
	
			generated_coords.append(np.vstack((x, y, z)).T)
			
	return generated_coords
