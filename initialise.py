from constants import *
import parse_data
import generate_isotropic_population

# Minimum number of satellites needed after the cuts (cubic spline fit + nearest-neighbour analysis)
MIN_SATELLITES = 5
QUIT_WORDS = ('q', 'quit', 'exit')
YES_NO = {"1": True, "2": False}

# input() wrapper: strips whitespace and lets the user quit at any prompt by typing q.
def read_input(prompt):
	try:
		text = input(prompt).strip()
	except (EOFError, KeyboardInterrupt):
		print("\nExiting.")
		sys.exit(0)
	if text.lower() in QUIT_WORDS:
		print("Exiting.")
		sys.exit(0)
	return text

# Keeps asking until the answer is one of the keys of choices and returns the mapped value
def ask_choice(prompt, choices):
	while True:
		text = read_input(prompt).rstrip('.').strip()
		if text in choices:
			return choices[text]
		print(f"\n⚠︎ Invalid option! Please enter one of: {', '.join(choices)} (or q to quit).\n")

# Keeps asking until the answer is a finite number of the requested type within [min_value, max_value]
def ask_number(prompt, cast, kind, min_value, max_value, integer=False):
	while True:
		text = read_input(prompt)
		try:
			value = cast(text)
		except ValueError:
			print(f"\n⚠︎ Invalid option! '{text}' is not a valid {kind}. Please try again (or q to quit).\n")
			continue
		if not np.isfinite(value):
			print(f"\n⚠︎ Invalid option! The value must be finite. Please try again.")
			continue
		if not (min_value <= value <= max_value):
			print(f"\n⚠︎ Invalid option! The value must be between {min_value} and {max_value}. Please try again.\n")
			continue
		if integer and not value.is_integer():
			print(f"\n⚠︎ Invalid option! The value must be a whole number. Please try again.\n")
			continue
		return int(value) if integer else value

# Asks for the optional cuts and returns them as a dictionary
def ask_filters():
	filters = {}

	print("\nWould you like to apply a stellar mass cut?")
	print("1. Yes")
	print("2. No")
	if ask_choice("Option: ", YES_NO):
		print("\nRemove dwarfs with log10(M*/M☉) < ...  (enter a number, e.g. 5 or 6.5)")
		filters['stellar_mass_min'] = ask_number("Option: ", float, "number", 0, 12)

	print("\nWould you like to apply any magnitude threshold?")
	print("1. Yes")
	print("2. No")
	if ask_choice("Option: ", YES_NO):
		print("\nRemove dwarfs with absolute magnitude M_V fainter than ...  (enter a number, e.g. -8 or -10.5)")
		filters['magnitude_threshold'] = ask_number("Option: ", float, "number", -25, 5)

	return filters

if __name__ == '__main__':
	options = {}

	try:
		with open('banner.txt', 'r') as banner_file:
			print(banner_file.read())
	except FileNotFoundError:
		pass

	read_input("Press enter to start (you can type q at any prompt to quit)")

	# Dataset
	print("\nWhich dataset do you want to work with?")
	print("1. Milky Way")
	print("2. Andromeda")
	options['host'] = ask_choice("Option: ", {"1": "mw", "2": "m31"})

	# Cuts: repeat until the user is happy and enough satellites remain
	while True:
		options['filters'] = ask_filters()

		try:
			df = parse_data.parse(options, obs_message=True)
		except Exception as e:
			print(f"\nCould not load the dataset ({e}). Please check your internet connection and try again.")
			sys.exit(1)

		if len(df) < MIN_SATELLITES:
			print(f"\nOnly {len(df)} satellite galaxies remain after these cuts, but at least {MIN_SATELLITES} are needed. Please choose less strict cuts.")
			continue

		print(f"\nAfter the cuts that were applied to the dataset, {len(df)} satellite galaxies remain.")
		print("1. Continue")
		print("2. Change the cuts")
		print("3. Quit")
		decision = ask_choice("Option: ", {"1": "continue", "2": "redo", "3": "quit"})
		if decision == "quit":
			print("Exiting.")
			sys.exit(0)
		if decision == "continue":
			break

	parse_data.plot_cdf(df[parse_data.DIST_HOST])
	plt.pause(0.1)  # makes sure the figures are drawn before the next prompt

	# Fitting method
	print("\nWhich fitting method for the radial distribution of the satellites would you like to use?")
	print("1. spline")
	print("2. kde")
	options['fitting_method'] = ask_choice("Option: ", {"1": "spline", "2": "kde"})

	# Number of populations
	print("\nHow many isotropic populations would you like to generate? (whole number, e.g. 1e3 or 1000)")
	options['populations'] = ask_number("Number of populations: ", float, "whole number", 1, 1e7, integer=True)

	generated_coords = generate_isotropic_population.generate_isotropic_population(options['host'], len(df), options['populations'], df[parse_data.DIST_HOST], options['fitting_method'])

	# Path to save results
	sim_name = parse_data.build_sim_name(options)
	output_folder = os.path.join(RESULTS_DIR, sim_name)
	os.makedirs(output_folder, exist_ok=True)
	with open(os.path.join(output_folder, "simulated_populations.data"), 'wb') as f:
		np.save(f, generated_coords)
	with open(os.path.join(output_folder, "simulated_populations.opts"), 'w') as f:
		json.dump(options, f)

	print(f"\nDone! {options['populations']} isotropic populations were simulated based on your choices.")
	print(f"Your simulated populations and options were stored in the folder {output_folder}.")
	print("You can now run nearest_neighbours.py with these files.")