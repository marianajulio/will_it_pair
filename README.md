# Will it pair? Satellite galaxies edition
**Will it pair?** finds the satellite galaxies of the Milky Way or Andromeda that are unusually close to their nearest neighbour. It compares each observed nearest-neighbour separation with the distribution expected from $N$ isotropic populations that have the same number of satellites and the same radial distribution, and flags the pairs that are unlikely to arise by chance as candidates. For each candidate pair, it then estimates the minimum mass the pair would need to be gravitationally bound and computes properties such as separations, velocity differences and stellar-mass ratio, together with diagnostic plots.

For a detailed description of the method, check the accompanying paper!

⚠️  If you use WIP?SGE, don't forget to cite the following papers:
+ Code implementation - [Júlio et al. submitted]
+ Original NN methodology - [Fattahi et al. 2013](https://ui.adsabs.harvard.edu/abs/2013MNRAS.431L..73F/abstract)
+ [Local Volume Database](https://github.com/apace7/local_volume_database) - [Pace 2025](https://ui.adsabs.harvard.edu/abs/2025OJAp....8E.142P/abstract)

## Getting started
To download WIP?SGE, please go to your preferred folder and clone the repository using:
```
git clone https://github.com/marianajulio/will_it_pair.git
```

Run all the scripts from inside the repository folder, since some files and the `results/` folder are accessed through relative paths.

### Requirements

Before you begin using WIP?SGE you will need to download and install some dependencies. The code has been verified to work with:
+ Local Volume Database 1.1.2 
+ python 3.9.20
+ numpy version 1.24.4
+ scipy version 1.13.1
+ matplotlib version 3.9.2
+ astropy version 6.0.1
+ pandas version 2.3.3

*(Creating a new environment is recommended, for example with `conda create -n wipsge python=3.9.20` followed by `conda activate wipsge`, and then installing the packages from inside the repository folder with `pip install -r requirements.txt`.)*

You will also need:
+ A working LaTeX installation, since the plots use LaTeX for the text (`usetex=True` in `constants.py`).
+ An internet connection, since the observational data are downloaded from the Local Volume Database every time the code is run.

**Alternatively, if you want to work without internet connection:** download the files `dwarf_mw.csv` and `dwarf_m31.csv` from the [data folder of the Local Volume Database](https://github.com/apace7/local_volume_database/tree/main/data), save them in a folder of your choice (for example, `data/observations/` inside the repository), and in `constants.py` change

```python
OBS_REPO = 'https://raw.githubusercontent.com/apace7/local_volume_database/main/data/'
```
to
```python
OBS_REPO = 'data/observations/'
```
The path is relative to the folder from which you run the scripts, so run them from the repository folder, or use an absolute path instead. The file names must be `dwarf_mw.csv` and `dwarf_m31.csv`.

Since the local files will not be updated automatically, it is a good idea to note the date (or version) of the database you downloaded, so that your results can be reproduced.

## Files

### Generating the isotropic populations

In order to run the code, first you need to generate the isotropic populations that mimic your observed distribution. To do that, you need to run:
```
python initialise.py
```
When you run this file, you have to answer some questions, such as:
+ which dataset do you wish to use (i.e., satellites of which host)
+ stellar mass threshold cut (if any)
+ magnitude threshold cut (if any)
+ the method used to model the radial distribution (spline or KDE)
+ the number of random isotropic populations to generate ($N$)

⚠️  For the menu questions, enter **only** the number of the option you want. For the stellar mass cut, the magnitude cut and the number of populations, enter a numerical value. If an answer is not valid, you will be asked again. You can type ```q``` at any prompt to quit.

❗  If your host is the Milky Way, an obscuration will be automatically applied, so no satellites will be generated within $|b| < 15^\circ$.

At the end of `initialise.py`, a folder will be created inside `results/` with a name that summarises your choices. This name will be printed when the code finishes running.

*(If you would like to implement additional selection criteria, you can modify `initialise.py` and update the `parse_data()` function in `parse_data.py` accordingly.)*

#### Selection of satellites
Besides the optional cuts above, the code always keeps only confirmed galaxies that have a distance measurement and lie within 350 kpc of their host. The optional cuts work as follows:
+ stellar mass cut: keeps galaxies with $\log_{10}(M_\star/M_\odot)$ greater than or equal to the value you enter;
+ magnitude cut: keeps galaxies with absolute magnitude $M_V$ less than or equal to the value you enter (i.e., removes the galaxies that are fainter).

❗ At least 5 satellites must remain after the cuts, otherwise `initialise.py` will ask you to choose less strict cuts.

**Distance uncertainties:** after the paper was finished, new satellites were added to the [Local Volume Database](https://github.com/apace7/local_volume_database) that do not yet have a measured distance uncertainty. The sample used in the paper does not contain any of these cases. To allow the code to handle them, we assume a fractional error of 3% on the heliocentric distance for any satellite without a distance uncertainty in the catalogue. The value is set in `constants.py` (`FRAC_DIST_ERROR`). These assumed errors only affect the uncertainties (the error band in the radial distance CDF, `p_low` and `p_high`, and the error bars in `metrics.txt`). The nominal p-values are not affected.

### Nearest neighbours

To get the distribution of expected nearest neighbours for your observed satellites, you can now run:
```
python nearest_neighbours.py results/<folder name>
```
*( ⚠️ Don't forget to change `<folder name>` to the name of the folder that was previously created - and printed - at the end of `initialise.py`.)*

The following files are stored in the same folder:

+ `pairs.txt`: the pairs of satellites that are unusually close;
+ `p_values.txt`: the sorted p-values of all satellites and the corresponding CDF (nominal, low and high estimates);
+ `nn_pvalues_summary.txt`: the p-values of each identified pair, with the IDs of the two satellites;
+ plots: the CDF and histogram of the p-values, the observed radial distribution compared with the simulated ones, the observed 3D distribution of the satellites, and the distribution of expected nearest-neighbour distances for each satellite.

### Check the properties of the pair candidates

You can now run:
```
python check_pairs.py results/<folder name>
```
*( ⚠️ Don't forget to change `<folder name>` to the name of the folder that was previously created by `initialise.py` and used to run `nearest_neighbours.py`.)*

This file will determine the minimum mass that each pair would need to have to be gravitationally bound, as described in the paper. It will also output some plots that are helpful to analyse the pairs (for the Milky Way, this includes sky maps in Galactic coordinates), and a file `metrics.txt` with several properties of each pair. They will also be stored in the folder that was created by `initialise.py`.

Quantities that cannot be computed because of missing measurements (e.g. line-of-sight velocities) are written as `nan`.

### Notes

In the paper, you have the option to generate isotropic populations for both the observed data of the Local Group and data from HESTIA simulations. Since the latter is not public, the code was slightly modified to account **only for the observations**.

The code is connected to the [Local Volume Database](https://github.com/apace7/local_volume_database), so when you run it, it always considers the most up-to-date version of the database, and results may change slightly if the database is updated. The results in the paper were obtained with the last access to the database on 12 March 2026.

If you want to use a different dataset, you need to change the variable `OBS_REPO` defined in the `constants.py` file and change the names of the quantities that come from your dataset accordingly. Note that `parse_data.py` also uses one column name directly (`confirmed_galaxy`) and builds the file name from the host (`dwarf_mw.csv` or `dwarf_m31.csv`), so these need to be updated as well.

#### Final remarks
If you find any bugs or have any questions, please [contact me](mailto:mpouseirojulio@aip.de).


Mariana P. Júlio

12 November 2026