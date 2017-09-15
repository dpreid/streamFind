README for the streamFind set of algorithms for identifying tidal streams in GAIA mock catalogues.
All files (Aquarius .pe.dat and .info.dat) are in this folder.
Each simulation has its own folder of idstream files. the filenameGen script in the corresponding folder
generates the file variable for use in rVrClusterFind4 and phaseSpace_GAIA

Both algorithms are Octave programs, so commands are run in Octave.
REQUIRED PACKAGES:
statistics
image
signal
fits

1) muVt algorithm:
Adapt the pmTanV_init.m script depending on which simulation you want to run data for. 
Run the pmTanV_init.m script which initialises filenames, packages and puts you in the correct folder.
phaseSpace_GAIA.m runs the plotting program and outputs density plots and spatial plots of streams in
	tangential velocity - proper motion space.
[phaseSpace_GAIAerrors.m runs the same plots but with GAIA errors
phaseSpace_GAIAerrors_nomet.m adds error cuts to the data if required]
Fits images of the density plots are saved in density_fits folder.
These files can be transferred to a linux machine to run the CLUMPFIND function in the CUPID program of the
STARLINK package using:

CUPID is initialised (in the command terminal) using:

export STARLINK_DIR=~/star-2016A/
source $STARLINK_DIR/etc/profile
cupid

The fits file must be converted to NDF format for use with CUPID:

convert
fits2ndf in out

Then the program is run, in CUPID, using:

findclumps in=density_plot out=clumps config = ^config.dat method=fellwalker

Where config.dat is a text file with:

fellwalker.minheight=10
fellwalker.minpix=5

etc as the configuration file parameters for the fellwalker algorithm.

Clumps can be plotted and the number of each plotted over them using the KAPPA program (run command 'kappa'):

display clump_file
listshow clump_file plot=Text

and follow the commands

The original density plot output files are stored in density_fits.
The cluster plots after transfer to linux and back are in cluster_plots.

2) rVr method

Run the rVr_init.m script for initialising the catalogue address, packages and file names. NOTE: must change the data depending on 
whether ERRORLESS or ERROR CONVOLVED data is being used as well as the magnitude limit (by default set to Grvs < 16).

Run gmagConv.m in order to convert the G and V magnitude values from the mock catalogues into Grvs data.

Run rVrClusterFind4.m which will load radial velocity and radius data depending on the conditions stated in rVr_init.m; plot a full
plot of vrad against r; plot a filtered density plot; plot the peak finding plot; output the average number of peaks found (not accurate particularly) and save the extracted data to a file within the rVr folder.