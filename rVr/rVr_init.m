## Author: dpryd <dpryd@HAL9000JNR>## Created: 2017-06-08## This script performs the setup for running the phaseSpace_GAIA.m script for## plotting in proper motion and tangential velocity space.
##======= USAGE===========## Running this script loads the required packages for the plotting programs,## initialises the file names for the appropriate catalogue.## INPUT: need to change the catalogue_address variable in setup variables to ## point at the appropriate catalogue address where a filenameGen script should## also exist which will generate the correct file names for each IDstream.## After this script is run, you can use rVrClusterFind4.m and phaseSpace_GAIA.m.
more off	%allows for disp function to output during program
##======= SETUP VARIABLES========
errors = "N";		%YES(Y) or NO(N)
sim_name = "Aquarius A2 KIII, Grvs < 16, no errors";
datafile = "AQA2KIII_Grvs16_r20100_noerrors.dat";
%N_prog = 22;	%number of progenitor idstreams in the simulation
catalogue_address = "C:/Users/dpryd/Desktop/streamFind/mock_catalogues_final/AquariusA2_KIII_pe_byIDstream/";
script_address = "C:/Users/dpryd/Desktop/streamFind/rVr/"; 
%filename = "AQA2KIII_r20_100_noerrors";
%ind_stream_file = "C:/Users/dpryd/Desktop/streamFind/Aquarius_plots_final/individual_streams_plots/AQE2_RRLS/AquariusE2RRLS.";
nbins=600;			%resolution of the density plots
r_min = 20;		%min and max distance for plotting
r_max = 100;
mag_limit = 16;		% Grvs mag limit
radVel_limit = 600;	%limit in measured radial velocity
%radVel_err = 0.5;	% error cuts for radial velocity and parallax
mub_err = 0.5;		


pkg load statistics		%loads the required packages
pkg load signal
pkg load image
pkg load fits

cd (catalogue_address)	% changes to appropriate folders for initialising the filenames
filenameGen
cd (script_address)

