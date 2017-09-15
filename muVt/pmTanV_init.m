## Copyright (C) 2017 dpryd
## 
## This program is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.


## Author: dpryd <dpryd@HAL9000JNR>
## Created: 2017-06-08

## This script performs the setup for running the phaseSpace_GAIA.m script for
## plotting in proper motion and tangential velocity space.

##======= USAGE===========
## Running this script loads the required packages for the plotting programs,
## initialises the file names for the appropriate catalogue.
## INPUT: need to change the catalogue_address variable in setup variables to 
## point at the appropriate catalogue address where a filenameGen script should
## also exist which will generate the correct file names for each IDstream.
## After this script is run, you can use rVrClusterFind4.m and phaseSpace_GAIA.m.

more off	%allows for disp function to output during program
##======= SETUP VARIABLES========
%N_prog = 43;	%AQA2 RRLS number of progenitor idstreams in the simulation
%N_prog = 30;	%AQB2
N_prog = 54;	%AQC2
%N_prog = 45;	%AQD2
%N_prog = 22;	%AQE2
%N_prog = 166;	%AQA2 KIII
%N_prog = 77;	%AQB2
%N_prog = 204;	%AQC2
%N_prog = 127;	%AQD2 KIII tracer stars
%N_prog = 83;	%AQE2

catalogue_address = "C:/Users/dpryd/Desktop/streamFind/mock_catalogues_final/AquariusC2_RRLS_pe_byIDstream/";
script_address = "C:/Users/dpryd/Desktop/streamFind/muVt/"; 
filename = "AQC2RRLS_noerrors";
ind_stream_file = "C:/Users/dpryd/Desktop/streamFind/Aquarius_plots_final/individual_streams_plots/AQC2_RRLS/AquariusC2RRLS";
nbins=600;			%resolution of the density plots
r_min = 20;		%min and max distance for plotting
r_max = 100;
mub_err = 0.7;	%relative error in proper motion
tan_vel_max = 700;

pkg load statistics		%loads the required packages
pkg load signal
pkg load image
pkg load fits

cd (catalogue_address)	% changes to appropriate folders for initialising the filenames
filenameGen
cd (script_address)

