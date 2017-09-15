## Copyright (C) 2016 dpryd
##
## Latest version: includes the conversions between G magnitude and Grvs, since the RVS
## limit is the limit we want Grvs < 17. I assume that V-Ic is positive for all KIII tracer stars
## NEED TO FIND BETTER EVIDENCE OF THIS.

## MUST RUN gmagConv.m prior to this program

##Found that the relative parallax error was a constant 0.2 in data.
## This version calculates a relative error based on the simulated data vs the GAIA
## error data and removes stars if too high.
## Also includes cut in radial velocity to below 500km/s


## filenameGen.m needs to be run first in order to get the appropriate
## names for the simulation files into the variable file.
## Three steps to the program:
## 1. Plotting radial velocity - radius plots 
## This program is a full implementation of my algorithm for identifying
## tidal streams in radial velocity - radius space. 
## INPUT: .dat files with position and radial velocity data as well as metallicity
## OUTPUT: Radial velocity - radius plots, density plots and number of tidal streams identified
##
## 2. Plotting density plots
## Plotting conditions added in order to only plot the clusters
## Unsharp mask filter applied
## 

## must load the statistics package for this to work correctly
## pkg load statistics
## pkg load image
## pkg load signal
##
## 3. Counting clusters 
## This program takes star density count data from densityPlot2
## and attempts to count the number of clusters in radial velocity space
## 

%======INITIALISATION======
clear rad r_v clusters;

% sim_type = "A";		% type H (HYDRO) or A (Aquarius)
% tracer = "RRLS";	%tracer type, either KIII or RRLS
% errors = "Y";		%YES(Y) or NO(N)
% r_min = 20;			%min and max radius for plotting
% r_max = 100;
% mag_limit = 16;		% Grvs mag limit
% error_limit = 0.3;	% largest relative error in parallax to be plotted

%nbins = 300;		%number of bins for density plot

%file = load("hydro009files.dat");		%list of filenames for the data files


for i=1:size(file)(1)
	clear rad r_v clusters;
	rad = [];
	data = load(file(i,:));		%load first file into data from variable file generated by filenameGen.m
	
	G = data(:,3);	%GAIA measured magnitude G < 20 observable by GAIA
	V = data(:,2);	% Johnson V magnitude
	GV = G-V;
	n = size(G);		%number of star particles
	GGrvs = polyval(ppos,GV);	% G-Grvs colour from polynomial fit run in gmagConv.m
	Grvs = G - GGrvs;			% Grvs magnitude for mag limit cut.
	%======NO ERRORS========
	radVel = data(:,18);
	radius = data(:,15);	%heliocentric distance
	parallax = data(:,12);	%simulated parallax
	%======GAIA ERRORS========
	gradVel = data(:,31);
	gradius = data(:,28); 	%error convolved heliocentric distance
	gparallax = data(:,25);	% GAIA measured parallax
	%=======RELATIVE ERRORS IN DATA========
	radVel_relerr = data(:,37);	%radial velocity relative error
	%radVel_limit = 500;	%limit in measured radial velocity
	par_relerr = data(:,5);		%parallax relative error
	relerr_mub = data(:,36);	%proper motion relative error
	%par_relerr = abs((gparallax - parallax)./parallax);	%absolute relative error
	%=======PLOTTING===========
	if (errors == "N")	%choose to plot error free data or not
		for i=1:n
			if ((radius(i) > r_min) && (radius(i) < r_max) && Grvs(i) <= mag_limit)	%plotting conditions
				rad(end+1) = radius(i);
				r_v(end+1) = radVel(i);
			endif
		endfor
	else
		for i=1:n
			if ((gradius(i) > r_min) && (gradius(i) < r_max) && (Grvs(i) <= mag_limit) && (relerr_mub(i) <= mub_err) && (abs(gradVel(i)) <= radVel_limit)) %(radVel_relerr(i) <= radVel_err)) %plotting conditions
				rad(end+1) = gradius(i);
				r_v(end+1) = gradVel(i);
			endif
		endfor
	endif
	if (rows(rad) > 0)
		figure (i)
		scatter(rad,r_v)
		xlabel("radius/kpc","FontSize",20)
		ylabel("radial velocity/km/s","FontSize",20)
		title(i,"FontSize",20)
	endif
endfor

% savedata = [rot90(rad,-1),rot90(r_v,-1)];	%data to be saved to file
% save("-ascii",datafile,"savedata")
