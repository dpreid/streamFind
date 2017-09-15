## Author: dpryd <dpryd@HAL9000JNR>
## Created: 2017-05-24
## Copyright (C) 2016 dpryd
## 
## This script calculates and plots multiple phase spaces
## This version is for the mock catalogues from the Aquarius and EAGLE simulations.
## 
## 
## This version plots with NO ERRORS in data.
## Use phaseSpace_GAIAerrors.m for ERROR CONVOLVED data plotting.
##
## USES PROPER MOTION VELOCITIES
## ============ USE =============
## Must run pmTanVinit.m prior to this script in order to initialise
## the file variable with filenames, as well as starting the appropriate packages.
## INPUT: Need to change the output filename, number of streams being plotted
## OUTPUT: DENSITY PLOTS which can be run in the CLUMPFIND program.

clear out_vrad_mean out_vrad_std out_metal out_infallz tanV metal stars_mul stars_mub L_x L_y L_z L_perp L2 L rad colour stars_vx stars_vy stars_vz stars_x stars_y stars_z starsA_x starsA_y starsA_z;

%======INITIALISATION=======
filename = "p000halo000_obs.dat"; 
%N_prog = 43;		% SET IN pmTanV_init.m number of files in the simulation, corresponding to the number of progenitor systems
nbins=600;			%resolution of the density plots
Ms = 1.989E30;	% solar mass in kg
mass = 0.7*Ms;	%mass of each RRLS particle

%r_min = 20;		%min and max distance for plotting
%r_max = 100;
colour = [];
out_metal = [];


%=======FILE LOAD AND DATA EXTRACTION=============

data = load (filename);		%load the required halo file

pm_l = data(:,1);		%spatial positions
pm_b = data(:,2);
parallax = data(:,3);
met = data(:,4);	%metallicities, Z
met_av = mean(met);	%mean metallicity of the progenitor 
particle_id = data(:,5);
mu = sqrt(pm_l.*pm_l + pm_b.*pm_b);	%overall proper motion
tan_vel = 4.75*mu./parallax;	%tangential velocity of each star
n=rows(pm_l);
% Select particle within peak in density
stream_ids = [];
for i=1:n
	if(tan_vel(i) > 150 && tan_vel(i) < 300)
		if(tan_vel(i) > (200/1.4)*pm_l(i) && tan_vel(i) < (210/1.1)*pm_l(i))
			stream_ids = [stream_ids;particle_id(i)];
		endif
	endif
endfor
disp("Number of stars in peak")
disp(rows(stream_ids))



%=========PLOTTING================


% figure 10
% scatter(stars_mul, stars_mub)
% xlabel("mu_lcosb/mas/yr","FontSize",20)
% ylabel("mu_b/mas/yr","FontSize",20)

% figure 11
% scatter3(rad,tanV,metal,"b")
% xlabel("heliocentric distance/kpc","FontSize",20)
% ylabel("tangential velocity/km/s","FontSize",20)
% zlabel("metallicity","FontSize",20)

% figure 12
% scatter(stars_mul,tanV)
% xlabel("mu_lcosb/mas/yr","FontSize",20)
% ylabel("tangential velocity/km/s","FontSize",20)


%============DENSITY PLOTS=================

%a = rot90(stars_mul,-1);		%gets the transpose vector (column to row)
%b = rot90(tanV,-1);
a=pm_l;
b=tan_vel;
M = [b,a];
density = hist3(M,[nbins,nbins]);	%count data stored in density
[X,Y] = meshgrid(linspace(min(a),max(a)),linspace(min(b),max(b)));	%meshgrid for plotting the colourmap


%====UNSHARP MASKING========
n_filt = 30;	% the array size for finding the median
smooth = medfilt2(density,true(n_filt));	% true(20) is a 20x20 matrix
mask = density - smooth;	% subtract the smoothed background from the original density plot
smask = mask./sqrt(smooth);	%sharpen the plot

%===SMOOTHING OF FINAL PLOT======
% mask image appears to be better than the smask (at least with low resolution data)
mask_filt = 2;
mask_median = medfilt2(mask,true(mask_filt));
mask_med2 = medfilt2(mask_median,true(mask_filt));

smooth2 = medfilt2(mask_median,true(n_filt));
mask_sub = mask_median - smooth2;


figure 18	%tangential velocity-proper motion plot
hold on
im = imagesc(X,Y,density)
set(gca(),"ydir","normal")
xlabel("mu_lcosb/mas/yr","FontSize",20)
ylabel("tangential velocity/km/s","FontSize",20)
title("full density plot","FontSize",20)


plot(pm_l,(200/1.4)*pm_l,"r")
plot(pm_l,(210/1.1)*pm_l,"r")
plot([-10,10],[150,150],"r")
plot([-10,10],[300,300],"r")
hold off

figure 19
im_median = imagesc(X,Y,smooth)
set(gca(),"ydir","normal")
xlabel("mu_lcosb/mas/yr","FontSize",20)
ylabel("tangential velocity/km/s","FontSize",20)
title("smoothed density plot","FontSize",20)

% figure 20
% im_mask = imagesc(X,Y,mask)
% set(gca(),"ydir","normal")
% xlabel("mu_lcosb/mas/yr","FontSize",20)
% ylabel("tangential velocity/km/s","FontSize",20)
% title("density - smooth","FontSize",20)

figure 21
hold on
im_mask_median = imagesc(X,Y,mask_median)
set(gca(),"ydir","normal")
xlabel("mu_lcosb/mas/yr","FontSize",20)
ylabel("tangential velocity/km/s","FontSize",20)
title("median filtered mask density plot","FontSize",20)
%plot(pm_l,(-300/1.05)*pm_l,"r")
%plot(pm_l,(-300/1.2)*pm_l,"r")
hold off
% figure 22					%proper motion density plot
% im_pm = imagesc(X_pm,Y_pm,density_pm,[0 5])
% set(gca(),"ydir","normal")
% xlabel("mu_lcosb/mas/yr","FontSize",20)
% ylabel("mu_b/mas/yr","FontSize",20)
% title("full density plot","FontSize",20)

% figure 23						% metallicity plot
% scatter(stars_mul,tanV,[],colour)
% xlabel("mu_lcosb/mas/yr","FontSize",20)
% ylabel("tangential velocity/km/s","FontSize",20)
% title("metallicity", "FontSize",20)

%======SAVE AND PRINT PLOTS AND FILES==========
% Outputs the density plots to screen and saves them to file.
%Converts to fits format for use with FINDCLUMPS
% Generates metallicity and infallz data in met_z and saves to file.
%print(18,strcat("C:/Users/dpryd/Desktop/streamFind/EAGLE_plots/",filename,"_tanV.png"),"-dpng")
%print(21,strcat("C:/Users/dpryd/Desktop/streamFind/EAGLE_plots/",filename,"_tanV_filtered.png"),"-dpng")
%print(23,strcat("C:/Users/dpryd/Desktop/streamFind/Aquarius_plots_final/",filename,"_tanV_metallicity.png"),"-dpng")

%save(strcat("C:/Users/dpryd/Desktop/streamFind/muVt/density_fits/",filename,"_tanV_density"),"density")
%save(strcat("C:/Users/dpryd/Desktop/streamFind/muVt/density_fits/",filename,"_tanV_filtered"),"mask_median")

%save_fits_image(strcat("C:/Users/dpryd/Desktop/streamFind/muVt/density_fits/",filename,"_tanV_density.fits"), density)
%save_fits_image(strcat("C:/Users/dpryd/Desktop/streamFind/muVt/density_fits/",filename,"_tanV_filtered.fits"), mask_median)

% mtl = rot90(out_metal,-1);
% iz = rot90(out_infallz,-1);
% o_vrad = rot90(out_vrad_mean,-1);
% o_vrad_std = rot90(out_vrad_std,-1);
% out_data =[mtl,iz,o_vrad,o_vrad_std];
% save(strcat("C:/Users/dpryd/Desktop/streamFind/muVt/stream_properties/",filename,"_out_data.dat"),"out_data")