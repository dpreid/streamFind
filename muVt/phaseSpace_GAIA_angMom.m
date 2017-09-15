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
%filename = "AQA2RRLS_r10_100_noerrors"; SET IN pmTanV_init.m
%N_prog = 43;		% SET IN pmTanV_init.m number of files in the simulation, corresponding to the number of progenitor systems
% nbins=600;			%resolution of the density plots
Ms = 1.989E30;	% solar mass in kg
mass = 0.7*Ms;	%mass of each RRLS particle

% r_min = 10;		%min and max distance for plotting
% r_max = 100;
colour = [];
out_metal = [];
out_infallz = [];
out_vrad_mean = [];
out_vrad_std = [];
% tot_mass = 0.0;					
% for i=1:n			%find total mass
	% tot_mass = tot_mass + mass(i);
% endfor
% velCoMx = 0.0;
% velCoMy = 0.0;
% velCoMz = 0.0;
% for i=1:n			%calculate the centre of mass velocity
	% velCoMx = velCoMx + mass(i)*v_x(i);
	% velCoMy = velCoMy + mass(i)*v_y(i);
	% velCoMz = velCoMz + mass(i)*v_z(i);
% endfor

% velCoMx = velCoMx/tot_mass;		%need to divide by total mass to get CoM velocity
% velCoMy = velCoMy/tot_mass;
% velCoMz = velCoMz/tot_mass;

%=======FILE LOAD AND DATA EXTRACTION=============

for j=1:N_prog
	disp("Plotting stream:")
	disp(j)
	data = load (file(j,:));		%load the required halo file
	%info = load(file_info(j,:));	% load the corresponding info file with metallicity etc
	
	x= data(:,6);		%spatial positions
	y = data(:,7);
	z = data(:,8);
	%met = info(:,5);	%metallicities, Z
	%met_av = mean(met);	%mean metallicity of the progenitor 
	%out_metal(end+1)=met_av;
	%infallz = info(:,8);	%redshift when the progenitor is considered to have fallen into the host galaxy.
	%out_infallz(end+1)=infallz(1);
	mul = data(:,16);	%proper motion in galactic longitude, mulcosb
	mub = data(:,17);	% proper motion in galactic latitude
	v_x = data(:,9);	%velocities/km/s
	v_y = data(:,10);
	v_z = data(:,11);
	v_rad = data(:,18);		%radial velocity
	out_vrad_mean(end+1) = mean(v_rad);	%calculate mean and dispersion of v_rad
	out_vrad_std(end+1) = std(v_rad);
	rhel = data(:,15);	%heliocentric distance/kpc
	parallax = data(:,12); % parallax
	mu = sqrt(mul.*mul + mub.*mub);	%overall proper motion
	tan_vel = 4.75*mu./parallax;	%tangential velocity of each star
	n=size(x);			%number of particles in each progenitor
	for i=1:n
		if(rhel(i) > r_min && rhel(i) < r_max)	%only plot particles within these spatial ranges
			rad(end+1) = rhel(i); 
			L_x(end+1) = (-y(i).*v_z(i) + z(i).*v_y(i));	%angular momentum components per unit mass
			L_y(end+1) = (-z(i).*v_x(i) + x(i).*v_z(i));
			L_z(end+1) = (-x(i).*v_y(i) + y(i).*v_x(i));
			stars_vx(end+1) = v_x(i);%-velCoMx;		%velocity of stars within the range
			stars_vy(end+1) = v_y(i);%-velCoMy;
			stars_vz(end+1) = v_z(i);%-velCoMz;
			stars_mul(end+1) = mul(i);
			stars_mub(end+1) = mub(i);
			
			stars_x(end+1) = x(i);					%positions of stars inside the ranges
			stars_y(end+1) = y(i);
			stars_z(end+1) = z(i);
			tanV(end+1) = tan_vel(i);
			% if(met(i) >= min(met) && met(i) < max(met)/24)
				% colour = [colour;[0,0,0]];
				% metal(end+1) = met(i);
			% elseif(met(i) >= max(met)/24 && met(i) < max(met)*2/24)
				% colour = [colour;[0,0,50]];
				% metal(end+1) = met(i);
			% elseif(met(i) >= max(met)*2/24 && met(i) < max(met)*3/24)
				% colour = [colour;[0,0,100]];
				% metal(end+1) = met(i);
			% elseif(met(i) >= max(met)*3/24 && met(i) < max(met)*4/24)
				% colour = [colour;[0,0,150]];
				% metal(end+1) = met(i);
			% elseif(met(i) >= max(met)*4/24 && met(i) < max(met)*5/24)
				% colour = [colour;[0,0,200]];
				% metal(end+1) = met(i);
			% elseif(met(i) >= max(met)*5/24 && met(i) < max(met)*6/24)
				% colour = [colour;[0,0,250]];
				% metal(end+1) = met(i);
			% elseif(met(i) >= max(met)*6/24 && met(i) < max(met)*7/24)
				% colour = [colour;[0,50,250]];
				% metal(end+1) = met(i);
			% elseif(met(i) >= max(met)*7/24 && met(i) < max(met)*8/24)
				% colour = [colour;[0,100,250]];
				% metal(end+1) = met(i);
			% elseif(met(i) >= max(met)*8/24 && met(i) < max(met)*9/24)
				% colour = [colour;[0,150,250]];
				% metal(end+1) = met(i);
			% elseif(met(i) >= max(met)*9/24 && met(i) < max(met)*10/24)
				% colour = [colour;[0,200,250]];
				% metal(end+1) = met(i);
			% elseif(met(i) >= max(met)*10/24 && met(i) < max(met)*11/24)
				% colour = [colour;[0,250,250]];
				% metal(end+1) = met(i);
			% elseif(met(i) >= max(met)*11/24 && met(i) < max(met)*12/24)
				% colour = [colour;[0,250,200]];
				% metal(end+1) = met(i);
			% elseif(met(i) >= max(met)*12/24 && met(i) < max(met)*13/24)
				% colour = [colour;[0,250,150]];
				% metal(end+1) = met(i);
			% elseif(met(i) >= max(met)*13/24 && met(i) < max(met)*14/24)
				% colour = [colour;[0,250,100]];
				% metal(end+1) = met(i);
			% elseif(met(i) >= max(met)*14/24 && met(i) < max(met)*15/24)
				% colour = [colour;[0,250,50]];
				% metal(end+1) = met(i);
			% elseif(met(i) >= max(met)*15/24 && met(i) < max(met)*16/24)
				% colour = [colour;[0,250,0]];
				% metal(end+1) = met(i);
			% elseif(met(i) >= max(met)*16/24 && met(i) < max(met)*17/24)
				% colour = [colour;[50,250,0]];
				% metal(end+1) = met(i);
			% elseif(met(i) >= max(met)*17/24 && met(i) < max(met)*18/24)
				% colour = [colour;[100,250,0]];
				% metal(end+1) = met(i);
			% elseif(met(i) >= max(met)*18/24 && met(i) < max(met)*19/24)
				% colour = [colour;[150,250,0]];
				% metal(end+1) = met(i);
			% elseif(met(i) >= max(met)*19/24 && met(i) < max(met)*20/24)
				% colour = [colour;[200,250,0]];
				% metal(end+1) = met(i);
			% elseif(met(i) >= max(met)*20/24 && met(i) < max(met)*21/24)
				% colour = [colour;[250,250,0]];
				% metal(end+1) = met(i);
			% elseif(met(i) >= max(met)*21/24 && met(i) < max(met)*22/24)
				% colour = [colour;[250,200,0]];
				% metal(end+1) = met(i);
			% elseif(met(i) >= max(met)*22/24 && met(i) < max(met)*23/24)
				% colour = [colour;[250,150,0]];
				% metal(end+1) = met(i);
			% else
				% colour = [colour;[250,0,0]];
				% metal(end+1) = met(i);
			% endif
		endif
	endfor
endfor

%========ADDITIONAL PARAMETER CALCULATIONS==========

L2 = L_x.*L_x + L_y.*L_y + L_z.*L_z;	%total angular momentum (per unit mass) squared
L = sqrt(L2);
L_perp = sqrt(L_x.*L_x + L_y.*L_y);		%angular momentum components perpendicular to L_z



%=========PLOTTING================


% figure 1
% scatter(rad,L)
% xlabel("radius/kpc","FontSize",20)
% ylabel("L/kgkpckm/s","FontSize",20)

% figure 2
% scatter(rad,L_x)
% xlabel("radius/kpc","FontSize",20)
% ylabel("L_x/kgkpckm/s","FontSize",20)

% figure 3
% scatter(rad,L_y)
% xlabel("radius/kpc","FontSize",20)
% ylabel("L_y/kgkpckm/s","FontSize",20)

% figure 4
% scatter(rad,L_z)
% xlabel("radius/kpc","FontSize",20)
% ylabel("L_z/kgkpckm/s","FontSize",20)

figure 5
scatter(L_z, L_perp)
xlabel("L_z/kpckm/s","FontSize",20)
ylabel("L_perp/kpckm/s","FontSize",20)

% figure 6
% scatter(stars_vx, stars_vy)
% xlabel("v_x/km/s","FontSize",20)
% ylabel("v_y/km/s","FontSize",20)

% figure 7
% scatter(stars_vx, stars_vz)
% xlabel("v_x/km/s","FontSize",20)
% ylabel("v_z/km/s","FontSize",20)

% figure 8
% scatter(stars_vy, stars_vz)
% xlabel("v_y/km/s","FontSize",20)
% ylabel("v_z/km/s","FontSize",20)

% figure 9
% scatter3(stars_x, stars_y, stars_z)
% xlabel("X","FontSize",20)
% ylabel("Y","FontSize",20)
% zlabel("Z","FontSize",20)

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

% figure 13
% scatter3(rad,stars_mul,metal,"c")
% xlabel("heliocentric distance/kpc","FontSize",20)
% ylabel("mu_lcosb/mas/yr","FontSize",20)
% zlabel("metallicity","FontSize",20)

% figure 14
% scatter(stars_mul,metal)
% xlabel("mu_lcosb/mas/yr","FontSize",20)
% ylabel("metallicity","FontSize",20)

% figure 15
% scatter(stars_x, stars_y)
% xlabel("X","FontSize",20)
% ylabel("Y","FontSize",20)

% figure 16
% scatter(stars_x, stars_z)
% xlabel("X","FontSize",20)
% ylabel("Z","FontSize",20)

% figure 17
% scatter(stars_y, stars_z)
% xlabel("Y","FontSize",20)
% ylabel("Z","FontSize",20)


%============DENSITY PLOTS=================

a = rot90(L_z,-1);		%gets the transpose vector (column to row)
b = rot90(L_perp,-1);
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

figure 18				%tangential velocity-proper motion plot
im = imagesc(X,Y,density,[0 5])
set(gca(),"ydir","normal")
xlabel("L_z/kpckm/s")
ylabel("L_p/kpckm/s")
title("full density plot","FontSize",20)

figure 19
im_median = imagesc(X,Y,smooth)
set(gca(),"ydir","normal")
xlabel("L_z/kpckm/s")
ylabel("L_p/kpckm/s")
title("smoothed density plot","FontSize",20)

% figure 20
% im_mask = imagesc(X,Y,mask)
% set(gca(),"ydir","normal")
% xlabel("mu_lcosb/mas/yr","FontSize",20)
% ylabel("tangential velocity/km/s","FontSize",20)
% title("density - smooth","FontSize",20)

figure 21
im_mask_median = imagesc(X,Y,mask_median)
set(gca(),"ydir","normal")
xlabel("L_z/kgkpckm/s")
ylabel("L_p/kgkpckm/s")
title("median filtered mask density plot","FontSize",20)


% figure 23						% metallicity plot
% scatter(stars_mul,tanV,[],colour)
% xlabel("mu_lcosb/mas/yr","FontSize",20)
% ylabel("tangential velocity/km/s","FontSize",20)
% title("metallicity", "FontSize",20)

%======SAVE AND PRINT PLOTS AND FILES==========
% Outputs the density plots to screen and saves them to file.
%Converts to fits format for use with FINDCLUMPS
% Generates metallicity and infallz data in met_z and saves to file.
print(18,strcat("C:/Users/dpryd/Desktop/streamFind/Aquarius_plots_final/",filename,"_LpLz.png"),"-dpng")
print(21,strcat("C:/Users/dpryd/Desktop/streamFind/Aquarius_plots_final/",filename,"_LpLz_filtered.png"),"-dpng")
%print(23,strcat("C:/Users/dpryd/Desktop/streamFind/Aquarius_plots_final/",filename,"_LpLz_metallicity.png"),"-dpng")

save(strcat("C:/Users/dpryd/Desktop/streamFind/muVt/density_fits/",filename,"_LpLz_density"),"density")
save(strcat("C:/Users/dpryd/Desktop/streamFind/muVt/density_fits/",filename,"_LpLz_filtered"),"mask_median")

save_fits_image(strcat("C:/Users/dpryd/Desktop/streamFind/muVt/density_fits/",filename,"_LpLz_density.fits"), density)
save_fits_image(strcat("C:/Users/dpryd/Desktop/streamFind/muVt/density_fits/",filename,"_LpLz_filtered.fits"), mask_median)


