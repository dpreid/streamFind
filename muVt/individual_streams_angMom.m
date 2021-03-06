## Author: dpryd <dpryd@HAL9000JNR>
## Created: 2017-05-24
## Copyright (C) 2016 dpryd
## 
## This script calculates and plots multiple phase spaces
## This version is for the mock catalogues from the Aquarius and EAGLE simulations.
## 
## 
## Currently using no errors data
##
## USES PROPER MOTION VELOCITIES



%======INITIALISATION=======
%filename = "AQA2RRLS_r10_100_noerrors";
N_prog = 54;		% SET IN pmTanV_init.m number of files in the simulation, corresponding to the number of progenitor systems
nbins=600;			%resolution of the density plots
Ms = 1.989E30;	% solar mass in kg
mass = 0.7*Ms;	%mass of each RRLS particle

r_min = 20;		%min and max distance for plotting
r_max = 100;



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

for j=54:N_prog
	clear tanV metal stars_mul stars_mub L_x L_y L_z L_perp L2 L rad colour stars_vx stars_vy stars_vz stars_x stars_y stars_z starsA_x starsA_y starsA_z;
	disp("Plotting stream:")
	disp(j)
	colour = [];
	stars_mul = [];
	data = load (file(j,:));		%load the required halo file
	%info = load(file_info(j,:));	% load the corresponding info file with metallicity etc
	stream = file(j,:);
	x= data(:,6);		%spatial positions
	y = data(:,7);
	z = data(:,8);
	v_x = data(:,9);	%velocities/km/s
	v_y = data(:,10);
	v_z = data(:,11);
	%met = info(:,5);	%metallicities, Z
	%met_av = mean(met);	%mean metallicity of the progenitor 
	%infallz = info(:,8);	%redshift when the progenitor is considered to have fallen into the host galaxy.
	mul = data(:,16);	%proper motion in galactic longitude, mulcosb
	mub = data(:,17);	% proper motion in galactic latitude
	
	rhel = data(:,15);	%heliocentric distance/kpc
	parallax = data(:,12); % parallax
	mu = sqrt(mul.*mul + mub.*mub);	%overall proper motion
	tan_vel = 4.75*mu./parallax;	%tangential velocity of each star
	n=size(x);			%number of particles in each progenitor
	for i=1:n
		if(rhel(i) > r_min && rhel(i) < r_max)	%only plot particles within these spatial ranges
			rad(end+1) = rhel(i); 
			stars_mul(end+1) = mul(i);
			stars_mub(end+1) = mub(i);
			L_x(end+1) = (-y(i).*v_z(i) + z(i).*v_y(i));	%angular momentum components per unit mass
			L_y(end+1) = (-z(i).*v_x(i) + x(i).*v_z(i));
			L_z(end+1) = (-x(i).*v_y(i) + y(i).*v_x(i));
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
	
	
	
	if(rows(stars_mul) > 0)
		%L = sqrt(L_x.*L_x + L_y.*L_y + L_z.*L_z);	%total angular momentum (per unit mass) squared
		L_perp = sqrt(L_x.*L_x + L_y.*L_y);		%angular momentum components perpendicular to L_z
		a = rot90(L_z,-1);		%gets the transpose vector (column to row)
		b = rot90(L_perp,-1);
		M = [b,a];
		density = hist3(M,[nbins,nbins]);	%count data stored in density
		[X,Y] = meshgrid(linspace(min(a),max(a)),linspace(min(b),max(b)));	%meshgrid for plotting the colourmap

	%============DENSITY PLOTS=================
		figure (j)				%tangential velocity-proper motion plot
		im = imagesc(X,Y,density,[0 5])
		set(gca(),"ydir","normal")
		xlabel("L_z/kpckm/s","FontSize",20)
		ylabel("L_perp/kpckm/s","FontSize",20)
		title("full density plot","FontSize",20)
		
		
		
		print(j,strcat(ind_stream_file, sprintf("%03d",j),"_LpLz_density.png"),"-dpng")
	endif
endfor
