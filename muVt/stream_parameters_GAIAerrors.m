## Copyright (C) 2017 dpryd
## 
## This program is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
## 
## This script calculates the average position of the particles making up
## a stream in the Aquarius mock catalogues for use with the galpy stream integrator
## function. Also, velocity and infall z.
## Calculates energy of each particle with assumed grav. potential

clear s_x s_y s_z s_vx s_vy s_vz s_infallz s_l s_b s_pm_l s_pm_b s_Rhel s_Vlos s_R
stream_id = "157";
st = load(strcat("C:/Users/dpryd/Desktop/streamFind/mock_catalogues_final/AquariusA2_RRLS_pe_byIDstream/AquariusA2RRLS.id",stream_id,".pe.dat"));
st_info = load(strcat("C:/Users/dpryd/Desktop/streamFind/mock_catalogues_final/AquariusA2_RRLS_info/AquariusA2RRLS.id",stream_id,".info.dat"));

x_min = -75;	%limit the volume to calculate the average position for input to galpy
x_max = -50;
y_min = 0;
y_max = 50;
z_min = 20;
z_max = 40;
vrmin = -2;
vrmax = 2;

stream_x = st(:,19);	%position of particles in stream
stream_y = st(:,20);
stream_z = st(:,21);

stream_vx = st(:,22);	%velocities of particles in stream
stream_vy = st(:,23);
stream_vz = st(:,24);

stream_R = sqrt(stream_x.*stream_x + stream_y.*stream_y + stream_z.*stream_z);	%distance from centre of galaxy in galactocentric coords
stream_Rc = sqrt(stream_x.*stream_x + stream_y.*stream_y);	%galactocentric cylindrical distance
%==== OBSERVABLE QUANTITIES========
% galpy orbit integration can be initialised in observable quantities, which is
% more appropriate for the data provided from the simulations.

stream_l = st(:,26);	%galactic longitude
stream_b = st(:,27);	% galactic latitude
stream_pm_l = st(:,29);		% proper motion in l including cos b term
stream_pm_b = st(:,30);
stream_Rhel = st(:,28);		%heliocentric distance
stream_Vlos = st(:,31);		% line of sight/radial velocity
stream_parallax = st(:,25);
stream_Vlos_norm = stream_Vlos/220.0;
stream_R_norm = stream_R/8.0;
mu = sqrt(stream_pm_l.*stream_pm_l + stream_pm_b.*stream_pm_b);	%overall proper motion
tan_vel = 4.75*mu./stream_parallax;	%tangential velocity of each star
infallz = st_info(:,8);	%infall z for stream particles

for i=1:rows(stream_x)
	if(stream_x(i) > x_min && stream_x(i) < x_max)
		if(stream_y(i) > y_min && stream_y(i) < y_max)
			if(stream_z(i) > z_min && stream_z(i) < z_max)
				if(stream_Vlos(i) > vrmin && stream_Vlos(i) < vrmax)
					s_x(end+1) = stream_x(i);
					s_y(end+1) = stream_y(i);
					s_z(end+1) = stream_z(i);
					s_vx(end+1) = stream_vx(i);
					s_vy(end+1) = stream_vy(i);
					s_vz(end+1) = stream_vz(i);
					s_l(end+1) = stream_l(i);
					s_b(end+1) = stream_b(i);
					s_pm_l(end+1) = stream_pm_l(i);
					s_pm_b(end+1) = stream_pm_b(i);
					s_Rhel(end+1) = stream_Rhel(i);
					s_Vlos(end+1) = stream_Vlos(i);
					s_R(end+1) = stream_R(i);
				endif
			endif
		endif
	endif
endfor
	
	
avg_x = mean(s_x);	%average values
avg_y = mean(s_y);
avg_z = mean(s_z);
avg_vx = mean(s_vx);
avg_vy = mean(s_vy);
avg_vz = mean(s_vz);
%========= GALPY ORBIT INPUT PARAMETERS=======
disp("Inputs for galpy: ")
avg_l = mean(s_l)
avg_b = mean(s_b)
avg_pm_l = mean(s_pm_l)
avg_pm_b = mean(s_pm_b)
%avg_Rhel = mean(s_Rhel)
avg_R = mean(s_R)
avg_Vlos = mean(s_Vlos)
disp("Stream velocity parameters:")
disp("Total IDstream average radial velocity")
avg_Vlos_tot = mean(stream_Vlos)
disp("Total IDstream dispersion")
disp_Vlos = std(stream_Vlos)

figure 1
hold on
scatter3(stream_x,stream_y,stream_z)
scatter3(avg_x,avg_y,avg_z,"r","x")
xlabel("X/kpc")
ylabel("Y/kpc")
zlabel("Z/kpc")

% figure 2
% hold on
% scatter(stream_pm_l,stream_pm_b)
% scatter(avg_pm_l,avg_pm_b,"r","x")
% xlabel("mu_l/mas/yr")
% ylabel("mu_b/mas/yr")

figure 3
hold on
scatter(stream_R,stream_Vlos)
scatter(avg_R,avg_Vlos,"r","x")
xlabel("R_hel/kpc")
ylabel("V_rad/km/s")

% figure 4
% scatter(stream_Rhel,stream_z)
% xlabel("R_hel/kpc")
% ylabel("Z/kpc")

% figure 5
% scatter(stream_R_norm,stream_Vlos_norm)
% xlabel("R_hel/kpc")
% ylabel("V_rad/km/s")

% figure 6
% scatter(stream_R,stream_Vlos)
% xlabel("R/kpc")
% ylabel("V_r/km/s")

% figure 7
% scatter(stream_Rc,stream_Vlos)
% xlabel("R_c/kpc")
% ylabel("V_r/km/s")

figure 8
scatter(stream_pm_l,tan_vel)
xlabel("mu_l/mas/yr")
ylabel("V_t/km/s")


%======= ENERGY CALCULATION=========
%====== FULL STREAM=======
parsec = 3.086E16;	%metres
G = 6.67408E-11 %m3kg-1s-2
H = 73*1000/(3.086E22);
deltaV = 2.06E4;	#for AQA2
rhoC = 3*H^2/(8*pi*G);	%critical density
rho = 7.213*deltaV*rhoC;
R200 = 245.88*1000*parsec;	%in metres AQA2
c = 16.19;	%AQA2
Rs = R200/c;	%characteristic scale factor
stream_Rm = stream_R*1000*parsec;	%distance values in metres.
pe = ((-4*pi*G)*rho*(Rs^3).*log(1+stream_Rm/Rs))./stream_Rm;
vel2 = stream_vx.*stream_vx + stream_vy.*stream_vy + stream_vz.*stream_vz;
ke = 0.5*vel2;	%in (km/s)^2
E = ke + pe/(1E6);	%energy in (km/s)^2
disp("Average energy: ")
disp(mean(E))
disp("Max energy: ")
disp(max(E))
disp("Min energy: ")
disp(min(E))
%======PROGENITOR ENERGY==========
s_Rm = rot90(s_R,-1)*1000*parsec;	%distance values in metres.
pe_s = ((-4*pi*G)*rho*(Rs^3).*log(1+s_Rm/Rs))./s_Rm;
vel2_s = rot90(s_vx,-1).*rot90(s_vx,-1) + rot90(s_vy,-1).*rot90(s_vy,-1) + rot90(s_vz,-1).*rot90(s_vz,-1);
ke_s = 0.5*vel2_s;	%in (km/s)^2
E_s = ke_s + (pe_s/(1E6));	%energy in (km/s)^2
disp("Average energy of progenitor: ")
disp(mean(E_s))
disp("Max energy of progenitor: ")
disp(max(E_s))
disp("Min energy of progenitor: ")
disp(min(E_s))
