## Copyright (C) 2016 dpryd
## plots radial velocity and radius with inclination and metallicity as additional parameters.
## v3 includes code to identify stars in elliptical bins (ie picking out tidal streams)
## Can also plot stars within a straight bin in order to capture other clusters in phase space.

data = load ("p178halo015.dat");
clear colours rad r_v inc theta metal stream_x stream_y stream_z stream_vr stream_rad stream_mass stars_id;

%===== ELLIPSE PARAMETERS=====
vrmax = 240;
rmax = 42;
vrmax2 = 260;
rmax2 = 48;
shift = 0; % positive for upwards shift
a = 1;
b1 = (vrmax^2)/(rmax^2);
b2 = (vrmax2^2)/(rmax2^2);

%=================================
% ======== STRAIGHT LINE BIN PARAMETERS==========
m = 1.6;
c1 = -180;
c2 = -250;
%=================================================
x = data(:,1);		%spatial positions
y = data(:,2);
z = data(:,3);
v_x = data(:,4);	%velocities
v_y = data(:,5);
v_z = data(:,6);
met = data(:,7);	%metallicity of each star relative to sun's metallicity
mass = data(:,8);
r2 = x.*x + y.*y + z.*z;	%square radius
radius = sqrt(r2);
n=size(x);
vel2 = v_x.*v_x + v_y.*v_y + v_z.*v_z;
vel = sqrt(vel2);
radVel = (v_x.*x + v_y.*y + v_z.*z)./radius;
theta = acos(abs(z./radius));	%inclination to z axis
colours = [];
stream_colour = [];
id = data(:,10);
for i=1:n
	if ((radius(i) > 20) & (radius(i) < 250))	%plotting conditions
		rad(end+1) = radius(i);
		r_v(end+1) = radVel(i);
		metal(end+1) = met(i);
		if (met(i) > 3.0)						%highest metallicity coloured RED
			colours = [colours;[1 0 0]];
		elseif ((met(i) > 2.0) & (met(i) <= 3.0))	% next highest metallicity coloured YELLOW
			colours = [colours;[1 1 0]];
		elseif ((met(i) > 1.0) & (met(i) <= 2.0))	% medium metallicity colour GREEN
			colours = [colours;[0 1 0]];
		elseif ((met(i) > 0.5) & (met(i) <= 1.0))
			colours = [colours;[0 1 1]];			% medium-low metallicity coloured CYAN
		elseif ((met(i) > 0.25) & (met(i) <= 0.5))
			colours = [colours;[0 0 1]];			% second lowest metallicity coloured BLUE
		else
			colours = [colours;[0 0 0]];			% lowest metallicity coloured BLACK
		endif
	endif
	%=========for ELLIPTICAL BINS, include the following code=====================
	% if(radius(i) > 20 & (r2(i) > rmax^2 -(radVel(i)-shift)*(radVel(i)-shift)/b1) & (r2(i) < rmax2^2 -(radVel(i)-shift)*(radVel(i)-shift)/b2))	%condition to be inside elliptical bin
	%if((radius(i) > 0) & (radius(i) < 200))	
		% stars_id(end+1) = id(i);
		% stream_x(end+1) = x(i);
		% stream_y(end+1) = y(i);
		% stream_z(end+1) = z(i);
		% stream_vx(end+1) = v_x(i);
		% stream_vy(end+1) = v_y(i);
		% stream_vz(end+1) = v_z(i);
		% stream_vr(end+1) = radVel(i);
		% stream_rad(end+1) = radius(i);
		% stream_mass(end+1) = mass(i);
		% if(radVel(i) > 0)
			% stream_colour = [stream_colour;[1 0 0]];	%moving away from centre coloured red
		% else
			% stream_colour = [stream_colour;[0 0 1]];	% moving towards centre coloured blue
		% endif
	% endif
	%========== for STRAIGHT BIN include following code ===============
	if((radius(i) > 20) & (radius(i) < 200) & (radVel(i) < radius(i)*m+c1) & (radVel(i) > radius(i)*m+c2)) % conditions to be inside straight gradient bin
		stars_id(end+1) = id(i);
		stream_x(end+1) = x(i);
		stream_y(end+1) = y(i);
		stream_z(end+1) = z(i);
		stream_vx(end+1) = v_x(i);
		stream_vy(end+1) = v_y(i);
		stream_vz(end+1) = v_z(i);
		stream_vr(end+1) = radVel(i);
		stream_rad(end+1) = radius(i);
		stream_mass(end+1) = mass(i);
	endif
endfor

figure 1
scatter(rad,r_v,[],colours)
xlabel("radius/kpc","FontSize",20)
ylabel("radial velocity/km/s","FontSize",20)
title("0p000halo024, Z > 1, theta > 3pi/8","FontSize",20)
% figure 2
% scatter(stream_rad,stream_vr)
% xlabel("radius/kpc","FontSize",20)
% ylabel("radial velocity/km/s","FontSize",20)
% title("0p000halo031, stream only","FontSize",20)
% figure 3
% scatter3(stream_x,stream_y,stream_z,[],stream_colour)
% xlabel("X/kpc","FontSize",20)
% ylabel("Y/kpc","FontSize",20)
% zlabel("Z/kpc","FontSize",20)
% title("halo031, tidal stream plot","FontSize",20)
% figure 4
% scatter(stream_vx,stream_vz)
% xlabel("V_X/km/s","FontSize",20)
% ylabel("V_Z/km/s","FontSize",20)
% title("halo031, tidal stream plot","FontSize",20)
% figure 5
% scatter(stream_vx,stream_vy)
% xlabel("V_X/km/s","FontSize",20)
% ylabel("V_Y/km/s","FontSize",20)
% title("halo031, tidal stream plot","FontSize",20)