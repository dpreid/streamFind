## Copyright (C) 2016 dpryd
## 
## Requires phaseSpace_EAGLE.m to find the ID of the stars in the peak.
##Requires stream_ids (an array containing all the peak particle_IDs for the peak).
## Since there are a large number of stars and stream stars to loop through,
## the script first orders the stars by id number and then searches for the first stream star before starting 
## the search again from the previous star's index.
## RUN phaseSpace_EAGLE.m FIRST.

## After running this, run countPeak.m in order to count the number of particles
## returned to the progenitor.

clear s_x s_y s_z s_vx s_vy s_vz rad s_id_sort p_sort s_radVel s_mass metal;
particles = load("p302halo000.dat");
s_id_sort = sort(stream_ids);
p_sort = sortrows(particles,10);
x = p_sort(:,1);
y = p_sort(:,2);
z = p_sort(:,3);
v_x = p_sort(:,4);
v_y = p_sort(:,5);
v_z = p_sort(:,6);
met = p_sort(:,7);
m = p_sort(:,8);
id = p_sort(:,10);
r2 = x.*x + y.*y + z.*z;
radius = sqrt(r2);
n=size(x);
colours = [];
num = rows(stream_ids);
index=1;
for i=1:num
	for j=index:n
		if(s_id_sort(i) == id(j))
			rad(end+1) = radius(j);
			s_x(end+1) = x(j);
			s_y(end+1) = y(j);
			s_z(end+1) = z(j);
			s_vx(end+1) = v_x(j);
			s_vy(end+1) = v_y(j);
			s_vz(end+1) = v_z(j);
			s_mass(end+1) = m(i);
			metal(end+1) = met(i);
			% if (met(i) > 3.0)						%highest metallicity coloured RED
				% colours = [colours;[1 0 0]];
			% elseif ((met(i) > 2.0) & (met(i) <= 3.0))	% next highest metallicity coloured YELLOW
				% colours = [colours;[1 1 0]];
			% elseif ((met(i) > 1.0) & (met(i) <= 2.0))	% medium metallicity colour GREEN
				% colours = [colours;[0 1 0]];
			% elseif ((met(i) > 0.5) & (met(i) <= 1.0))
				% colours = [colours;[0 1 1]];			% medium-low metallicity coloured CYAN
			% elseif ((met(i) > 0.25) & (met(i) <= 0.5))
				% colours = [colours;[0 0 1]];			% second lowest metallicity coloured BLUE
			% else
				% colours = [colours;[0 0 0]];			% lowest metallicity coloured BLACK
			% endif
			% if(radius(j) < 100)
				% met_group_1 = [met_group_1;met(i)];
			% endif
			% if(radius(j) >= 100)
				% met_group_2 = [met_group_2;met(i)];
			% endif
			index = j+1;
			break;
		endif
	endfor
endfor

s_radVel = (s_vx.*s_x + s_vy.*s_y + s_vz.*s_z)./rad;
figure 1
scatter3(s_x,s_y,s_z)
xlabel("X/kpc","FontSize",20)
ylabel("Y/kpc","FontSize",20)
zlabel("Z/kpc","FontSize",20)
title("halo015, peak particles only, z = 0.178","FontSize",20)
axis equal
figure 2
scatter(rad,s_radVel)
xlabel("radius/kpc","FontSize",20)
ylabel("radial velocity/km/s","FontSize",20)
title("halo015 peak particles only, z=0.178","FontSize",20)
print(1,strcat("C:/Users/dpryd/Desktop/streamFind/EAGLE_plots/",filename,"_xyz_peak5.png"),"-dpng")
print(2,strcat("C:/Users/dpryd/Desktop/streamFind/EAGLE_plots/",filename,"_rVr_peak5.png"),"-dpng")