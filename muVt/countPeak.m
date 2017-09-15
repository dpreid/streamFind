% FIRSTLY count the number of particles which were part of the peak
% which return to the progenitor. Also find the number of additional, non-progenitor
% particles which were within the peak as well.

n=columns(rad);
countProg = 0;	%count of stars returned to same progenitor
countNProg = 0;	% count of other particles within peak
for i=1:n
	if(rad(i) > 100) %s_radVel(i) > rad(i) - 200 & rad(i) > 50)
		countProg +=1;
	else
		countNProg +=1;
	endif
endfor

disp("number of stars returned to progenitor")
disp(countProg)
disp("number of additional stars within peak")
disp(countNProg)
disp("sum")
disp(countNProg+countProg)

% NEXT: find the number of particles within the original progenitor
countTotal = 0;
%data = load ("p302halo000.dat");
rad_particles = [p_sort,radius];
%rad_sort = sortrows(rad_particles,11);	%sort the data by radius
con1 = rad_particles(:,11) < 140;	%remove the particles well outwith the range of the progenitor
rad_particles(con1,:) = [];
con2 = rad_particles(:,11) > 160;
rad_particles(con2,:) = [];
x= rad_particles(:,1);		%spatial positions
y = rad_particles(:,2);
z = rad_particles(:,3);
v_x = rad_particles(:,4);	%velocities
v_y = rad_particles(:,5);
v_z = rad_particles(:,6);
rad = rad_particles(:,11);
radVel = (v_x.*x + v_y.*y + v_z.*z)./rad;

for i=1:rows(rad)
	if(rad(i) > 140) %radVel(i) > rad(i)-200)
		countTotal += 1;
	endif
endfor
disp("Total number of stars in original progenitor")
disp(countTotal)
		

figure 1
scatter(rad,radVel)
xlabel("distance/kpc","FontSize",20)
ylabel("radial velocity/km/s","FontSize",20)
title("p302halo000 all particles","FontSize",20)
print(1,strcat("C:/Users/dpryd/Desktop/streamFind/EAGLE_plots/",filename,"_rVrFULL.png"),"-dpng")
%hold on
%plot(rad,rad-200,"r")
%hold off

figure 2
scatter3(x,y,z)
xlabel("x/kpc","FontSize",20)
ylabel("y/kpc","FontSize",20)
zlabel("z/kpc","FontSize",20)
title("p302halo000 all particles","FontSize",20)

