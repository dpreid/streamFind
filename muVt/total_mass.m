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
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.


## Author: dpryd <dpryd@HAL9000JNR>
## Created: 2017-06-28

##Uses the total mass of each simulation to associate each RRLS tracer
## a certain total overall mass, so that each ID stream can be associated
## a mass individually.

clear count;

%particle_mass = 1.37E4;	%AQA2 Solar masses
%M200 = 1.842E12;	%AQA2 solarmasses
%R200 = 245.88;		%AQA2 kpc
%particle_mass = 6.447E3; %AQB2
M200 = 8.194E11;	%AQB2
% R200 = 187.70;		%AQB2 kpc
%particle_mass = 1.399E4; %AQC2
%M200 = 1.774E12;	%AQC2
% R200 = 242.82;		%AQC2 kpc
%particle_mass = 1.397E4; %AQD2
%M200 = 1.774E12;	%AQD2
% R200 = 242.85;		%AQD2 kpc
% particle_mass = 9.593E3; %AQE2
%M200 = 1.185E12;	%AQE2
% R200 = 212.28;		%AQE2 kpc
m = 0.7;		%assumed mass of RRLS stars in Msun
num_dG = 0;		%number of dwarf galaxies
count=[];		%number of stars in each idstream
num_stars = 0;	%total number of stars
n_mp = 0;	%number of major progenitors (with more than 200 tracer particles)
star_count = [];	%number of stars in each progenitor
met_mean = [];
met_min = [];		%min values of metallicity in each progenitor
met_max = [];		% max values of metallicity in each progenitor
met_disp = [];		% standard deviation for metallicity
for j=1:N_prog
	d = load (file(j,:));
	di = load(file_info(j,:));
	met = di(:,5);	%metallicities, Z
	disp("stream:")
	disp(j)
	if(rows(d) >= 200)	%only get data for larger progenitors
		star_count = [star_count;rows(d)];
		met_mean = [met_mean;mean(met)];
		met_min = [met_min;min(met)];
		met_max = [met_max;max(met)];
		met_disp = [met_disp;std(met)];
	endif
	count(end+1) = rows(d);
	disp("Total number of stars:")
	num_stars += rows(d)
	disp("Furthest distance:")
	max(d(:,15))
	if(rows(d) >= 200)
		n_mp += 1;
	endif
endfor
disp("Mass associated with each RRLS star:")
mass = M200/num_stars
massIDs = count*mass;
massRRLS = num_stars*m;
disp("Fraction of mass in RRLS stars:")
fracRRLS = massRRLS/M200
for i=1:columns(count)
	if(massIDs(i) > 100000000)
		num_dG += 1;
	endif
endfor

massMet = [star_count,met_mean,met_min,met_max,met_disp];

disp("The number of streams with mass greater than 10E8 Msun:")
disp(num_dG)
disp("The number of major progenitors (>=200 tracer particles) =")
disp(n_mp)

%save(strcat("C:/Users/dpryd/Desktop/streamFind/muVt/stream_properties/",filename,"_massMet.dat"),"massMet")
figure 1
lmet_disp = log(met_disp);
lstar_count = log(star_count);
p=polyfit(lmet_disp,lstar_count,1)
x=linspace(min(lmet_disp),max(lmet_disp));
y=polyval(p,x);
scatter(lmet_disp,lstar_count)
hold on
plot(x,y)
title("AQE2 Star count - metallicity dispersion relation","FontSize",15)
xlabel("log(metallicity dispersion)","FontSize",20)
ylabel("log(star count)","FontSize",20)
hold off

figure 2
lmet_mean = log(met_mean);
[pp,s]=polyfit(lmet_mean,lstar_count,1)
w=linspace(min(lmet_mean),max(lmet_mean));
z=polyval(pp,w);
scatter(lmet_mean,lstar_count)
hold on
plot(w,z)
title("AQE2 Star count - metallicity relation","FontSize",15)
xlabel("log(mean metallicity)","FontSize",20)
ylabel("log(star count)","FontSize",20)
sqrt(diag(s.C)/s.df)*s.normr
dev=lstar_count-s.yf;
sterr = sqrt(sum(dev.*dev)/(rows(dev)-1))

