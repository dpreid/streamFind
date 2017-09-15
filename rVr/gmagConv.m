## Copyright (C) 2016 David
## 
## The plot produced by this program converts GAIA G magnitudes into Grvs (radial velocity spectrometer).
## Based upon the conversion tables (Table 3) in the paper GAIA broadband photometry Jordi et al 2010.
##
## Use: run this program before running rVrClusterFind4.m. The polynomial conversion is
## found here and then used in rVrClusterFind4 in order to convert G-V into Grvs

## A different polynomial is used for V-Ic < 0 and V-Ic > 0 else we have two possible values
## for some G-V


xneg = [-10:0.01:0];	% x is V-Ic Johnson-Cousins colour magnitude
xpos = [0:0.01:2];		% V-Ic > 0
GVneg = -0.0257-0.0924*xneg-0.1623*xneg.*xneg+0.0090*xneg.*xneg.*xneg;	%G-V colour magnitude
GGrvsneg = -0.0138+1.1168*xneg-0.1811*xneg.*xneg+0.0085*xneg.*xneg.*xneg;	%G-Grvs colour magnitude
GVpos = -0.0257-0.0924*xpos-0.1623*xpos.*xpos+0.0090*xpos.*xpos.*xpos;	%G-V colour magnitude
GGrvspos = -0.0138+1.1168*xpos-0.1811*xpos.*xpos+0.0085*xpos.*xpos.*xpos;	%G-Grvs colour magnitude

pneg = polyfit(GVneg,GGrvsneg,5);	% finds the polynomial of order n which best fits the data
ppos = polyfit(GVpos,GGrvspos,5);
fneg = polyval(pneg,GVneg);			%calculates the value of the polynomial at specific values
fpos = polyval(ppos,GVpos);
hold on
%scatter(GVneg,GGrvsneg)
%scatter(GVneg,fneg,"r")
scatter(GVpos,GGrvspos)
scatter(GVpos,fpos,"r")
