## Copyright (C) 2016 dpryd
## Creates the list of filenames for reading into rVrClusterFind2.m 
#Currently for Aquarius A2 simulation
% missing_files = [18,58,102,115,119,126,127,133,139,141,149,176];	%for AQA2
% n_files = 177;	%for AQA2
##missing_files = [2];	%for hydro008
##n_files = 26;			 %for hydro008
##missing_files = [];	%for hydro009
##n_files = 32;			 %for hydro009
##missing_files = [];	%for AQD2
##n_files = 127;			 %for AQD2
##n_files = 78;			 %for AQB2
##missing_files = [1,3,4,41,58,78];	%for AQE2
##n_files = 91;			 %for AQE2
missing_files = [];
n_id = 177;		%for AQA2RRLS

% file = [];
% file_path_hydro = "D:/MSc_Astrophysics/MSci_analysis/GAIA_mock_catalogues/HYDRO_ZOOM/HYDRO009_accreted/";
% file_path_aq = "C:/Users/dpryd/Documents/MSci_analysis/GAIA_mock_catalogues/Aquarius/Aquarius_A2_RRLS/";
% for i=1:n_files
	% if(any(missing_files==i)==0)
		% if(i < 10)
			% name = strcat(file_path_aq, "AquariusA2RRLS (",num2str(i),").dat");
		% endif
		% if(i >= 10 && i < 100)
			% name = strcat(file_path_aq, "AquariusA2RRLS (",num2str(i),").dat");
		% endif
		% if(i > 100)
			% name = strcat(file_path_aq, "AquariusA2RRLS (",num2str(i),").dat");
		% endif
	% file = [file;name];	
	% endif
% endfor

file = [];
file_path_hydro = "D:/MSc_Astrophysics/MSci_analysis/GAIA_mock_catalogues/HYDRO_ZOOM/HYDRO009_accreted/";
file_path_aq = "C:/Users/dpryd/OneDrive/MSc_astrophysics_data/mock_catalogues_final/AquariusA2_RRLS_pe_byIDstream/";
% for i=1:n_files
	% if(any(missing_files==i)==0)
		% if(i < 10)
			% name = strcat(file_path_aq, "AquariusA2RRLS.id00",num2str(i),".pe.dat");
		% endif
		% if(i >= 10 && i < 100)
			% name = strcat(file_path_aq, "AquariusA2RRLS.id0",num2str(i),".pe.dat");
		% endif
		% if(i > 100)
			% name = strcat(file_path_aq, "AquariusA2RRLS.id",num2str(i),".pe.dat");
		% endif
	% file = [file;name];	
	% endif
% endfor

for i=1:n_id
	name = strcat(file_path_aq, "AquariusA2RRLS.id", sprintf("%3d",i),".pe.dat");
	if(exist(name) == 2)
		file = [file;name];
	endif
endfor