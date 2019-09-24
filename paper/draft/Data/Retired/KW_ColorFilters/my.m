% Results Format:
%
% Names_KW_ColFilters: name of data folder set during measurements
%
% Lambda: 380nm to 780nm in steps of 10nm
%
% T_Spectro_m: temporal mean over 10 measured spectra
% 41 wavelengths x 5 color filters
%
% T_Spectro_s: temporal std dev over 10 measured spectra, total uncertainty
% 41 wavelengths x 5 color filters
%
% T_Spectro_cmp: T Spectro data, Lambda then (value, 1 X sigma, 2 X sigma) x 5 (15 col)
%
% Cam_m: camera spatial averages, mean values over 10 measured spectra
% 41 wavelengths x 5 color filters
%
% Cam_s: camera spatial averages, std dev values over 10 measured spectra, total uncertainty
% 41 wavelengths x 5 color filters
%
% Cam_T_cmp: cam data, Lambda then (value, 1 X sigma, 2 X sigma) x 5 (15 col)
%
% LAB_spectro_Kodak_ND, LAB_cam_Kodak_ND: 5 color filters x 15 columns
% col 1:3: L, a, b
% col 4:6: uncertainty on L, a, b (repetability)
% col 7:9: reproducibility uncertainties on L, a, b
% col 10:12: total uncertainty (sqrt of sum of square of col 4:6 and col 7:9)
% col 13:15: expanded uncertainty: 2 * col 10:12
%
% DE_Kodak_ND: 5 color filters x 5 columns
% col 1: Delta E
% col 2: uncertainty on Delta E
% col 3: reproducibility uncertainties on Delta E
% col 4: total uncertainty (sqrt of sum of square of col 2 and col 3)
% col 5: expanded uncertainty: 2 * col 4

lambda = csvread('lambda.txt')
Cam_T_m = csvread('Cam_T_m.txt')
Cam_T_s = csvread('Cam_T_s.txt')
Cam_T_cmp = csvread('Cam_T_cmp.txt')
T_Spectro_m = csvread('T_Spectro_m.txt')
T_Spectro_s = csvread('T_Spectro_s.txt')
T_Spectro_cmp = csvread('True_T_cmp.txt')
LAB_cam_KW_ColFilters = csvread('LAB_cam_KW_ColFilters.txt')
LAB_spectro_KW_ColFilters = csvread('LAB_spectro_KW_ColFilters.txt')

fid = fopen('Names_KW_ColFilters.txt');
txt = textscan(fid,'%s','delimiter','\n'); 
fclose(fid)
str = txt{1}

% figure
hold on
for filter_id = 1:5
    plot(lambda,Cam_T_cmp(:,2+3*(filter_id-1)))
end
legend(str')
xlabel('Wavelength (nm)')
ylabel('Transmittance')
title('Kodak Wratten Color Filters')
axis([380  780 0 1])