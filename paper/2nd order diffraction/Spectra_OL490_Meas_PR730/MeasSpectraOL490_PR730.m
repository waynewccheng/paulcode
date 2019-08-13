%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display spectra of Ol490 light source %
%

close all;

% Load data
load('spectra','spectra');
load('spectra_m','spectra_m');
load('spectra_s','spectra_s');

% Graphics
figure('units','normalized','outerposition',[0 0 1 1]);
lambda = 380:780;
wl_array = 380:10:780;
for i = 1:41
    subplot(6,7,i)
    plot(lambda, spectra_m(i, :))
    title(sprintf('%d',wl_array(i)))
end

figure('units','normalized','outerposition',[0 0 1 1]);
wl_array = 380:10:780;
for i = 1:41
    subplot(6,7,i)
    errorbar(lambda, spectra_m(i, :), spectra_s(i, :))
    title(sprintf('%d',wl_array(i)))
end
