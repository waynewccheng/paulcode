%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display spectra of Ol490 light source %
%

close all;

% Load data
load('lambda','lambda');
load('spectra','spectra');
load('spectra_m','spectra_m');
load('spectra_s','spectra_s');

% Graphics
wl_range = 380:390;
n_lambda = size(wl_range, 2);
n_rows = double(uint8(sqrt(n_lambda)));
n_cols = double(uint8(n_lambda/floor(sqrt(n_lambda))));

figure('units','normalized','outerposition',[0 0 1 1]);
for i = 1:n_lambda 
    subplot(n_rows,n_cols,i)
    plot(lambda, spectra_m(i, :))
    title(sprintf('%d',wl_range(i)))
end

figure('units','normalized','outerposition',[0 0 1 1]);
for i = 1:n_lambda 
    subplot(n_rows,n_cols,i)
    errorbar(lambda, spectra_m(i, :), spectra_s(i, :))
    title(sprintf('%d',wl_range(i)))
end
