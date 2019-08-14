%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display spectra of Ol490 light source %
%
% Question: does OL490 generate 2nd order diffraction at 2x frequency?
% Method: use PR730 to detect any peaks in [380,780] that, generate stimulus at
% [380,780]/2
% Result: another peak shown at 2x wavelength with 50% power

close all;

% Load data
load('lambda','lambda');
load('spectra','spectra');
load('spectra_m','spectra_m');
load('spectra_s','spectra_s');

% 2nd order analysis
for i = 1:n_lambda 

    % get peaks
    spd = spectra_m(i, :);
    spd_first_half = spd(1,[1:200]);
    spd_second_half = spd(1,[201:end]);
    
    max1 = max(spd_first_half);
    max2 = max(spd_second_half);
    
end

% Graphics
wl_range = 380:390;
n_lambda = size(wl_range, 2);
n_rows = double(uint8(sqrt(n_lambda)));
n_cols = double(uint8(n_lambda/floor(sqrt(n_lambda))));

figure('units','normalized','outerposition',[0 0 1 1]);
for i = 1:n_lambda 

    % get peaks
    spd = spectra_m(i, :);
    spd_first_half = spd(1,[1:200]);
    spd_second_half = spd(1,[201:end]);
    
    max1 = max(spd_first_half);
    max2 = max(spd_second_half);
    
    subplot(n_rows,n_cols,i)
    hold on
    plot(lambda, spectra_m(i, :),'-b')
    
    % predict 2nd order diffraction
    % use only 380:1:390
    % should appear at 780:2:780
    % divided by 2
    plot(2*lambda(1:11),spectra_m(i,1:11)/2,'or')
    
    title(sprintf('%d nm',wl_range(i)))
    xlabel('Wavelength (nm)')
    ylabel('Power')
    
end

snapnow

return

figure('units','normalized','outerposition',[0 0 1 1]);
for i = 1:n_lambda 

    % get peaks
    spd = spectra_m(i, :);
    spd_first_half = spd(1,[1:200]);
    spd_second_half = spd(1,[201:end]);
    
    max1 = max(spd_first_half);
    max2 = max(spd_second_half);
    
    subplot(n_rows,n_cols,i)
    bar(lambda, spectra_m(i, :))
    title(sprintf('%d nm - %f %f %f',wl_range(i),max1/max2,max1,max2))
    xlabel('Wavelength (nm)')
    ylabel('Power')
end

return

figure('units','normalized','outerposition',[0 0 1 1]);
for i = 1:n_lambda 
    subplot(n_rows,n_cols,i)
    errorbar(lambda, spectra_m(i, :), spectra_s(i, :))
    title(sprintf('%d nm',wl_range(i)))
end
