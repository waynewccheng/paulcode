%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compare WSI to Truth using the CIELAB from measurements
%%%%%

% 03-29-2020: First version

clearvars;
close all;

% File names 
imgTruthName = 'truth.tif';
labTruthName = 'LAB_array';
imgTestName = 'BiomaxOrgan10_UterineCervix_B10_13_5x_reg.tif';
transName = 'trans_mean_camera';

% Load files
Truth = imread(imgTruthName);
RegImg = imread(imgTestName);
temp = load(labTruthName);
LAB_array = temp.LAB_array;
temp = load(transName);
trans_array_m = temp.trans_array_m;
clearvars temp;

%% Displays Truth and WSi image
figure('units','normalized','OuterPosition',...
     [0.2453125 0.178703703703704 0.5484375 0.45]);
subplot(1, 2, 1);
imshow(Truth);
xlabel('Truth');

subplot(1, 2, 2);
imshow(RegImg);
xlabel('WSI registered img');
print('Truth_vs_WSI','-dpng');

%% Heatmap and boxplot using the HIMS sRGB values
diffMatrix = compareColorFunction_LAB(imgTruthName, imgTestName); % From Sam's code

figure;
imagesc(diffMatrix);
colorbar;
title('dE (Truth CIELAB based on sRGB values)')
print('dE_heatmap','-dpng');

diffMatrix_rs = reshape(diffMatrix, size(diffMatrix, 1) * size(diffMatrix, 2), 1);
figure; 
boxplot(diffMatrix_rs);
title(sprintf('dE (Truth CIELAB based on sRGB values), Median = %.2f', median(diffMatrix_rs)))
print('boxplot','-dpng');

%% Heatmap and boxplot using the HIMS CIELAB values
diffMatrix = compareColorFunction_LAB2(imgTestName, LAB_array); % Modif of Sam's code to accomodate using the CIELAB measurements

figure;
imagesc(diffMatrix);
colorbar;
title('dE (Truth CIELAB based on transmittance values)')
print('dE_2_heatmap','-dpng');

diffMatrix_rs = reshape(diffMatrix, size(diffMatrix, 1) * size(diffMatrix, 2), 1);
figure; 
boxplot(diffMatrix_rs);
title(sprintf('dE (Truth CIELAB based on transmittance values), Median = %.2f', median(diffMatrix_rs)))
print('boxplot_2','-dpng');

%% Principal Component Analysis
% Following https://stats.stackexchange.com/questions/134282/relationship-between-svd-and-pca-how-to-use-svd-to-perform-pca
b_s = 2; % First wavelength is 390nm = first band, no signal at 380 nm
data = trans_array_m(b_s:end-1, :)';
b_f = b_s + size(data, 2) - 1; % last position spectrum is 770nm = last band, no signal at 780 nm
n = size(data,1);

data_mean = mean(data,2);   % find mean
A = data - data_mean;       % subtract mean from each point
rho = norm(A,'fro')^2;      % total variation of data, norm is sqrt(trace(A_t A)) with A centered

[U,S,V] = svd(A,'econ');    % find singular value decomposition
sigma = diag(S);            % singular values

% Reconstruct original data, limiting the number of eigenvectors
nc = 5;
q = norm(sigma(1:nc))^2/rho;
% prct_var = 100*sum(sigma(1:nc).^2)/sum(sigma.^2); Same as q
A_nc = U(:, 1:nc)*S(1:nc, 1:nc)*V(:, 1:nc)'; 

% figure;
% step = 10000;
% for i = 1:step:size(A_nc, 1)
%     plot(A_nc(i, :), 'r'); hold on; 
%     plot(A(i, :), 'b');
% end
% title(sprintf('Transmittance, %u components, captures %.4g%% of total variation (%% var)',nc, 100*q));

data_nc = A_nc + data_mean;
lambda = 380:10:780;
figure;
step = 10000;
for i = 1:step:size(data_nc, 1)
    plot(lambda(b_s: b_f), data_nc(i, :), 'r'); hold on; 
    plot(lambda(b_s: b_f), data(i, :), 'b');
end
title(sprintf('Transmittance, %u components, captures %.4g%% of total variation (%% var)',nc, 100*q));
print('T_after_PCA','-dpng');

i = 11000;
% figure;
% plot(A_nc(i, :), 'r'); hold on;
% plot(A(i, :), 'b');
% title(sprintf('%u components, captures %.4g%% of total variation (prct variance)',nc, 100*q));

figure;
plot(lambda(b_s: b_f),data_nc(i, :), 'r'); hold on;
plot(lambda(b_s: b_f),data(i, :), 'b');
title(sprintf('T (pixel %u), %u components, captures %.4g%% of total var (%% var)', i, nc, 100*q));
print('T_one_pix_after_PCA','-dpng');

% Mean data vector projected on eigenvector base
spectral_mean = mean(data, 1);

b_nc = spectral_mean*V(:, 1:nc);
figure; 
plot(lambda(b_s: b_f),spectral_mean, 'b'); hold on; 
plot(lambda(b_s: b_f),mean(spectral_mean) + b_nc*V(:, 1:nc)', 'r');
title(sprintf('Spectral mean, %u components, captures %.4g%% of total var (%%var)',nc, 100*q));
print('Spectral_mean_after_PCA','-dpng');

% Plot the 5 first eigenvectors with the meand spectra
figure('units','normalized','OuterPosition',...
    [0.2453125 0.178703703703704 0.5484375 0.821296296296296]);

for i = 1:5
    subplot(2, 3, i);
    plot(lambda(b_s: b_f),V(:, i)', 'b');
    title(sprintf('Eigenvector %u', i));
    
end
subplot(2, 3, 6);
plot(lambda(b_s: b_f),spectral_mean-mean(spectral_mean), 'b');
title('Spectral mean');
print('Eigenvector_SpectralMean','-dpng');

%% Convert transmittance data to CIEXYZ, CIELAB
    
% Prepares the illuminant
load ('C:\Users\Paul_Lemaillet\Documents\WholeSlideImaging\ProgMatlab\DataIlluminants\spec_cied65','spec');
ls = spec(1:10:401,2);

% Trans -> XYZ -> LAB with saving XYZ and LAB (values + uncertainties)
sizey = size(Truth, 1);
sizex = size(Truth, 2);
[LAB_array, ~, XYZ_array, ~] = f_transmittance2LAB(trans_array_m, [], sizey, sizex, ls, 'y'); % 'y' top trim the max tranmsittance to 1
[LAB_array_nc, ~, XYZ_array_nc, ~] = f_transmittance2LAB_pca(data_nc', [], sizey, sizex, ls, b_s, b_f, 'y'); % 'y' top trim the max tranmsittance to 1
    
% Based on Sam's code, I did it fast and dirty (no function)
diff2D = (LAB_array - LAB_array_nc).^2;
diff1D = sqrt(sum(diff2D,2));
diffMatrix = reshape(diff1D, sizey, sizex);
save('dE_after_PCA.mat','diffMatrix');
 
figure;
imagesc(diffMatrix);
colorbar;
title('dE for Truth after PCA')
print('dE_heatmap_after_PCA','-dpng');

diffMatrix_rs = reshape(diffMatrix, sizex * sizey, 1);
figure; 
boxplot(diffMatrix_rs);
title(sprintf('dE for Truth after PCA, Median = %.2f', median(diffMatrix_rs)))
print('boxplot_after_PCA','-dpng');

% Heatmap for each of the 5 first eigenvectors
for i = 1:5
    A_nc = U(:, i)*S(i, i)*V(:, i)';
    data_nc = A_nc + data_mean;
    [LAB_array_nc, ~, XYZ_array_nc, ~] = f_transmittance2LAB_pca(data_nc', [], sizey, sizex, ls, b_s, b_f, 'y'); % 'y' top trim the max tranmsittance to 1
    diffMatrix(:, :, i) = compareColorFunction_LAB2(imgTestName, LAB_array_nc);
end

figure('units','normalized','outerposition',[0 0 1 1]);
for i = 1:5
    subplot(2, 3, i);
    imagesc(diffMatrix(:, :, i));
    colorbar;
    title(sprintf('Eigenvector %u', i));
end
sgtitle('dE of Truth vs WSI');
print('dE_heatmap_eigenvector_1_5','-dpng');

figure('units','normalized','outerposition',[0 0 1 1]);
for i = 1:5
    subplot(2, 3, i);
    diffMatrix_rs = reshape(diffMatrix(:, :, i), sizex * sizey, 1);
    boxplot(diffMatrix_rs);
    title(sprintf('Eigenvector %u, Median = %.2f', i, median(diffMatrix_rs)));
end
print('boxplot_eigenvector_1_5','-dpng');