%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compare WSI to Truth using the CIELAB from measurements
%%%%%

% 03-29-2020: First version

clearvars;
close all;

% File names 
imgTruthName = 'truth.tif';
labTruthName = 'LAB_array';
imgTestName = 'BiomaxOrgan10_Liver_H9_13_5x_reg.tif';
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
save('dE_imgs.mat','diffMatrix');

figure;
imagesc(diffMatrix);
colorbar;
title('dE (Truth CIELAB based on sRGB values)')
print('dE_heatmap_imgs','-dpng');

diffMatrix_rs = reshape(diffMatrix, size(diffMatrix, 1) * size(diffMatrix, 2), 1);
figure; 
boxplot(diffMatrix_rs);
title(sprintf('dE (Truth CIELAB based on sRGB values), Median = %.2f', median(diffMatrix_rs)))
print('boxplot_imgs','-dpng');

%% Heatmap and boxplot using the HIMS CIELAB values
diffMatrix = compareColorFunction_LAB2(imgTestName, LAB_array); % Modif of Sam's code to accomodate using the CIELAB measurements
save('dE_img_LAB.mat','diffMatrix');

figure;
imagesc(diffMatrix);
colorbar;
title('dE (Truth CIELAB based on transmittance values)')
print('dE_heatmap_img_LAB','-dpng');

diffMatrix_rs = reshape(diffMatrix, size(diffMatrix, 1) * size(diffMatrix, 2), 1);
figure; 
boxplot(diffMatrix_rs);
title(sprintf('dE (Truth CIELAB based on transmittance values), Median = %.2f', median(diffMatrix_rs)))
print('boxplot_img_LAB','-dpng');

%% Principal Component Analysis
% Following https://stats.stackexchange.com/questions/134282/relationship-between-svd-and-pca-how-to-use-svd-to-perform-pca
b_s = 2; % First wavelength is 390nm = first band, no signal at 380 nm
data = trans_array_m(b_s:end-1, :)';
b_f = b_s + size(data, 2) - 1; % last position spectrum is 770nm = last band, no signal at 780 nm
n = size(data,1);

data_mean = mean(data,2);   % find mean
A = data - data_mean;       % subtract mean from each point, centralize data
rho = norm(A,'fro')^2;      % total variation of data, norm is sqrt(trace(A_t A)) with A centered

[U,S,V] = svd(A,'econ');    % find singular value decomposition
sigma = diag(S);            % singular values

% Reconstruct original data, limiting the number of eigenvectors
npc = 5; % Number of components
q = norm(sigma(1:npc))^2/rho;
A_pc = U(:, 1:npc)*S(1:npc, 1:npc)*V(:, 1:npc)'; 

%% Plot transmittance spectra for a series of pixels
data_pc = A_pc + data_mean;
lambda = 380:10:780;
figure;
step = 10000;
for i = 1:step:size(data_pc, 1)
    plot(lambda(b_s: b_f), data_pc(i, :), 'r'); hold on; 
    plot(lambda(b_s: b_f), data(i, :), 'b');
end
title(sprintf('Transmittance, %u components, captures %.4g%% of total variation (%% var)',npc, 100*q));
print('T_after_PCA','-dpng');

%% Plot transmittance spectrum for a one of pixel
i = 11000;
figure;
plot(lambda(b_s: b_f),data_pc(i, :), 'r'); hold on;
plot(lambda(b_s: b_f),data(i, :), 'b');
title(sprintf('T (pixel %u), %u components, captures %.4g%% of total var (%% var)', i, npc, 100*q));
print('T_one_pix_after_PCA','-dpng');

%% Mean data vector projected on eigenvector base
spectral_mean = mean(data, 1);

b_pc = spectral_mean*V(:, 1:npc); % Coefficients of the linear projection on the reduce eigenvectors base
figure; 
plot(lambda(b_s: b_f), spectral_mean, 'b'); hold on; 
plot(lambda(b_s: b_f), mean(spectral_mean) + b_pc*V(:, 1:npc)', 'r');
title(sprintf('Spectral mean, %u components, captures %.4g%% of total var (%%var)',npc, 100*q));
print('Spectral_mean_after_PCA','-dpng');

% Plot the 5 first eigenvectors with the mean spectra
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

% Trans -> XYZ -> LAB, total data set and reduced one through PCA 
sizey = size(Truth, 1);
sizex = size(Truth, 2);
[LAB_array, ~, ~, ~] = f_transmittance2LAB(trans_array_m, [], sizey, sizex, ls, 'y'); % 'y' top trim the max tranmsittance to 1
[LAB_array_5c, ~, ~, ~] = f_transmittance2LAB_pca(data_pc', [], sizey, sizex, ls, b_s, b_f, 'y'); % 'y' top trim the max tranmsittance to 1
    
% Based on Sam's code, I did it fast and dirty (no function) to compare LAB
% values of the whole data set to the one corresponding to the reconstructed data set
diff2D = (LAB_array - LAB_array_5c).^2;
diff1D = sqrt(sum(diff2D,2));
diffMatrix = reshape(diff1D, sizey, sizex);
save('dE_Truth_after_PCA.mat','diffMatrix');
 
figure;
imagesc(diffMatrix);
colorbar;
title('dE for Truth after PCA')
print('dE_Truth_heatmap_after_PCA','-dpng');

diffMatrix_rs = reshape(diffMatrix, sizex * sizey, 1);
figure; 
boxplot(diffMatrix_rs);
title(sprintf('dE for Truth after PCA, Median = %.2f', median(diffMatrix_rs)))
print('boxplot_Truth_after_PCA','-dpng');

%% Heatmap and boxplot using the HIMS CIELAB values, PCA-reconstructed
diffMatrix = compareColorFunction_LAB2(imgTestName, LAB_array_5c); % Modif of Sam's code to accomodate using the CIELAB measurements
save('dE_img_LAB_PCA_Red.mat','diffMatrix');

figure;
imagesc(diffMatrix);
colorbar;
title('dE (Truth CIELAB based on T values, PCA reduced)')
print('dE_heatmap_PCA_Red','-dpng');

diffMatrix_rs = reshape(diffMatrix, size(diffMatrix, 1) * size(diffMatrix, 2), 1);
figure; 
boxplot(diffMatrix_rs);
title(sprintf('dE (Truth CIELAB based on T values, PCA reduced), Median = %.2f', median(diffMatrix_rs)))
print('boxplot_PCA_Red','-dpng');

%% Heatmap for pilling up the 5 first eigenvectors, with data_mean + data
for i = 1:5
    A_pc = U(:, 1:i)*S(1:i, 1:i)*V(:, 1:i)';
    data_pc = A_pc + data_mean;
    [LAB_array_nc, ~, ~, ~] = f_transmittance2LAB_pca(data_pc', [], sizey, sizex, ls, b_s, b_f, 'y'); % 'y' top trim the max tranmsittance to 1
    diffMatrix(:, :, i) = compareColorFunction_LAB2(imgTestName, LAB_array_nc);
end

figure('units','normalized','outerposition',[0 0 1 1]);
for i = 1:5
    subplot(2, 3, i);
    imagesc(diffMatrix(:, :, i));
    colorbar;
    if i == 1
        title('Eigenvector 1');
    else
        title(sprintf('Eigenvector 1-%u', i));
    end
end
sgtitle('dE of Truth vs WSI');
print('dE_heatmap_eigenvector_1_5','-dpng');

figure('units','normalized','outerposition',[0 0 1 1]);
for i = 1:5
    subplot(2, 3, i);
    diffMatrix_rs = reshape(diffMatrix(:, :, i), sizex * sizey, 1);
    boxplot(diffMatrix_rs);
    if i == 1
        title(sprintf('Eigenvector 1, Median = %.2f', median(diffMatrix_rs)));
    else
        title(sprintf('Eigenvector 1-%u, Median = %.2f', i, median(diffMatrix_rs)));
    end
end
print('boxplot_eigenvector_1_5','-dpng');
 
%% Truth image using the 5 first eigenvectors and the mean vector, no trim
% for CIE coordinates computation
% XYZ, data mean    
[~, ~, XYZ_array_each_nc(:, :, 1), ~] = f_transmittance2LAB_pca(data_mean', [], sizey, sizex, ls, b_s, b_f, 'n'); % 'y' top trim the max tranmsittance to 1

% XYZ, LAB for the first 5 eigenvectors
for i = 1:5
    A_pc = U(:, i)*S(i, i)*V(:, i)';
    [~, ~, XYZ_array_each_nc(:, :, i+1), ~] = f_transmittance2LAB_pca(A_pc', [], sizey, sizex, ls, b_s, b_f, 'n'); % 'y' top trim the max tranmsittance to 1
end

% Verify by comparing the XYZ for data_nc = A_nc + data_mean to
% sum(XYZ_data_mean + XYZ_array_each_nc)
A_pc = U(:, 1:npc)*S(1:npc, 1:npc)*V(:, 1:npc)'; 
data_pc = A_pc + data_mean;
[~, ~, XYZ_array_5nc, ~] = f_transmittance2LAB_pca(data_pc', [], sizey, sizex, ls, b_s, b_f, 'n'); % 'y' top trim the max tranmsittance to 1
comp_XYZ = sum(XYZ_array_each_nc, 3) -  XYZ_array_5nc; % Small if transmittance is not trimmed in CIE computations, otherwise can be as big as 20
max(max(comp_XYZ))

% Convert XYZ to sRGB image
% Rescale XYZ so that Y of D65 illuminant is 1
Y0 = 100;

figure('units','normalized','outerposition',[0 0 1 1]);
for i = 1:6
    subplot(2, 3, i);
    % Convert to sRGB
    rgb = f_XYZ2sRGB(XYZ_array_each_nc(:, :, i)/Y0);

    % Tiff
    im = reshape(rgb,sizey,sizex,3);

    % Visualize
    image(im);
    axis image;
    
    if i == 1
        title('Data mean');
    else
        title(sprintf('Eigenvector %u', i-1));
    end
      
end
print('Truth_PCA_1_5','-dpng');

%% Reconstruct the truth image by summation individual XYZ array, convert to rgb linear then gamma correction
% Linearize
XYZ_sum = sum(XYZ_array_each_nc, 3);
rgb_l = f_rgb_linear(XYZ_sum/Y0);

% Summation
rgb_sum = rgb_l;
    
% Gamma correction
rgb = f_gamma(rgb_sum);

% Tiff
im = reshape(rgb,sizey,sizex,3);

% Visualize
figure
image(im);
axis image
title(sprintf('Truth reconstructed using %u eigenvectors', npc));
print('Truth_PCA_Reconstructed','-dpng');
