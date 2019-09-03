%% For verification of the data
% rawdata (measurements) are in input folder
% processed data are in output folder

% ouputs the transmittance spectra, CIEXYZ, CIELAB and DeltaE values

% 08-29-19: First code

function t_dE_KWColFiltwcc(filter_id)

    close all;
    
    %% 1: Initialization
    % Folders for camera measurements
    filter_name = ['Filter_' filter_id];
    foldername_sample = [filter_name '_sample'];  % For the filter spectra
    foldername_white = [filter_name '_white'];    % For the 100% tranmittance

    % Folder for spectro-radiometer measurements
    foldername_spectro = [filter_name '_spectro_8xFast'];

    %% 2: Spectrometer: computes transmittance
    % Number of repeated measurements
    n_meas = 10;

    % Spectrometer spectrum range
    lambda_s = 380:780;

    % Load the spectrometer data
    load(['input\' filter_name '\' foldername_spectro '\spectro_meas'],'spectra')
    load(['input\' filter_name '\' foldername_spectro '\spectro_meas_background'],'spectra_background')

    % Compute mean value and error with OL490 on
    s_filter = spectra(:, :, 1);
    s_white = spectra(:, :, 2);

    s_filter_m = mean(s_filter);
    s_white_m = mean(s_white);
    s_filter_s = std(s_filter)./sqrt(n_meas);
    s_white_s = std(s_white)./sqrt(n_meas);

    % Compute mean value and error with OL490 off (background)
    s_filter_background = spectra_background(:, :, 1);
    s_white_background = spectra_background(:, :, 2);

    s_filter_background_m = mean(s_filter_background);
    s_white_background_m = mean(s_white_background);
    s_filter_background_s = std(s_filter_background)./sqrt(n_meas);
    s_white_background_s = std(s_white_background)./sqrt(n_meas);

    % Compute the transmittance
    [t_mean_spectro, t_std_spectro] = f_transmittance_PL3(s_filter_m, s_white_m,...
        s_filter_background_m, s_white_background_m, s_filter_s, s_white_s,...
        s_filter_background_s, s_white_background_s);

    trans_spectro = [lambda_s; t_mean_spectro; t_std_spectro]';
    %% 3: Spectrometer: computes LAB (T -> XYZ -> LAB)
    % Prepares the illuminant
    load ('input\DataIlluminants\spec_cied65','spec');
    ls = spec(1:10:401,2);

    % Compute LAB for the spectro
    [LAB_spectro, CovLAB_spectro, XYZ_spectro, CovXYZ_spectro] = transmittance2LAB_PL2(t_mean_spectro(1:10:401)', t_std_spectro(1:10:401)', 1, 41, ls);

    %% 6: Camera: computes transmittance
    % Compute the tranmittance based on the spatial average of numberofshots
    % images and the corresponding stat stored in img_ms
    % , trans_ms is the spatial average with temporal mean and std
    % dev of the transmittance, trans_array_m and trans_array_s are pixel by
    % pixel values
    numberofshots = 10;
    f_sample = ['input\' filter_name '\' foldername_sample];
    f_white = ['input\' filter_name '\' foldername_white];

    [trans_cam_ms, trans_array_m, trans_array_s, sizey, sizex] = frame2transmittance_white_PL4(f_sample, f_white, numberofshots);

    %% 7: Camera: computes LAB
    % Trans -> XYZ -> LAB with saving XYZ and LAB (values + uncertainties)
    [LAB_cam, CovLAB_cam, XYZ_cam, CovXYZ_cam] = transmittance2LAB_PL2(trans_cam_ms(:, 2), trans_cam_ms(:, 3), 1, 41, ls);
    [LAB_array, CovLAB_array, XYZ_array, CovXYZ_array] = transmittance2LAB_PL2(trans_array_m, trans_array_s, sizey, sizex, ls);

    %% 8: Graphics
    figure;
    errorbar(trans_spectro(1:10:end, 1), trans_spectro(1:10:end, 2), 2 * trans_spectro(1:10:end, 3)); hold on;
    errorbar(trans_cam_ms(:, 1), trans_cam_ms(:, 2) , 2 * trans_cam_ms(:, 3) , 'k');
    axis([350 800 -0.1 1]);
    legend('Spectro', 'Whole img');
    title('Error bars at k = 2');

    figure;
    errorbar(trans_spectro(1:10:end, 1), trans_spectro(1:10:end, 2), 2 * trans_spectro(1:10:end, 3)); hold on;
    errorbar(trans_spectro(1:10:end, 1), trans_array_m(:, 1*1) , 2 * trans_array_s(:, 1*1) , 'k');
    legend('Spectro', 'one pixel');
    title('Error bars at k = 2');

    % CIELAB space
    figure;
    step = 500;
    [DE(1, 1), DE(1, 2)] = f_deltaE(LAB_spectro, LAB_cam, CovLAB_spectro, CovLAB_cam);
    scatter3(LAB_array(1:step:end, 3), LAB_array(1:step:end, 2), LAB_array(1:step:end, 1), '.b'); hold on;
    scatter3(LAB_cam(3), LAB_cam(2), LAB_cam(1), 'r', 'Filled');
    scatter3(LAB_spectro(3), LAB_spectro(2), LAB_spectro(1), 'k', 'LineWidth', 2);
    xlabel('b^*'); ylabel('a^*'); zlabel('L^*');
    legend('Pixel', 'Img mean', 'Spectro');
    title(['\Delta E^*_{ab} = ' num2str(DE(1, 1)) ' \sigma_{\Delta E^*_{ab}} = ' num2str(DE(1, 2))]);

    % Beam profile
    fnin_m = sprintf('%s/vim_mean_array',f_white);
    load(fnin_m,'vim_mean_array');

    vt = vim_mean_array;
    figure('units','normalized','outerposition',[0 0 1 1]);
    sgtitle('Beam profile');
    wl_array = 380:10:780;
    for wl = 1:41
        subplot(6,7,wl)

        vvname = sprintf('%s(wl,:,:)','vt');
        vv = eval(vvname);
        im = squeeze(vv);
        imagesc(im)
        axis off
        axis image
        colorbar
        title(sprintf('%d',wl_array(wl)))
    end


    % Tranmittance images
    vt = reshape(trans_array_m, size(trans_array_m, 1), sizey, sizex);
    figure('units','normalized','outerposition',[0 0 1 1]);
    sgtitle('Transmittance images');
    wl_array = 380:10:780;
    for wl = 1:41
        subplot(6,7,wl)

        vvname = sprintf('%s(wl,:,:)','vt');
        vv = eval(vvname);
        im = squeeze(vv);
        imagesc(im)
        axis off
        axis image
        colorbar
        title(sprintf('%d',wl_array(wl)))
    end

    %% 9: Save the output
    % Save the data
    fld_name = ['output\' filter_name];
    mkdir(fld_name);

    save([fld_name '\trans_spectro'],'trans_spectro');
    save([fld_name '\trans_cam_ms'],'trans_cam_ms');
    save([fld_name '\trans_mean_camera'],'trans_array_m','sizey','sizex','-v7.3');
    save([fld_name '\trans_std_camera'],'trans_array_s','sizey','sizex','-v7.3');

    % CIELAB coordinates
    save([fld_name '\LAB_cam'],'LAB_cam');
    save([fld_name '\CovLAB_cam'],'CovLAB_cam');
    save([fld_name '\XYZ_cam'],'XYZ_cam');
    save([fld_name '\CovXYZ_cam'],'CovXYZ_cam');

    save([fld_name '\LAB_spectro'],'LAB_spectro');
    save([fld_name '\CovLAB_spectro'],'CovLAB_spectro');
    save([fld_name '\XYZ_spectro'],'XYZ_spectro');
    save([fld_name '\CovXYZ_spectro'],'CovXYZ_spectro');

    save([fld_name '\LAB_array'],'LAB_array');
    save([fld_name '\CovLAB_array'],'CovLAB_array');
    save([fld_name '\XYZ_array'],'XYZ_array');
    save([fld_name '\CovXYZ_array'],'CovXYZ_array');

    save([fld_name '\DeltaE'],'DE');
    
end