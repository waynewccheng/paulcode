function [t_spectro_tbl_m, t_spectro_tbl_s, t_cam_tbl_m, t_cam_tbl_s, lab_spectro_tbl, lab_cam_tbl, DE_tbl, n_filters] = f_gather_results(filter_list)

    % Storage
    n_filters = size(filter_list, 2);
    t_spectro_tbl_m = zeros(401, n_filters+1);
    t_spectro_tbl_s = zeros(401, n_filters+1);
    t_cam_tbl_m = zeros(41, n_filters+1);
    t_cam_tbl_s = zeros(41, n_filters+1);
    lab_spectro_tbl = zeros(n_filters, 6);
    lab_cam_tbl = zeros(n_filters, 6);
    DE_tbl = zeros(n_filters, 2);
    
    for i = 1:n_filters

        % Which folder?
        filter_name = filter_list{i};

        % Load the spectro data
        load(['output\' filter_name '\trans_spectro'],'trans_spectro');
        t_spectro_tbl_m(:, i+1) = trans_spectro(:, 2);
        t_spectro_tbl_s(:, i+1) = trans_spectro(:, 3);

        % Load the camera data
        load(['output\' filter_name '\trans_cam_ms'],'trans_cam_ms');
        t_cam_tbl_m(:, i+1) = trans_cam_ms(:, 2);
        t_cam_tbl_s(:, i+1) = trans_cam_ms(:, 3);

        % Load CIE data
        load(['output\' filter_name '\LAB_Spectro'],'LAB_spectro');
        load(['output\' filter_name '\CovLAB_Spectro'],'CovLAB_spectro');
        load(['output\' filter_name '\LAB_cam'],'LAB_cam');
        load(['output\' filter_name '\CovLAB_cam'],'CovLAB_cam');

        lab_spectro_tbl(i, 1:3) = LAB_spectro;
        lab_spectro_tbl(i, 4:6) = sqrt(diag(CovLAB_spectro))';

        lab_cam_tbl(i, 1:3) = LAB_cam;
        lab_cam_tbl(i, 4:6) = sqrt(diag(CovLAB_cam))';

        [DE_tbl(i, 1), DE_tbl(i, 2)] = f_deltaE(LAB_spectro, LAB_cam, CovLAB_spectro, CovLAB_cam);

    end

    t_spectro_tbl_m(:, 1) = trans_spectro(:, 1);
    t_cam_tbl_m(:, 1) = trans_cam_ms(:, 1);
    t_spectro_tbl_s(:, 1) = trans_spectro(:, 1);
    t_cam_tbl_s(:, 1) = trans_cam_ms(:, 1);
    
end

