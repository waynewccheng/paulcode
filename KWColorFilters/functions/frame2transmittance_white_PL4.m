% 08-20-19
% No black tape measurements

% 08-19-19
% Measurements of the light background taken into account

% 7-30-2015
% convert frames (DDL) to transmittance by using reference white background
%

% function [transmittance_array, sizey, sizex] = frame2transmittance_white_PL (foldername_sample, foldername_white, nshots)
function [trans_ms, trans_array_m, trans_array_s, sizey, sizex] = frame2transmittance_white_PL4(foldername_sample, foldername_white, nshots)

    % Compute the spatial + temporal mean and std dev
    
    disp('Combining spatial + temporal mean and std dev into transmittance...')
    
    % 100% transmittance
    fnin = sprintf('%s/img_ms',foldername_white);
    load(fnin,'img_ms');
    fnin = sprintf('%s/img_background_ms',foldername_white);
    load(fnin,'img_background_ms');
    
    % Copy
    img_ms_w = img_ms;
    img_background_ms_w = img_background_ms;
    
    % Sample transmittance
    fnin = sprintf('%s/img_ms',foldername_sample);
    load(fnin,'img_ms');
    fnin = sprintf('%s/img_background_ms',foldername_sample);
    load(fnin,'img_background_ms');
    
    % Transmittance
    trans_ms(:, 1) = img_ms(:, 1);
    [trans_ms(:, 2), trans_ms(:, 3)] = f_transmittance_PL3(img_ms(:, 2), img_ms_w(:, 2), ...
        img_background_ms(:, 2), img_background_ms_w(:, 2),...
        img_ms(:, 3), img_ms_w(:, 3),...
        img_background_ms(:, 3), img_background_ms_w(:, 3));

    disp('Combining frames into transmittance...')

    % Load images mean values for the 100% transmittance
    fnin_m = sprintf('%s/vim_mean_array',foldername_white);
    load(fnin_m,'vim_mean_array');
    fnin_s = sprintf('%s/vim_std_array',foldername_white);
    load(fnin_s,'vim_std_array');

    fnin_m = sprintf('%s/vim_background_mean_array',foldername_white);
    load(fnin_m,'vim_background_mean_array');
    fnin_s = sprintf('%s/vim_background_std_array',foldername_white);
    load(fnin_s,'vim_background_std_array');
    
    % Copy to vim_mean_array_w, to save it and get dimensions
    vim_mean_array_w = vim_mean_array;
    vim_std_array_w = vim_std_array;
    vim_background_mean_array_w = vim_background_mean_array;
    vim_background_std_array_w = vim_background_std_array;
    [sizewl sizey sizex] = size(vim_mean_array_w);

    % Load images mean values for the sample
    fnin_m = sprintf('%s/vim_mean_array',foldername_sample);
    load(fnin_m,'vim_mean_array');
    fnin_s = sprintf('%s/vim_std_array',foldername_sample);
    load(fnin_s,'vim_std_array');

    fnin_m = sprintf('%s/vim_background_mean_array',foldername_sample);
    load(fnin_m,'vim_background_mean_array');
    fnin_s = sprintf('%s/vim_background_std_array',foldername_sample);
    load(fnin_s,'vim_background_std_array');

    % calculate the reflectance
    ddl_array_m = reshape(vim_mean_array, sizewl, sizey*sizex);
    ddl_white_array_m = reshape(vim_mean_array_w, sizewl, sizey*sizex);

    ddl_background_array_m = reshape(vim_background_mean_array, sizewl, sizey*sizex);
    ddl_background_white_array_m = reshape(vim_background_mean_array_w, sizewl, sizey*sizex);
    
    ddl_array_s = reshape(vim_std_array, sizewl, sizey*sizex);
    ddl_white_array_s = reshape(vim_std_array_w, sizewl,sizey*sizex);

    ddl_background_array_s = reshape(vim_background_std_array, sizewl, sizey*sizex);
    ddl_background_white_array_s = reshape(vim_background_std_array_w, sizewl,sizey*sizex);
    
    [trans_array_m , trans_array_s] = f_transmittance_PL3(ddl_array_m, ddl_white_array_m,...
        ddl_background_array_m, ddl_background_white_array_m,...
        ddl_array_s, ddl_white_array_s,...
        ddl_background_array_s, ddl_background_white_array_s);

    % ----------------------------------
    return
end
