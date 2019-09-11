%% Outputs the compound tables of the transmittance, CIELAB coordinates and DeltaE values
% Tables with values, std dev, reproducibility uncertainties, total type A
% and expanded uncertaintites

% 08-30-19: First code

function f_expanded_uncert(filter_list_repro, filter_list, f_name)

    %% 1: Gather results
    % Add 'Filter_' to filter_list to accomodate for f_gather_results
    % function
    filter_list = cellfun(@(a,b)[a b], repmat({'Filter_'}, 1, size(filter_list, 2)), filter_list, 'Uniform', 0);
    [t_spectro_tbl_m, t_spectro_tbl_s, t_cam_tbl_m, t_cam_tbl_s, lab_spectro_tbl, lab_cam_tbl, DE_tbl, n_filters] = f_gather_results(filter_list);
    
    % Get reproducibility uncertainties
    for i = 1:n_filters
        
        % Load
        fld_name = ['output\Repro_Filter_' filter_list_repro{i}];
        load([fld_name '\t_spectro_repro'], 't_spectro_repro');
        load([fld_name '\t_cam_repro'], 't_cam_repro');
        load([fld_name '\lab_spectro_repro'], 'lab_spectro_repro');
        load([fld_name '\lab_cam_repro'], 'lab_cam_repro');
        load([fld_name '\DE_repro'], 'DE_repro');
    
        % Store
        t_spectro_repro_tbl(:, i+1) = t_spectro_repro(:, 2);
        t_cam_repro_tbl(:, i+1) = t_cam_repro(:, 2);
        lab_spectro_repro_tbl(i, :) = lab_spectro_repro;
        lab_cam_repro_tbl(i, :) = lab_cam_repro;
        DE_repro_tbl(i, :) = DE_repro;
        
    end
    
    %% 2: Create the compound tables with the extended uncertainties
    % Spectro: Transmittance table with uncertaintities, Lambda, value, std dev,
    % repro uncertainty, total Type A, Expanded (k=2)
    t_spectro_cmp = zeros(size(t_spectro_tbl_m, 1), 5*n_filters+1);
    t_spectro_cmp(:, 1) = t_spectro_tbl_m(:, 1);
    t_spectro_cmp(:, 2:5:5*n_filters-3) = t_spectro_tbl_m(:, 2:end);
    t_spectro_cmp(:, 3:5:5*n_filters-2) = t_spectro_tbl_s(:, 2:end);
    t_spectro_cmp(:, 4:5:5*n_filters-1) = t_spectro_repro_tbl(:, 2:end);
    t_spectro_cmp(:, 5:5:5*n_filters) = sqrt(t_spectro_tbl_s(:, 2:end).^2+t_spectro_repro_tbl(:, 2:end).^2);
    t_spectro_cmp(:, 6:5:5*n_filters+1) = 2*t_spectro_cmp(:, 5:5:5*n_filters);

    % Camera: Transmittance table with uncertaintities, Lambda, value, std dev,
    % repro uncertainty, total Type A, Expanded (k=2)
    t_cam_cmp = zeros(size(t_cam_tbl_m, 1), 5*n_filters+1);
    t_cam_cmp(:, 1) = t_cam_tbl_m(:, 1);
    t_cam_cmp(:, 2:5:5*n_filters-3) = t_cam_tbl_m(:, 2:end);
    t_cam_cmp(:, 3:5:5*n_filters-2) = t_cam_tbl_s(:, 2:end);
    t_cam_cmp(:, 4:5:5*n_filters-1) = t_cam_repro_tbl(:, 2:end);
    t_cam_cmp(:, 5:5:5*n_filters) = sqrt(t_cam_tbl_s(:, 2:end).^2+t_cam_repro_tbl(:, 2:end).^2);
    t_cam_cmp(:, 6:5:5*n_filters+1) = 2*t_cam_cmp(:, 5:5:5*n_filters);
    
    % Spectro: LAB coordinates value, std dev, repro uncertainty, 
    % total Type A, Expanded (k=2)
    lab_spectro_cmp(:, 1:5:15-4) = lab_spectro_tbl(:, 1:3);
    lab_spectro_cmp(:, 2:5:15-3) = lab_spectro_tbl(:, 4:6);
    lab_spectro_cmp(:, 3:5:15-2) = lab_spectro_repro_tbl(:, 1:3);
    lab_spectro_cmp(:, 4:5:15-1) = sqrt(lab_spectro_tbl(:, 4:6).^2+lab_spectro_repro_tbl(:, 1:3).^2);
    lab_spectro_cmp(:, 5:5:15) = 2*lab_spectro_cmp(:, 4:5:15-1);
    
    % Camera: LAB coordinates value, std dev, repro uncertainty,
    % total Type A, Expanded (k=2)
    lab_cam_cmp(:, 1:5:15-4) = lab_cam_tbl(:, 1:3);
    lab_cam_cmp(:, 2:5:15-3) = lab_cam_tbl(:, 4:6);
    lab_cam_cmp(:, 3:5:15-2) = lab_cam_repro_tbl(:, 1:3);
    lab_cam_cmp(:, 4:5:15-1) = sqrt(lab_cam_tbl(:, 4:6).^2+lab_cam_repro_tbl(:, 1:3).^2);
    lab_cam_cmp(:, 5:5:15) = 2*lab_cam_cmp(:, 4:5:15-1);
    
    % Delta E: value, std dev, repro uncertainty,
    % total Type A, Expanded (k=2)
    DE_cmp(:, 1) = DE_tbl(:, 1);
    DE_cmp(:, 2) = DE_tbl(:, 2);
    DE_cmp(:, 3) = DE_repro_tbl(:, 1);
    DE_cmp(:, 4) = sqrt(DE_tbl(:, 2).^2+DE_repro_tbl(:, 1).^2);
    DE_cmp(:, 5) = 2*DE_cmp(:, 4);
    
    
    %% 3: Save the output
    fld_name = ['output\' f_name];
    mkdir(fld_name);

    % Names of filters
    name = f_name(9:end-8);
    str_tbl = string(filter_list);

    filepath = fopen([fld_name '\Names_' name '_ColFilters.txt'],'w');
    fprintf(filepath,'%s,\n',str_tbl{:});
    fclose(filepath);
    
    save([fld_name '\Names_' name '_ColFilters'],'filter_list');
    
    % Transmittances
    csvwrite([fld_name '\t_spectro_cmp.txt'], t_spectro_cmp);
    csvwrite([fld_name '\t_cam_cmp.txt'], t_cam_cmp);
    
    save([fld_name '\t_spectro_cmp'],'t_spectro_cmp');
    save([fld_name '\t_cam_cmp'],'t_cam_cmp');
    
    % CIELAB results
    csvwrite([fld_name '\DE_' name '_ColFilters.txt'], DE_cmp);
    csvwrite([fld_name '\LAB_spectro_' name '_ColFilters.txt'], lab_spectro_cmp);
    csvwrite([fld_name '\LAB_cam_' name '_ColFilters.txt'], lab_cam_cmp);
    
    save([fld_name '\DE_cmp'], 'DE_cmp');
    save([fld_name '\lab_spectro_cmp'], 'lab_spectro_cmp');
    save([fld_name '\lab_cam_cmp'], 'lab_cam_cmp');
    
end
