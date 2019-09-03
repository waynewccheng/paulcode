function f_lab_array_tbl(filter_list, f_name)

    %% 1: Gather results and store
    % Add 'Filter_' to filter_list
    filter_list = cellfun(@(a,b)[a b], repmat({'Filter_'}, 1, size(filter_list, 2)), filter_list, 'Uniform', 0);
    
    % Storage
    n_filters = size(filter_list, 2);
    lab_array_tbl = zeros(676*844, 3, n_filters);
    cov_lab_array_tbl = zeros(3, 3, 676*844, n_filters);
        
    % Get lab_array values
    for i = 1:n_filters
        
        % Which folder?
        filter_name = filter_list{i};
        
        % Load CIE data
        load(['output\' filter_name '\LAB_array'],'LAB_array');
        load(['output\' filter_name '\CovLAB_array'],'CovLAB_array');
        lab_array_tbl(:, :, i) = LAB_array;
        cov_lab_array_tbl(:, :, :, i) = CovLAB_array;
    end
    
    %% 2: Save the output
    fld_name = ['output\' f_name];
    mkdir(fld_name);

    % CIELAB results
    csvwrite([fld_name '\lab_array_tbl.txt'], lab_array_tbl);
    save([fld_name '\lab_array_tbl'], 'lab_array_tbl');
    
    csvwrite([fld_name '\cov_lab_array_tbl.txt'], cov_lab_array_tbl);
    save([fld_name '\cov_lab_array_tbl'], 'cov_lab_array_tbl');

end

