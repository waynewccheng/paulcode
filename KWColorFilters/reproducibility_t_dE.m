%% Outputs the uncertainty results due to the reproducibility experiments
% in xxx_repro

% Ouputs tables of the transmittance, LAB coordinates and DeltaE over the
% reproduced experiments in xxx_tbl_x

% data are in output folder
% results are in the output folder

% 08-30-19: f_gather_results function

% 08-29-19: First code

function reproducibility_t_dE(filter_id)

    %% 1: Gather the results  
    tmp = struct2cell(dir(['output\Filter_' filter_id '*'])); % Transform to cell array to accomodate f_gather_results function
    filter_list = tmp(1, :);

    [t_spectro_tbl_m, t_spectro_tbl_s, t_cam_tbl_m, t_cam_tbl_s, lab_spectro_tbl, lab_cam_tbl, DE_tbl, nb_exp] = f_gather_results(filter_list);
    
    %% 2: Standard deviation over the experiments
    t_spectro_repro(:, 1) = t_spectro_tbl_m(:, 1);
    t_cam_repro(:, 1) = t_cam_tbl_m(:, 1);

    t_spectro_repro(:, 2) = std(t_spectro_tbl_m(:, 2:nb_exp+1), [], 2);
    t_cam_repro(:, 2) = std(t_cam_tbl_m(:, 2:nb_exp+1), [], 2);
    DE_repro(1, 1) = std(DE_tbl(1:nb_exp, 1));
    lab_spectro_repro(1, 1:3) = std(lab_spectro_tbl(:, 1:3),1);
    lab_cam_repro(1, 1:3) = std(lab_cam_tbl(:, 1:3),1);

    %% 3: Save the output
    fld_name = ['output\Repro_' filter_list{1}];
    mkdir(fld_name);

    save([fld_name '\t_spectro_tbl_m'],'t_spectro_tbl_m');
    save([fld_name '\t_spectro_tbl_s'],'t_spectro_tbl_s');
    save([fld_name '\t_cam_tbl_m'],'t_cam_tbl_m');
    save([fld_name '\t_cam_tbl_s'],'t_cam_tbl_s');

    save([fld_name '\lab_spectro_tbl'],'lab_spectro_tbl');
    save([fld_name '\lab_cam_tbl'],'lab_cam_tbl');
    save([fld_name '\DE_tbl'],'DE_tbl');

    save([fld_name '\t_spectro_repro'],'t_spectro_repro');
    save([fld_name '\t_cam_repro'],'t_cam_repro');

    save([fld_name '\lab_spectro_repro'],'lab_spectro_repro');
    save([fld_name '\lab_cam_repro'],'lab_cam_repro');
    save([fld_name '\DE_repro'],'DE_repro')

end
