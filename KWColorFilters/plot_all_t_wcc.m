% 09-03-19: add lab_array_tbl

% 08-30-19: first code

close all;
clearvars;
%% 1: Load filter results, includes Perkin Elmer Lambda 1050 measurements
% Filter names
fld_name = 'output\Results_KW_Filters';
load([fld_name '\Names_KW_ColFilters'],'filter_list');
n_filters = size(filter_list, 2);

% Transmittance Spectro
load([fld_name '\t_spectro_cmp'],'t_spectro_cmp');

% Transmittance Camera
load([fld_name '\t_cam_cmp'],'t_cam_cmp');

% LAB Spectro
load([fld_name '\lab_spectro_cmp'],'lab_spectro_cmp');

% LAB Camera
load([fld_name '\lab_cam_cmp'],'lab_cam_cmp');
load([fld_name '\lab_array_tbl'],'lab_array_tbl');
load([fld_name '\cov_lab_array_tbl'],'cov_lab_array_tbl');

% Delta E
load([fld_name '\DE_cmp'],'DE_cmp');

% Perkin Elmer measurements
filter_list_sp1050 = {'kw12.Sample.Raw.csv', 'kw25.Sample.Raw.csv', 'kw32.Sample.Raw.csv', 'kw47.Sample.Raw.csv', 'kw58.Sample.Raw.csv'};
pe_spectro = zeros(431, 6);

for i = 1:size(filter_list_sp1050, 2)
    tmp = readmatrix(['input\KodakWrattenFilter_PEMeas\' filter_list_sp1050{i}]);
    pe_spectro(:, i+1) = tmp(:, 2);
end

pe_spectro(:, 1) = tmp(:, 1);

%% 2: Graphics

% Unceratinty choice: 'sigma' or 'expanded'
% Save curves of not: 'on' or 'off'

% All filters
plot_all_t(t_spectro_cmp, t_cam_cmp, pe_spectro, lab_spectro_cmp, 'expanded', 'on', fld_name, n_filters); % unceratinty choice: 'sigma' or 'expanded'

% KW32
plot_t_lab(t_spectro_cmp, t_cam_cmp, lab_spectro_cmp, lab_cam_cmp, lab_array_tbl, DE_cmp, 'expanded', 'on', 'Filter_KW32eBW10', filter_list, fld_name);

% KW47
plot_t_lab(t_spectro_cmp, t_cam_cmp, lab_spectro_cmp, lab_cam_cmp, lab_array_tbl, DE_cmp, 'expanded', 'on', 'Filter_KW47eBW10', filter_list, fld_name);

% Boxplot
statgraph_dE(lab_array_tbl, lab_spectro_cmp, 'on', filter_list, fld_name);


%% 3: Functions
function plot_all_t(t_spectro_cmp, t_cam_cmp, pe_spectro, lab_spectro_cmp, uncert_opt, save_opt, fld_name, n_filters)

    fig1 = figure('units','normalized','outerposition',[0 0 1 1]);
    fig2 = figure('units','normalized','outerposition',[0 0 1 1]);  
    
    j = 2;
    for i = 1:n_filters

        % Color
        c = double(lab2rgb(lab_spectro_cmp(i, 1:5:11),'OutputType','uint8'))/255;

        % Graphic
        figure(fig1); axis([350 800 -0.1 1]);
        plot(t_spectro_cmp(:, 1), t_spectro_cmp(:, j), '.-', 'Color', c); hold on;
        plot(pe_spectro(:, 1), pe_spectro(:, i+1)/100, '-', 'Color', c);
        switch uncert_opt
            case 'sigma'
                errorbar(t_cam_cmp(:, 1), t_cam_cmp(:, j), 2*t_cam_cmp(:, j+1), '--', 'Color', c);
            case 'expanded'
                errorbar(t_cam_cmp(:, 1), t_cam_cmp(:, j), t_cam_cmp(:, j+4), '--', 'Color', c); % k = 2 already included
        end

        xlabel('\lambda (nm)'); ylabel('T');

        figure(fig2); axis([350 800 -0.1 1.1]);
        errorbar(t_spectro_cmp(1:10:end, 1), t_spectro_cmp(1:10:end, j), 2*t_spectro_cmp(1:10:end, j+1), '.-', 'Color', c); 
        hold on;
        switch uncert_opt
            case 'sigma'
                errorbar(t_cam_cmp(:, 1), t_cam_cmp(:, j), 2*t_cam_cmp(:, j+1), '--', 'Color', c);
            case 'expanded'
                errorbar(t_cam_cmp(:, 1), t_cam_cmp(:, j), t_cam_cmp(:, j+4), '--', 'Color', c); % k = 2 already included
        end
        xlabel('\lambda (nm)'); ylabel('T');
        
        j = j + 5;
    end
    
    % Save figures in tif format
    switch save_opt
        case 'on'
            saveas(fig1,[fld_name '\T_KW_Filters' uncert_opt '.tif']);
            saveas(fig2,[fld_name '\T_KW_Filters_NoPE_L1050' uncert_opt '.tif']);
        case 'off'
    end
    
end

function plot_t_lab(t_spectro_cmp, t_cam_cmp, lab_spectro_cmp, lab_cam_cmp, lab_array_tbl, DE_cmp, uncert_opt, save_opt, filter_name, filter_list, fld_name)

    fig1 = figure; 
    fig2 = figure;
    i = find(ismember(filter_list,filter_name));
    j = 2 + 5*(i-1);
    
    % Color
    c = double(lab2rgb(lab_spectro_cmp(i, 1:5:11),'OutputType','uint8'))/255;

    % Graphics
    figure(fig1);
    errorbar(t_spectro_cmp(1:10:end, 1), t_spectro_cmp(1:10:end, j), 2*t_spectro_cmp(1:10:end, j+1), '.-', 'Color', 'b'); hold on;
    
    switch uncert_opt
        case 'sigma'
            errorbar(t_cam_cmp(:, 1), t_cam_cmp(:, j), 2*t_cam_cmp(:, j+1), '--', 'Color', 'r');
        case 'expanded'
            errorbar(t_cam_cmp(:, 1), t_cam_cmp(:, j), t_cam_cmp(:, j+4), '--', 'Color', 'r'); % k = 2 already included
    end
    xlabel('\lambda (nm)'); ylabel('T');
%     legend('T_S', 'T_{SA}', 'Location','northwest');
    legend('T_{Spectro}', 'T_{Cam}', 'Location','northwest');
    axis([350 800 -0.1 1.1]);
    
    figure(fig2);
    step = 500;
    scatter3(lab_array_tbl(1:step:end, (i-1)*3+3), lab_array_tbl(1:step:end, (i-1)*3+2), lab_array_tbl(1:step:end, (i-1)*3+1),...
        'MarkerEdgeColor', c, 'MarkerEdgeAlpha', 0.25); hold on;
    scatter3(lab_cam_cmp(i, 11), lab_cam_cmp(i, 6), lab_cam_cmp(i, 1), 'r', 'Filled', 'LineWidth', 0.6);
    scatter3(lab_spectro_cmp(i, 11), lab_spectro_cmp(i, 6), lab_spectro_cmp(i, 1), 'k', 'Filled', 'LineWidth', 0.6);
    xlabel('b^*'); ylabel('a^*'); zlabel('L^*');
    legend('Pixels', 'Img Spatial Avg', 'Spectroradiometer', 'Location','northwest');
    title(['\Delta E_{ab}^* = ' sprintf('%0.2f',DE_cmp(i, 1)) '; U_{\Delta E_{ab}^*} = ' sprintf('%0.2f', DE_cmp(i, 5))]);

    % Save figures in tif format
    switch save_opt
        case 'on'
            saveas(fig1,[fld_name '\T_' filter_name ' ' uncert_opt '.tif']);
            saveas(fig2,[fld_name '\LAB_' filter_name ' ' uncert_opt '.tif']);
        case 'off'
    end

end

function statgraph_dE(lab_array_tbl, lab_spectro_cmp,  save_opt, filter_list, fld_name)

    % Delta E, a bit different than in f_deltaE, can treat lists
    DeltaE = @(LAB_1, LAB_2) sqrt((LAB_1(:, 1) - LAB_2(:, 1)).^2 + ...    % L
    (LAB_1(:, 2) - LAB_2(:, 2)).^2 + ...                                  % a
    (LAB_1(:, 3) - LAB_2(:, 3)).^2);                                      % b

    % Storage
    n_filters = size(filter_list, 2);
    dE_array = zeros(676*844, n_filters);
    
    % Delta E
    tmp = lab_spectro_cmp(:, 1:5:11);
    lab_spectro_array = permute(repmat(tmp', 1, 1, size(lab_array_tbl, 1)), [3 1 2]);
    
    for i = 1:n_filters
        dE_array(:, i) = DeltaE(lab_array_tbl(:, :, i), lab_spectro_array(:, :, i));
    end

    % Graphics
    fig1 = figure('units','normalized','outerposition',[0 0 1 1]);
    figure(fig1); boxplot(dE_array, filter_list);
    
    fig2 = figure('units','normalized','outerposition',[0 0 1 1]);
    figure(fig2); 
    
    for i = 1:n_filters
        subplot(2, 3, i);
        histogram(dE_array(:, i));
    end
    
    % Save figures in tif format
    switch save_opt
        case 'on'
            saveas(fig1,[fld_name '\KW_Filters_DeltaE_Bxplt.tif']);
            saveas(fig2,[fld_name '\KW_Filters_DeltaE_Histo.tif']);
        case 'off'
    end
    
end

