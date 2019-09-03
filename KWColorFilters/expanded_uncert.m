%% Call f_expanded_uncert
filter_list = {'KW12BW10', 'KW25BW10', 'KW32eBW10', 'KW47eBW10', 'KW58BW10'};
filter_list_repro = {'KW32BW10', 'KW32BW10', 'KW32BW10', 'KW47BW10', 'KW32BW10'};
fld_name = 'Results_KW_Filters';
    
f_expanded_uncert(filter_list_repro, filter_list, fld_name);

%% Call f_lab_array_tbl to get the images LAB values
f_lab_array_tbl(filter_list, fld_name);
