Kodak Wratten color filters measurements on the lab computer:
KW #12, #25, #32, #47, #58

Bandwidth of OL490: 10nm

Lab computer folders
Rawdata:
\\192.168.44.167\f\Data_Paul\RawData\082019
\\192.168.44.167\f\Data_Paul\RawData\082119
copied in input folder

ProcessedData: (recomputed here and paced in output folder)
\\192.168.44.167\f\Data_Paul\ProcessedData\082019
\\192.168.44.167\f\Data_Paul\ProcessedData\082119

PerkinElmer Lambda 1050 measurements:
in input folder

The folders are named as:
Filter_KW12BW10, Filter_KW25BW10, Filter_KW32BW10, Filter_KW47BW10, Filter_KW58BW10

Reproducibility measurements on #32 and #47 (10 experiments), the folders are named as:

Filter_KW32BW10, Filter_KW32bBW10, Filter_KW32cBW10, ..Filter_KW32jBW10
Filter_KW47BW10, Filter_KW47bBW10, Filter_KW47cBW10, ..Filter_KW47jBW10

* t_dE_KWColFiltwcc(filter_id): run for each input folder

	filter_id = 'KW12BW10' etc..

	output files in dedicated folders:

   	trans_spectro: 401 x 3 (Lambda, mean T, std dev T)
    	trans_cam_ms: 41 x 3 (Lambda, mean T, std dev T)
    	trans_array_m (trans_mean_camera): 41 X (676x844)
    	trans_array_s (trans_std_camera): 41 X (676x844)

    	LAB_cam: 1 x 3
    	CovLAB_cam: 3 x 3
    	XYZ_cam: 1 x 3
    	CovXYZ_cam: 3: 3

    	LAB_spectro: 1 x 3
    	CovLAB_spectro: 3 x 3
    	XYZ_spectro: 1 x 3
    	CovXYZ_spectro: 3 x 3

    	LAB_array: (676x844) x 3
    	CovLAB_array: 3 x 3 x (676x844)
    	XYZ_array: (676x844) x 3
    	CovXYZ_array: 3 x 3 x (676x844)

    	DE (DeltaE): 1 x 2

* reproducibility_t_dE(filter_id)

	filter_id = 'KW32' and 'KW47'

	The reproducibility results are in \output\Repro_Filter_KW32BW10 and Repro_Filter_KW47BW10:

	t_spectro_tbl_m: mean values for each sample (10 measurements each time), 10 samples
	401 rows x 11 cols, Lambda (col 1), T values (col 2:11)

	t_spectro_tbl_s: std dev values for each sample (10 measurements each time), 10 samples
	401 rows x 11 cols, Lambda (col 1), T values (col 2:11)

	t_cam_tbl_m: mean values for each sample (10 measurements each time), 10 samples
	41 rows x 11 cols, Lambda (col 1), T values (col 2:11)

	t_cam_tbl_s: std dev values for each sample (10 measurements each time), 10 samples
	41 rows x 11 cols, Lambda (col 1), T values (col 2:11)

	lab_spectro_tbl: 10 samples
	10 rows X 6 cols: LAB values (col 1:3) Std dev (col 4:6)

	lab_cam_tbl: 10 samples
	10 rows X 6 cols: LAB values (col 1:3) Std dev (col 4:6)

	DE_tbl: 10 samples
	10 rows x 2 cols: DE and std dev on DE

	t_spectro_repro:
	401 rows X 2 col: Lambda (col 1), Std dev on reproducted experiments

	t_cam_repro:
	41 rows X 2 col: Lambda (col 1), Std dev on reproducted experiments

	lab_spectro_repro:
	1 row x 3 cols: Std dev on reproducted experiments

	lab_cam_repro:
	1 row x 3 cols: Std dev on reproducted experiments	

	DE_repro:
	1 row x 1 col: Std dev on reproducted experiments

* expanded_uncert calls f_expanded_uncert and  f_lab_array_tbl

	Applies the unceratinty issued from the reproducibility experiements to the filters
	save results in "Results_KW_Filters" in csv and matlab formats

	Names_KW_ColFilters.txt
	filter_list (Names_KW_ColFilters)
        
	t_spectro_cmp(.txt): 401 x 26; Lambda, (T, std dev, repro uncertainty, total Type A, Expanded (k=2)) x 5 (samples)
    	t_cam_cmp(.txt): 401 x 26; Lambda, (T, std dev, repro uncertainty, total Type A, Expanded (k=2)) x 5 (samples)
    	
	DE_cmp (DE_KW_ColFilters.txt): 5 (samples) x 5; value, std dev, repro uncertainty, total Type A, Expanded (k=2)
    	lab_spectro_cmp (LAB_spectro_KW_ColFilters.txt): 5 (samples) x 15 (value, std dev, repro uncertainty, total Type A, Expanded (k=2)) x 3 (LAB coord)

    	lab_cam_cmp (LAB_cam_KW_ColFilters.txt): 5 (samples) x 15 (value, std dev, repro uncertainty, total Type A, Expanded (k=2)) x 3 (LAB coord)

	lab_array_tbl(.txt): (676x844) x 3 x 5, LAB_array x n_filters

	cov_lab_array_tbl(.txt): 3 x 3 x (676x844) x 5, CovLAB_array x n_filters

* plot_all_t_wcc: 
	plots all transmittance curves, T and CIELAB for KW 32 and KW 47, boxplot of DeltaE values (pixel vs spectro) for all filters

