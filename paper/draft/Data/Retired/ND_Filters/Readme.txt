Results Format:

Names_Kodak_ND: name of data folder set during measurements

Lambda: 380nm to 780nm in steps of 10nm

T_Spectro_m: temporal mean over 10 measured spectra
41 wavelengths x 6 neutral density filters

T_Spectro_s: temporal std dev over 10 measured spectra, total uncertainty
41 wavelengths x 6 neutral density filters

T_Spectro_cmp: T Spectro data, Lambda then (value, 1 X sigma, 2 X sigma) x 6 (18 col)

Cam_m: camera spatial averages, mean values over 10 measured spectra
41 wavelengths x 6 neutral density filters

Cam_s: camera spatial averages, std dev values over 10 measured spectra, total uncertainty
41 wavelengths x 6 neutral density filters

Cam_T_cmp: cam data, Lambda then (value, 1 X sigma, 2 X sigma) x 6 (18 col)

LAB_spectro_Kodak_ND, LAB_cam_Kodak_ND: 6 neutral density filters x 15 columns
col 1:3: L, a, b
col 4:6: uncertainty on L, a, b (repetability)
col 7:9: reproducibility uncertainties on L, a, b
col 10:12: total uncertainty (sqrt of sum of square of col 4:6 and col 7:9)
col 13:15: expanded uncertainty: 2 * col 10:12

DE_Kodak_ND: 6 neutral density filters x 5 columns
col 1: Delta E
col 2: uncertainty on Delta E
col 3: reproducibility uncertainties on Delta E
col 4: total uncertainty (sqrt of sum of square of col 2 and col 3)
col 5: expanded uncertainty: 2 * col 4

Lin_reg_cam: linear regression parameters for each wavelength, the spectrometer mesasurements is the truth, it's a weighed linear interpolation using the total uncertainty on the camera measurements
col 1: slope
col 2: intercept
col 3: rmse



