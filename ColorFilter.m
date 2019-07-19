classdef ColorFilter < handle
    
    properties
        id
        spectra
        spectra_mean
        transmittance
        spd
        spd_d65
        XYZ
        XYZ_d65
        lab
        wavelength_range_pr730 = 380:1:780;
        vim_mean_array_sample
        vim_mean_array_white
        vim_mean_array_black
        transmittanceMSI

        % where to save the files?
        finding_path = 'finding';
    end
    
    properties (Constant)
        list_patch = [3 24 34 39 40 43 46 48 51 52 56 57 59 60 68 71 77 85 88 91 93 95 97 99 310 316 336 337 342 347 356 360 389 398];
    end

    methods
        
        function obj = ColorFilter (patch_id)
            obj.id = patch_id;
            
            obj.getDataPr730(patch_id);
            
            
            obj.getDataMSI(patch_id);
            
            % obj.show3Spectra();
                        
            
            obj.calcSpectraMean();
            
            % obj.showSpectraMean();
            
            
            obj.calcTransmittance();
            
            % obj.showTransmittance();
            
            
            obj.calcSPD();
            
            % obj.showSPD();
            
            obj.calcXYZ();
            
            obj.calcLAB();
            
        end

        function do_msi (obj)
            
            obj.showDataMSI('transmittanceMSI')
            obj.showDataMSI('vim_mean_array_sample')
            obj.showDataMSI('vim_mean_array_white')
            obj.showDataMSI('vim_mean_array_black')
            
        end
                
        function do_truth (obj)
           
            figure('units','normalized','outerposition',[0 0 1 1])
            
            subplot(2,2,1)
            obj.show3Spectra();
            title('All, Digital')
            
            subplot(2,2,2)
            obj.showSpectraMean();
            title('Mean, Digital')

            subplot(2,2,3)
            obj.showSPD();
            title('SPD by D65')

            subplot(2,2,4)
            obj.showTransmittance();
            title('Transmittance')

            saveas(gcf,sprintf('%s/%s.tif',obj.finding_path,'pr730'))            
            
            close all
        
        end
        
        function showSPD (obj)
            hold on
            plot(obj.wavelength_range_pr730,obj.spd)
            axis([380 780 0 100])
        end
        
        function calcSPD (obj)
            obj.spd_d65 = obj.getD65();
            
            obj.spd = obj.transmittance .* obj.spd_d65;
        end
        
        function spec = getD65 (obj)
            load('spec_cied65','spec')
            spec = spec(:,2);
        end
        
        function showSpectraMean (obj)
            hold on
            plot(obj.wavelength_range_pr730,obj.spectra_mean(:,1),'-r')
            plot(obj.wavelength_range_pr730,obj.spectra_mean(:,2),'-g')
            plot(obj.wavelength_range_pr730,obj.spectra_mean(:,3),'-b')
            % too large for 34 subplots
            % legend('Sample','100%','0%')
            %axis([380 780 0 8e-4])
        end
        
        function calcTransmittance (obj)
            sample = obj.spectra_mean(:,1);
            white = obj.spectra_mean(:,2);
            black = obj.spectra_mean(:,3);
            
            obj.transmittance = (sample-black) ./ (white-black);
        end
        
        function showTransmittance (obj)
            plot(obj.wavelength_range_pr730,obj.transmittance)
            axis([obj.wavelength_range_pr730(1) obj.wavelength_range_pr730(end) 0 1])
        end
        
        function calcSpectraMean (obj)
            obj.spectra_mean = squeeze(mean(obj.spectra,1));
        end

        function getDataMSI (obj,patch_id)
            filename = sprintf('./input/Filter_%02d/Filter_%02d_sample/vim_mean_array.mat',patch_id,patch_id)
            load(filename,'vim_mean_array')
            obj.vim_mean_array_sample = vim_mean_array;

            filename = sprintf('./input/Filter_%02d/Filter_%02d_white/vim_mean_array.mat',patch_id,patch_id)
            load(filename,'vim_mean_array')
            obj.vim_mean_array_white = vim_mean_array;

            filename = sprintf('./input/Filter_%02d/Filter_%02d_black/vim_mean_array.mat',patch_id,patch_id)
            load(filename,'vim_mean_array')
            obj.vim_mean_array_black = vim_mean_array;
            
            obj.transmittanceMSI = ...
                100 * ...
                (obj.vim_mean_array_sample - obj.vim_mean_array_black)./ ...
                (obj.vim_mean_array_white - obj.vim_mean_array_black);
            
        end

        %% show heatmap
        function showDataMSI (obj,property_name)
            
            figure('units','normalized','outerposition',[0 0 1 1])

            wl_array = 380:10:780;
            wl_weighting = ColorFilter.colorMatchingFunctionWeighting();
            
            for wl = 1:41
                subplot(6,7,wl)
            
                vvname = sprintf('obj.%s(wl,:,:)',property_name);
                vv = eval(vvname);
                im = squeeze(vv);
                imagesc(im)

                % calculate the mean and std
                vv1 = reshape(im,size(im,1)*size(im,2),1);
                vv_mean = mean(vv1);
                vv_std = std(vv1);
                
                % axis off
                set(gca,'xtick',[])
                set(gca,'xticklabel',[])
                set(gca,'ytick',[])
                set(gca,'yticklabel',[])
                axis image
                title(sprintf('%d (%.2f%%)',wl_array(wl),wl_weighting(wl)*100))
                xlabel(sprintf('%.4f, %.4f',vv_mean,vv_std))
                colorbar
            end
            
            % show title of titles
            subplot(6,7,42)
            title(sprintf('%s - %d',property_name,obj.id),'Interpreter','none')
            axis off
            
            saveas(gcf,sprintf('%s/%s.tif',obj.finding_path,property_name))            
            close all
        end
        
        function getDataPr730 (obj,patch_id)
            filename = sprintf('input/Filter_%02d/Filter_%02d_spectro_8xFast/spectro_meas',patch_id,patch_id);
            load(filename,'spectra')
            obj.spectra = spectra;
        end
        
        function show3Spectra (obj)
            hold on
            for shot = 1:10
                plot(obj.wavelength_range_pr730,obj.spectra(shot,:,1),'-r')
                plot(obj.wavelength_range_pr730,obj.spectra(shot,:,2),'-g')
                plot(obj.wavelength_range_pr730,obj.spectra(shot,:,3),'-b')
                % too large when showing 34 subplots
                % legend('Sample','100%','0%')
                %axis([380 780 0 8e-4])
            end
        end
        
        %% convert spectrum into XYZxyz
        % spd is a 41xN matrix
        % each spd is 41x1; no wavelength
        %
        function calcXYZ (obj)
            
            spd = obj.spd(1:10:end);
            
            % CIEXYZ 1931
            cmf = [
                380.0 0.001368 0.000039 0.006450;
                390.0 0.004243 0.000120 0.020050;
                400.0 0.014310 0.000396 0.067850;
                410.0 0.043510 0.001210 0.207400;
                420.0 0.134380 0.004000 0.645600;
                430.0 0.283900 0.011600 1.385600;
                440.0 0.348280 0.023000 1.747060;
                450.0 0.336200 0.038000 1.772110;
                460.0 0.290800 0.060000 1.669200;
                470.0 0.195360 0.090980 1.287640;
                480.0 0.095640 0.139020 0.812950;
                490.0 0.032010 0.208020 0.465180;
                500.0 0.004900 0.323000 0.272000;
                510.0 0.009300 0.503000 0.158200;
                520.0 0.063270 0.710000 0.078250;
                530.0 0.165500 0.862000 0.042160;
                540.0 0.290400 0.954000 0.020300;
                550.0 0.433450 0.994950 0.008750;
                560.0 0.594500 0.995000 0.003900;
                570.0 0.762100 0.952000 0.002100;
                580.0 0.916300 0.870000 0.001650;
                590.0 1.026300 0.757000 0.001100;
                600.0 1.062200 0.631000 0.000800;
                610.0 1.002600 0.503000 0.000340;
                620.0 0.854450 0.381000 0.000190;
                630.0 0.642400 0.265000 0.000050;
                640.0 0.447900 0.175000 0.000020;
                650.0 0.283500 0.107000 0.000000;
                660.0 0.164900 0.061000 0.000000;
                670.0 0.087400 0.032000 0.000000;
                680.0 0.046770 0.017000 0.000000;
                690.0 0.022700 0.008210 0.000000;
                700.0 0.011359 0.004102 0.000000;
                710.0 0.005790 0.002091 0.000000;
                720.0 0.002899 0.001047 0.000000;
                730.0 0.001440 0.000520 0.000000;
                740.0 0.000690 0.000249 0.000000;
                750.0 0.000332 0.000120 0.000000;
                760.0 0.000166 0.000060 0.000000;
                770.0 0.000083 0.000030 0.000000;
                780.0 0.000042 0.000015 0.000000;
                ];
            
            % show color matching functions
            if 0
                clf
                hold on;
                plot(cmf(:,1),cmf(:,2),'Color','r');
                plot(cmf(:,1),cmf(:,3),'Color','g');
                plot(cmf(:,1),cmf(:,4),'Color','b');
                title('Color Matching Functions')
            end
            
            % extend x_bar from 41x1 to 41xN
            input_n = size(spd,2);
            x_bar = repmat(cmf(:,2),1,input_n);
            y_bar = repmat(cmf(:,3),1,input_n);
            z_bar = repmat(cmf(:,4),1,input_n);
            
            X = sum(spd .* x_bar);
            Y = sum(spd .* y_bar);
            Z = sum(spd .* z_bar);
            
            % output Nx3
            XYZ = [X' Y' Z'];
            
            obj.XYZ = XYZ;
            
            %% d65
            % extend x_bar from 41x1 to 41xN
            spd_d65 = obj.spd_d65(1:10:end);
            
            input_n = size(spd_d65,2);
            x_bar = repmat(cmf(:,2),1,input_n);
            y_bar = repmat(cmf(:,3),1,input_n);
            z_bar = repmat(cmf(:,4),1,input_n);
            
            X_d65 = sum(spd_d65 .* x_bar);
            Y_d65 = sum(spd_d65 .* y_bar);
            Z_d65 = sum(spd_d65 .* z_bar);
            
            % output Nx3
            obj.XYZ_d65 = [X_d65' Y_d65' Z_d65'];
            
        end
        
        function calcLAB (obj)
            obj.lab = obj.XYZ2lab (obj.XYZ,obj.XYZ_d65);
        end
        
        %%
        %
        % CIEXYZ to CIELAB
        %
        function ret = XYZ2lab (obj,XYZ, XYZ_white)
            % XYZ is k-by-3
            % XYZ_white is 1-by-3
            k = size(XYZ,1);
            XYZn = repmat(XYZ_white,k,1);
            XYZ_over_XYZn = XYZ./XYZn;
            
            lstar = 116 * helpf(XYZ_over_XYZn(:,2)) - 16;
            astar = 500 * (helpf(XYZ_over_XYZn(:,1)) - helpf(XYZ_over_XYZn(:,2)));
            bstar = 200 * (helpf(XYZ_over_XYZn(:,2)) - helpf(XYZ_over_XYZn(:,3)));
            
            ret = [lstar astar bstar];
            %     if lstar > 100
            %         ['exceeding in xyz2lab']
            %         [x y z xn yn zn lstar astar bstar]
            %         lstar = 100;
            %     end
            
            return
            
            function ys = helpf (t)
                % conditional mask
                t_greater = (t > power(6/29,3));
                
                % conditional assignment
                t(t_greater) = t(t_greater) .^ (1/3);
                t(~t_greater) = t(~t_greater) * (((29/6)^2)/3) + 4/29;
                
                ys = t;
            end
        end
        
        %
        % lan1, lab2: 3xN matrix
        %
        function dEab = LAB2dEab (obj, lab1, lab2)
            % sum in the vertical direction
            dEab = sum((lab2-lab1) .^ 2,1) .^ 0.5;
        end
    end
    
    %% use these static methods to go through the 34 patches
    methods (Static)

        
        %% generate 5 plots for each patch
        function show5plots
            
            ColorFilter.do_it_for_all('show3Spectra');
            ColorFilter.do_it_for_all('showSpectraMean');
            ColorFilter.do_it_for_all('showTransmittance');
            ColorFilter.do_it_for_all('showSPD');
            
            ColorFilter.show34colors_in_CIELAB();
        end
        

        %% generate the specific result for all patches
        function do_it_for_all (funcname)
            
            % where to save the files?
            finding_path = 'finding';
            
            %input_file_list = ColorFilter.getInputFileList();
            input_file_list = ColorFilter.list_patch;
            
            figure('units','normalized','outerposition',[0 0 1 1])
            hold on
            kk = 1;
            for i = input_file_list
                cf = ColorFilter(i);
                
                subplot(5,7,kk)
                
                % things like:
                % cf.show3Spectra();
                % cf.showSpectraMean();
                cmd =sprintf('cf.%s();',funcname);
                eval(cmd);
                title(sprintf('%d',i))
                
                kk = kk + 1;
            end
            
            % show title of titles
            subplot(5,7,kk)
            title(funcname)
            axis off
            
            saveas(gcf,sprintf('%s/%s.tif',finding_path,funcname))
            
            close all;
        end
        
        %% show the color patches in 3D CIELAB
        function show34colors_in_CIELAB
            %% show the color filters in CIELAB
            
            % where to save the files?
            finding_path = 'finding';
            
            input_file_list = ColorFilter.getInputFileList();
            
            labarray = zeros(34,3);
            kk = 1;
            for i = input_file_list
                cf = ColorFilter(i);
                
                labarray(kk,:) = cf.lab;
                
                kk = kk + 1;
            end
            save(sprintf('%s/CIELAB',finding_path),'labarray')
            
            figure('units','normalized','outerposition',[0 0 1 1])
            hold on
            rgb = lab2rgb(labarray);
            rgb = min(rgb,1);
            rgb = max(rgb,0);
            
            for i = 1:size(labarray,1)
                plot3(labarray(i,2),labarray(i,3),labarray(i,1),...
                    'o','Markersize',10,...
                    'MarkerFaceColor',rgb(i,:)',...
                    'MarkerEdgeColor',rgb(i,:)')
            end
            
            xlabel('CIE a*')
            ylabel('CIE b*')
            zlabel('CIE L*')
            grid on
            axis equal
            
            view(25,25)
            
            saveas(gcf,sprintf('%s/34 colors in CIELAB.tif',finding_path))
            
            close all;
        end
        
        function input_file_list = getInputFileList
            % create a list based on the file in the folder
            % or just enter the list
            %    [3 24 34 39 40 43 46 48 51 52 56 57 59 60 68 71 77 85 88 91 93 95 97 99 310 316 336 337 342 347 356 360 389 398]
            dirdata = dir('input');                            % the subfolder
            input_file_list = [];
            k = 0;
            for i = 1:size(dirdata,1)                          % for each file
                filename = dirdata(i).name;                    % get the filename part
                if strfind(filename,'Filter_') == 1            % if it is Filter_XX
                    k = k + 1;
                    str_pos = strfind(filename,'_');           % get the XX part
                    patch_label = filename(str_pos+1:end);
                    patch_no = str2num(patch_label);           % store as a number
                    input_file_list(k) = patch_no;
                end
            end
        end
        
        
        function s = colorMatchingFunctionWeighting
            
            % CIEXYZ 1931
            cmf = [
                380.0 0.001368 0.000039 0.006450;
                390.0 0.004243 0.000120 0.020050;
                400.0 0.014310 0.000396 0.067850;
                410.0 0.043510 0.001210 0.207400;
                420.0 0.134380 0.004000 0.645600;
                430.0 0.283900 0.011600 1.385600;
                440.0 0.348280 0.023000 1.747060;
                450.0 0.336200 0.038000 1.772110;
                460.0 0.290800 0.060000 1.669200;
                470.0 0.195360 0.090980 1.287640;
                480.0 0.095640 0.139020 0.812950;
                490.0 0.032010 0.208020 0.465180;
                500.0 0.004900 0.323000 0.272000;
                510.0 0.009300 0.503000 0.158200;
                520.0 0.063270 0.710000 0.078250;
                530.0 0.165500 0.862000 0.042160;
                540.0 0.290400 0.954000 0.020300;
                550.0 0.433450 0.994950 0.008750;
                560.0 0.594500 0.995000 0.003900;
                570.0 0.762100 0.952000 0.002100;
                580.0 0.916300 0.870000 0.001650;
                590.0 1.026300 0.757000 0.001100;
                600.0 1.062200 0.631000 0.000800;
                610.0 1.002600 0.503000 0.000340;
                620.0 0.854450 0.381000 0.000190;
                630.0 0.642400 0.265000 0.000050;
                640.0 0.447900 0.175000 0.000020;
                650.0 0.283500 0.107000 0.000000;
                660.0 0.164900 0.061000 0.000000;
                670.0 0.087400 0.032000 0.000000;
                680.0 0.046770 0.017000 0.000000;
                690.0 0.022700 0.008210 0.000000;
                700.0 0.011359 0.004102 0.000000;
                710.0 0.005790 0.002091 0.000000;
                720.0 0.002899 0.001047 0.000000;
                730.0 0.001440 0.000520 0.000000;
                740.0 0.000690 0.000249 0.000000;
                750.0 0.000332 0.000120 0.000000;
                760.0 0.000166 0.000060 0.000000;
                770.0 0.000083 0.000030 0.000000;
                780.0 0.000042 0.000015 0.000000;
                ];
            
            
            curves = cmf(:,2:4);
            curves_auc = sum(curves,1);
            curves_auc_rep = repmat(curves_auc,size(cmf,1),1);
            curves_normalized = curves ./ curves_auc_rep;
            
            curves_normalized_sum = sum(curves_normalized,2); % the sum of 3 curves
            curves_normalized_sum_auc = sum(curves_normalized_sum,1);
            s = curves_normalized_sum ./ curves_normalized_sum_auc;
            
            %             clf
            %             plot(cmf(:,1),curves_normalized(:,1),'r')
            %             hold on
            %             plot(cmf(:,1),curves_normalized(:,2),'g')
            %             plot(cmf(:,1),curves_normalized(:,3),'b')
            %
            %             plot(cmf(:,1),s,'m')
        end
    end
    
end
