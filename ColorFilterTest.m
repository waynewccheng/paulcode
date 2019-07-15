
function ColorFilterTest

% create a list based on the file in the folder
% or just enter the list
%    [3 24 34 39 40 43 46 48 51 52 56 57 59 60 68 71 77 85 88 91 93 95 97 99 310 316 336 337 342 347 356 360 389 398]
dirdata = dir('input');
input_file_list = [];
k = 0;
for i = 1:size(dirdata,1)
    filename = dirdata(i).name;
    if strfind(filename,'Filter_') == 1
        k = k + 1;
        str_pos = strfind(filename,'_');
        patch_label = filename(str_pos+1:end);
        patch_no = str2num(patch_label);
        input_file_list(k) = patch_no;
    end
end

%do_it('show3Spectra');

%do_it('showSpectraMean');

% do_it('showTransmittance');

% do_it('showSPD');

show34colors();


return

    function show34colors
        %% show the color filters in CIELAB
        
        labarray = zeros(34,3);
        kk = 1;
        for i = input_file_list
            cf = ColorFilter(i);
            
            labarray(kk,:) = cf.lab;
            
            kk = kk + 1;
        end
        save('CIELAB','labarray')
        
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
        
        saveas(gcf,sprintf('34 colors in CIELAB.tif'))
        
        close all;        
    end

    function do_it (funcname)
        
        figure('units','normalized','outerposition',[0 0 1 1])
        hold on
        kk = 1;
        for i = input_file_list
            cf = ColorFilter(i);
            
            subplot(5,7,kk)
            %    cf.show3Spectra();
            % cf.showSpectraMean();
            cmd =sprintf('cf.%s();',funcname);
            eval(cmd);
            title(sprintf('%d',i))
            
            kk = kk + 1;
        end
        
        saveas(gcf,sprintf('%s.tif',funcname))
        
        close all;
    end


end

