% 7-19-2016
% D65
% 7-23-2015
% convert reflectance into RGB using D65
% usage: rgb = reflectance2D65(reflectance_array7);

function [LAB, CovLAB, XYZ, CovXYZ] = transmittance2LAB_PL2(trans_array_m, trans_array_s, sizey, sizex, ls)

    disp('Combining reflectance and illuminant into LAB...')

    % reference white
    XYZ0 = spd2XYZ_PL2(1, 0, ls, 'white');

    % Image dimensions here
    if sizey ~= 1
        ls_array = repmat(ls,1,sizey*sizex);
    else
        ls_array = ls;
    end
    
    disp('  calculate XYZ...')
    [XYZ, CovXYZ] = spd2XYZ_PL2(trans_array_m, trans_array_s, ls_array, 'sample');

    disp('  calculate LAB...')
    [LAB, CovLAB] = XYZ2lab(XYZ, XYZ0, CovXYZ);

end