%% input 
% image1Name: ROI for viewer 1
% image2Name: ROI for viewer 2

% 03-26-20: use real measurement of CIELAB for Truth inmstead of conversion
% from sRGB to lab

%% output
% diffMatix: Euclidean distance in the CIELAB space 
function [ diffMatrix ] = compareColorFunction_LAB2( imgTestName, labTruth)

    % Read image
    image1 = imread(imgTestName);
    
    % Convert Test img to CIELAB
    CIELAB1 = rgb2lab(image1); 
    
    % Compare CIELAB values
    [row,col,channel] = size(CIELAB1);
    CIELAB1vector = reshape(CIELAB1,row*col,[]);
     
    diff2D = (CIELAB1vector - labTruth).^2;
    diff1D = sqrt(sum(diff2D,2));
    diffMatrix = reshape(diff1D,row,col);
%     save('dE_2.mat','diffMatrix');

end

