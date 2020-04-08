% Gamma correction on  rgb_linear images
function rgb = f_gamma(rgb_l)
    % Gamma correction
    a=0.055;

    % conditional mask
    rgb_lessorequal = (rgb_l <= 0.0031308);

    % conditional assignment
    rgb_l(rgb_lessorequal) = rgb_l(rgb_lessorequal) * 12.92;
    rgb_l(~rgb_lessorequal) = (1+a)*(rgb_l(~rgb_lessorequal).^(1/2.4)) - a;

    % comply with the old form
    rgb = rgb_l';
end

