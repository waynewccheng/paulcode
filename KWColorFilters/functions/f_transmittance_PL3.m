% Substract the background light from the sample and the 100% transmittance
function [t_m, t_s] = f_transmittance_PL3(s_m, w_m, s_bg_m, w_bg_m, s_s, w_s, s_bg_s, w_bg_s)
    % Compute the transmittace value and the uncertainty
    
    % Transmittance
    t_m = (s_m - s_bg_m )./ (w_m - w_bg_m);
    
    % Uncertainty
    t_s = sqrt((1./(w_m - w_bg_m)).^2.* s_s.^2 + ...
        (-1./(w_m - w_bg_m)).^2.* s_bg_s.^2 + ...
        ((s_m - s_bg_m)./(w_m - w_bg_m).^2).^2.* w_bg_s.^2 + ...
        (-(s_m - s_bg_m)./(w_m - w_bg_m).^2).^2 .* w_s.^2);
    
end

