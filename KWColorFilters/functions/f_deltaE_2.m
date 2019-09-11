% 09-11-19: Log-normal distribution for DeltaE

function [dE, s_dE] = f_deltaE_2(LAB_1, LAB_2, CovLAB_1, CovLAB_2)
    
    % Computes DeltaE and its uncertainty
    DeltaE = @(LAB_1, LAB_2) sqrt((LAB_1(1) - LAB_2(1)).^2 + ...    % L
    (LAB_1(2) - LAB_2(2)).^2 + ...                                  % a
    (LAB_1(3) - LAB_2(3)).^2);                                      % b

    dE = DeltaE(LAB_1, LAB_2);
    
    mu_log_dE = log(dE);

    v_log_dE = 1./dE.^4*[(LAB_1(1) - LAB_2(1)) (LAB_1(2) - LAB_2(2)) (LAB_1(3) -...
        LAB_2(3)) -(LAB_1(1) - LAB_2(1)) -(LAB_1(2) - LAB_2(2)) -(LAB_1(3) - LAB_2(3))] *...
        blkdiag(CovLAB_1, CovLAB_2)*...
        [(LAB_1(1) - LAB_2(1)) (LAB_1(2) - LAB_2(2)) (LAB_1(3) - LAB_2(3)) -...
        (LAB_1(1) - LAB_2(1)) -(LAB_1(2) - LAB_2(2)) -(LAB_1(3) - LAB_2(3))]';
    
    v_dE = (exp(v_log_dE) -1)*exp(2*mu_log_dE+v_log_dE);
    s_dE = sqrt(v_dE);

end

