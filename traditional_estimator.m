function [alpha_hat, gamma_hat, fail] = traditional_estimator(k1, k2, L, distribution, options)
% TRADITIONAL_ESTIMATOR  estimates alpha and gamma according to the
% traditional Log Cumulant based estimation, which relies on inverting teh
% trigamma function using an out-of-the-shelf solver. A boolean is returned 
% in case of an estimation failure.

    try
        if(isequal(distribution, 'Gi'))
            [alpha_hat, ~, flag] = fsolve(@(alpha) LCum_Gi(alpha, k2, L), -1.00001, options);
            gamma_hat = exp(k1 - psi(0, L) + psi(0,-alpha_hat))*L;

        else
            [alpha_hat, ~, flag] = fsolve(@(alpha) LCum_Ga(alpha, k2, L), -1.00001, options);
            gamma_hat = exp(2*k1 - psi(0, L) + psi(0,-alpha_hat))*L;
        end
        
        if (flag <=0 || alpha_hat > 0 || alpha_hat < -15)
             fail = 1;
             alpha_hat = 0;
             gamma_hat = 0;
        else
            fail = 0;
        end
        
    catch
        fail = 1;
        alpha_hat = 0;
        gamma_hat = 0;
    end

end

function F = LCum_Gi(x, k2, L)
    F = psi(1, L) + psi(1, -x) - k2;
end

function F = LCum_Ga(x, k2, L)
    F = psi(1, L) + psi(1, -x(1)) - 4*k2;
end