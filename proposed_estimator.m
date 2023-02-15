function [alpha_hat, gamma_hat, fail] = proposed_estimator(k1, k2, m4, n, L, c_alpha, do_correction)
% PROPOSED_ESTIMATOR  estimates alpha and gamma accordingto the
% appprximation proposed in the paper. The approximation may be accompanied
% by a bayesian correction. A boolean is returned in case of an estimation
% failure.
    
    % Computei initial estimation for eta
    eta_hat = c_alpha * k2 - psi(1,L);
    
    % Apply correction or skip
    if do_correction
        sigma_normal = sqrt((c_alpha^2)/(n)*(m4 - ((n-3)/(n-1)) * k2^2)); % k2 == variance
        
        c = normcdf(eta_hat/sigma_normal);
        eta_post = eta_hat + sigma_normal*exp(-(eta_hat^2)/(2*sigma_normal^2))/...
                             (sqrt(2*pi)*c);
    else
        eta_post = eta_hat;
    end
    
    % Detect failure (eta < 0)
    if eta_post < 0
        fail = 1;
        alpha_hat = 0; 
        gamma_hat = 0;
        return
    end
    
    % Invert trigamma function
    alpha_hat = poly(eta_post);
    
    % Detect failures according to paper
    if(~isreal(alpha_hat) || alpha_hat >= 0 || alpha_hat <= -15)
        fail = 1;
        alpha_hat = 0;
        gamma_hat = 0;
    else
        % If not a failure, estimate gamma
        fail = 0;
        gamma_hat = exp(sqrt(c_alpha) * k1 - psi(0, L) + psi(0,-alpha_hat))*L;
    end
end

function alpha_hat = poly(eta)
    p = [-eta*210 210 105 35 0 -7 0 5];
    r = roots(p);
    alpha_hat = -r(1);
end