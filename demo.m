%% Experimental parameters

num_montecalo_samples = 1000;

sample_size = 9;

distribution = "Gi"; % Choose between 'Ga' and 'Gi'
alpha = -1.5;
gamma = -alpha-1;
L = 1;

%% Functions to generate samples form Ga and Gi (ref. paper)

Gi = @(alpha, gamma, L, n_samp)...
    -(gamma/alpha)*finv(rand(n_samp,1),2*L,-2*alpha);

Ga = @(alpha, gamma, L, n_samp)...
    sqrt(-(gamma/alpha)*finv(rand(n_samp,1),2*L,-2*alpha));

%% Run Montecarlo experiments for each available estimator

estimator = ["TradLCum", "PropNonCorrected", "PropCorrected"]';

mse_alpha = zeros(size(estimator));
failure_rate = zeros(size(estimator));
average_time = zeros(size(estimator));

fsolve_options = optimoptions('fsolve','Display','none');
for i_est = 1:length(estimator)    
    c = 0; alpha_hat = zeros(num_montecalo_samples,1); gamma_hat = zeros(num_montecalo_samples,1); time = 0;
    for i_mc = 1:num_montecalo_samples
        % Generate samples and set c_alpha accoridng to the distribution
        switch distribution
            case "Gi"
                samples = Gi(alpha, gamma, L, sample_size);
                c_alpha = 1;
            case "Ga"
                samples = Ga(alpha, gamma, L, sample_size);
                c_alpha = 4;
            otherwise
                error("Wrong distribution name. Choose 'Ga' or 'Gi'.")
        end

        tic
        
        % Compute log-cumulats from samples
        z_samples = log(samples);
        k1 = mean(z_samples);
        k2 = mean((z_samples - k1).^2);
        
        % Estimate alpha and gamma accoridng to the desired estimator
        switch estimator(i_est)
            case "TradLCum"
                [a, g, fail] = traditional_estimator(k1, k2, L, distribution, fsolve_options);
            case "PropNonCorrected"
                m4 = mean(z_samples.^4); % Compute fouth order moment
                do_correction = false;
                
                % Run proposed estimator without correction
                [a, g, fail] = proposed_estimator(k1, k2, m4, sample_size, L, c_alpha, do_correction);
            case "PropCorrected"
                m4 = mean(z_samples.^4); % Compute fouth order moment
                do_correction = true;
                
                % Run proposed estimator with correction
                [a, g, fail] = proposed_estimator(k1, k2, m4, sample_size, L, c_alpha, do_correction);
            otherwise
                error("Wrong estimator name. Choose 'Pptim', 'PropNonCorrected' or 'PropCorrected'.")
        end

        time = time + toc;
        
        % Only keep the non-failure estimates (for computing MSE)
        if (~fail)
            c = c + 1;
            alpha_hat(c) = a;
            gamma_hat(c) = g;
        end
    end
    alpha_hat = alpha_hat(1:c); gamma_hat = gamma_hat(1:c);
    
    mse_alpha(i_est) = mean((alpha_hat-alpha).^2);
    failure_rate(i_est) = 1-(c/num_montecalo_samples);
    average_time(i_est) = time/num_montecalo_samples;
end

%% Display MSE, failure rate and average runtime

fprintf("<strong> Estimation results for: %s distribution, alpha = %.1f, L = %d, sample size = %d\n\n </strong>", ...
    distribution, alpha, L, sample_size)
disp(table(estimator, mse_alpha, failure_rate, average_time))
