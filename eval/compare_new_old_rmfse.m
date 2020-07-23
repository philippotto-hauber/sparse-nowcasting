clear; close all; clc; 
%------------------------------------------------------------------------ %
%- This function compares two ways of calculating the RMSFE. The old way 
%- calculated the average squared error for each period, then averaged
%- over the evaluation period and then takes the square root. 
%- The "new" way calculates the RMSFE for each draw and then takes the 
%- average of that distribution as the reported value in the tables.
%- Furthermore, it also plots the resulting distribution of the RMSFE for
%- a given model spec and the benchmark B-AR(1).
%------------------------------------------------------------------------ %

% spec details
flag_release = 'first'; 
flag_sample = 'rec'; 
flag_survey = 'level'; 
flag_country = 'GER'; 
Np = 1; % number of lags in factor VAR
Nprior = 4; % prior
Nh = 2; % horizon = 2 months
Nr = 8; % number of factors

% load results_eval mat file
dirname = 'C:\Users\Philipp\Documents\Dissertation\sparse nowcasting\documentation\results_eval mat files\';
load([dirname 'results_eval_' flag_country '_' flag_sample '_' flag_survey '_Np' num2str(Np) '_' flag_release '.mat'])

% calculate old and new RMSFE for benchmark B-AR and model, plotting the
% two distributions of the RMSFE
figure; 
[rmsfe_bar_old, rmsfe_bar_new, rmsfe_bar_distr] = f_calcRMSFE(results_eval, Nh, Nprior, Nr, 1, 0, 1);
hold on
[rmsfe_old, rmsfe_new, rmsfe_distr] = f_calcRMSFE(results_eval, Nh, Nprior, Nr, 0, 0, 1);

% add legend and title to plot
legend('B-AR', results_eval.priors(Nprior).name)
title(['posterior RMSFE: ' flag_country ',' flag_sample ',' flag_survey ',R=' num2str(Nr)]) 

% calculate relative RMSFE
rel_rmsfe_old = rmsfe_old / rmsfe_bar_old; 
rel_rmsfe_new = rmsfe_new / rmsfe_bar_new; 


function [rmsfe_old, rmsfe_new, rmsfe_distr] = f_calcRMSFE(results_eval, Nh, Nprior, Nr, flag_bar, flag_pool, flag_plot)
    % get draws of squared errors
    if flag_bar == 1
        sfe_draws = results_eval.benchmark_BAR.horizon(Nh).sfe;
    elseif flag_pool == 1
        sfe_draws = results_eval.priors(Nprior).pool.horizon(Nh).sfe;
    else
        sfe_draws = results_eval.priors(Nprior).R(Nr).horizon(Nh).sfe;
    end

    % calculate "old" RMSFE
    mean_sfe = mean(sfe_draws, 2); % mean squared error across draws for all period
    mse = mean(mean_sfe); 
    rmsfe_old = sqrt(mse);

    % new RMSFE
    rmsfe_distr = sqrt(mean(sfe_draws, 1)); % mean across forecast periods, then square root
    rmsfe_new = mean(rmsfe_distr); % mean over draws -> posterior mean!

    % plot
    if flag_plot == 1
        histogram(rmsfe_distr, 100)               
    end    
end





