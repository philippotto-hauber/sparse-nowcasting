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
Nr = 5; % number of factors

% load results_eval mat file
load(['results_eval_' flag_country '_' flag_sample '_' flag_survey '_Np' num2str(Np) '_' flag_release '.mat'])

% calculate old and new RMSFE for benchmark B-AR and model, plotting the
% two distributions of the RMSFE
figure; 
subplot(1,2,1)
plot_title = ['B-AR'];
[rmsfe_bar_old, rmsfe_bar_new, rmsfe_bar_distr] = f_calcRMSFE(results_eval, Nh, Nprior, Nr, 1, 0, 1, plot_title);

subplot(1,2,2)
plot_title = [results_eval.priors(Nprior).name]; 
[rmsfe_old, rmsfe_new, rmsfe_distr] = f_calcRMSFE(results_eval, Nh, Nprior, Nr, 0, 0, 1, plot_title);

function [rmsfe_old, rmsfe_new, rmsfe_distr] = f_calcRMSFE(results_eval, Nh, Nprior, Nr, flag_bar, flag_pool, flag_plot, plot_title)
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
    rmsfe_distr = sqrt(mean(sfe_draws, 1)); % mean across forecast periods
    rmsfe_new = mean(rmsfe_distr); 

    % plot
    if flag_plot == 1
        hist(rmsfe_distr, 100)
        hold on
        line([mean(rmsfe_distr) mean(rmsfe_distr)], ylim, 'Color','black','LineStyle','-', 'LineWidth', 3) 
        line([median(rmsfe_distr) median(rmsfe_distr)], ylim, 'Color','black','LineStyle',':', 'LineWidth', 3)
        title(plot_title)
        legend('RMSFE', 'mean', 'median')
    end    
end





