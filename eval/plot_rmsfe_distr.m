clear; close all; clc;
% This code plots density-estimates of the distributions of the RMSFE
% for the benchmark B-AR, the NG and a sparse prior. This is to illustrate
% that sparsity does not appear to play a large role in nowcasting
% applications in terms of performance. Code-wise there is an overlap with
% compare_new_old_rmsfe.m! 

% spec details
flag_release = 'first'; 
flag_sample = 'rec'; 
flag_survey = 'level'; 
flag_countries = {'GER', 'US'}; 
Np = 1; % number of lags in factor VAR
Nh = 2; % horizon = 2 months
Nrs = [1, 4, 8]; % number of factors
flag_pool = 0; % plot pool instead of Nr
Npriors = [1, 4, 5]; % priors 

figure; 
subplot_counter = 1;
for f = 1:length(flag_countries)
    flag_country = flag_countries{f};
    
    % load results_eval mat file
    dirname = 'C:\Users\Philipp\Documents\Dissertation\sparse nowcasting\documentation\results_eval mat files\';
    load([dirname 'results_eval_' flag_country '_' flag_sample '_' flag_survey '_Np' num2str(Np) '_' flag_release '.mat'])
    
    flag_pool = 0;
    for Nr = Nrs
        % plot
        subplot(length(flag_countries), length(Nrs)+1, subplot_counter)
        f_plot_dens(results_eval, Nh, Npriors, Nr, flag_country, flag_pool)
        subplot_counter = subplot_counter + 1;
    end
    
    flag_pool = 1;
    % plot
    subplot(length(flag_countries), length(Nrs)+1, subplot_counter)
    f_plot_dens(results_eval, Nh, Npriors, Nr, flag_country, flag_pool)
    subplot_counter = subplot_counter + 1;
end
    
function f_plot_dens(results_eval, Nh, Npriors, Nr, flag_country, flag_pool)
    % set-up
    n = 2^6 ; 
    priors_names = {'NIG', 'MG', 'PMNM', 'HS+', 'Nd'};
    colors = [[1 .5 0]; [0.6 0 0.6]; [0 0 1];[0 0.5 0.5]; [0.9 0 0]];    
    legend_entries = {'B-AR'};

    % calc rmsfe distr
    rmsfe_bar = f_draws_rmsfe(results_eval, Nh, [], Nr, 1, flag_pool);    
    
    for p = 1:length(priors_names) % load all priors, then select when plotting
        rmsfe_prior(p, :) = f_draws_rmsfe(results_eval, Nh, p, Nr, 0, flag_pool);
    end
    
    % range for densities
    max_global = max([max(rmsfe_bar), max(rmsfe_prior,[], 'all')]); 
    min_global = min([min(rmsfe_bar), min(rmsfe_prior,[], 'all')]); 
    range_global = max_global - min_global ; 
    
    % plot kernel density esimates
    [ ~ , ys , xs , ~ ] = kde(rmsfe_bar, n, min_global-range_global/10, max_global+range_global/10);
    plot(xs, ys, 'Color',[0 0 0], 'LineWidth', 1.2)    
    hold on
    for p = 1:length(priors_names)
        if any(Npriors == p)
            [ ~ , ys , xs , ~ ] = kde(rmsfe_prior(p, :), n, min_global-range_global/10, max_global+range_global/10);
            plot(xs, ys, 'Color', colors(p,:), 'LineWidth', 1.2)
            legend_entries = [legend_entries, priors_names{p}];
        end
    end
    legend(legend_entries)
    if flag_pool == 1
        title([flag_country, ', pool'])
    else
        title([flag_country, ', R=' num2str(Nr)])
    end
end
    
function rmsfe_distr = f_draws_rmsfe(results_eval, Nh, Nprior, Nr, flag_bar, flag_pool)
    % get draws of squared errors
    if flag_bar == 1
        sfe_draws = results_eval.benchmark_BAR.horizon(Nh).sfe;
    elseif flag_pool == 1
        sfe_draws = results_eval.priors(Nprior).pool.horizon(Nh).sfe;
    else
        sfe_draws = results_eval.priors(Nprior).R(Nr).horizon(Nh).sfe;
    end
    % RMSFE
    rmsfe_distr = sqrt(mean(sfe_draws, 1)); % mean across forecast periods, then square root
end