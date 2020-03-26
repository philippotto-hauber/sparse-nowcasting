clear; close all; clc; 

% ----------------------------
% - spec details -------------
% ----------------------------
Nps = [1] ; 
countries = {'GER', 'US'} ; 
samples = {'rec'} ; 
surveys = {'level'} ; 
truegdps = {'first'} ; 


% ----------------------------
% - switches -----------------
% ----------------------------
switch_eval_matfiles = 'off' ; 
switch_eval_tables = 'off' ; 
switch_eval_pointnowcasts = 'on' ; 
switch_eval_cumsum = 'off' ; 
switch_eval_densities = 'off' ; 
switch_eval_latextables = 'off'; 

% ----------------------------
% - loop -----------------
% ----------------------------
for i1 = 1:length(samples) 
    flag_sample = samples{i1} ; 
    for i2 = 1:length(surveys) 
        flag_survey = surveys{i2} ; 
        for Np = Nps
            for i3 = 1:length(countries)
                flag_country = countries{i3} ; 
                for i4 = 1:length(truegdps) 
                    flag_truegdp = truegdps{i4} ; 
                    
                    % call f_eval
                    if strcmp(switch_eval_matfiles, 'on') 
                        f_eval(flag_survey, flag_sample, flag_truegdp, Np, flag_country) ;
                    end
                    
                    % call f_eval_tables
                    if strcmp(switch_eval_tables, 'on') 
                        f_eval_tables(flag_survey, flag_sample, flag_truegdp, Np, flag_country) 
                    end
                    
                    % call f_eval_plot_pointnowcasts
                    if strcmp(switch_eval_pointnowcasts, 'on') 
                        f_eval_plot_pointnowcasts(flag_survey, flag_sample, flag_truegdp, Np, flag_country) 
                    end
                    
                    % call f_eval_plot_cumsum
                    if strcmp(switch_eval_cumsum, 'on') 
                        f_eval_plot_cumsum(flag_survey, flag_sample, flag_truegdp, Np, flag_country) 
                    end
                    
                    % call f_eval_plot_cumsum
                    if strcmp(switch_eval_densities, 'on') 
                        f_eval_plot_densities(flag_survey, flag_sample, flag_truegdp, Np, flag_country) 
                    end
                    
                    % call f_eval_latex_table
                    if strcmp(switch_eval_latextables, 'on')
                        f_eval_latex_table(flag_survey, flag_sample, flag_truegdp, Np, flag_country)
                    end
                end
            end
        end
    end
end