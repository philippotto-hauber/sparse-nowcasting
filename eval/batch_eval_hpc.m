function batch_eval_hpc(Nspec, country)
if (isdeployed)
    maxNumCompThreads(1);
    Nspec = str2double(Nspec) ;
    country = str2double(country);
end

% ----------------------------
% - spec details -------------
% ----------------------------
% back out flag_sample, survey and truegdp from Nspec
Nsurveys = {'level', 'diff'};
Nsamples = {'rec', 'rolling'};
Ntruegdps = {'first', 'second', 'final'};
Nps = [1, 3];

tmp = repmat(Nsurveys, length(Nsamples) * length(Ntruegdps) * length(Nps), 1); tmp_survey = tmp(:);
tmp = repmat(Nsamples, length(Ntruegdps) * length(Nps), length(Nsurveys)); tmp_sample = tmp(:);
tmp = repmat(Ntruegdps, length(Nps), length(Nsurveys) * length(Nsamples)); tmp_truegdp = tmp(:);
tmp = repmat(Nps, 1, length(Ntruegdps) * length(Nps) *length(Nsurveys) * length(Nsamples)); tmp_Np = tmp(:); 

% back out flags
flag_survey = tmp_survey{Nspec}; 
flag_sample = tmp_sample{Nspec};
flag_truegdp = tmp_truegdp{Nspec}; 
Np = tmp_Np(Nspec); 

% flag_country
if country == 0
    flag_country = 'GER'; 
elseif country == 1
    flag_country = 'US'; 
end
set(0,'DefaultFigureVisible','off') 


% ----------------------------
% - switches -----------------
% ----------------------------
switch_eval_matfiles = 'off' ; 
switch_eval_tables = 'off' ; 
switch_eval_pointnowcasts = 'off' ; 
switch_eval_cumsum = 'off' ;
switch_eval_densities = 'off' ; 
switch_eval_latextables = 'on'; 

% ----------------------------
% - call functions -----------
% ----------------------------
                    
% call f_eval
if strcmp(switch_eval_matfiles, 'on') 
f_eval_hpc(flag_survey, flag_sample, flag_truegdp, Np, flag_country) ;
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

% call f_eval_plot_densities
if strcmp(switch_eval_densities, 'on') 
f_eval_plot_densities(flag_survey, flag_sample, flag_truegdp, Np, flag_country) 
end

% call f_eval_latex_table
if strcmp(switch_eval_latextables, 'on') 
f_eval_latex_table(flag_survey, flag_sample, flag_truegdp, Np, flag_country)
end
end
