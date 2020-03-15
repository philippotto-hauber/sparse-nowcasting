clear; close all; clc; 

% ----------------------------
% - spec details -------------
% ----------------------------
Nps = [1] ; 
countries = {'GER', 'US'} ; 
samples = {'rec'} ; 
surveys = {'level'} ; 
truegdps = {'first', 'second'} ; 


% ----------------------------
% - switches -----------------
% ----------------------------
switch_eval_matfiles = 'off' ; 
swtich_eval_tables = 'off' ; 

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
                end
            end
        end
    end
end