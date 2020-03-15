clear; close all; clc; 

Nps = [1] ; 
countries = {'GER', 'US'} ; 
samples = {'rec'} ; 
surveys = {'level'} ; 
truegdps = {'first', 'second'} ; 

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
                    f_eval(flag_survey, flag_sample, flag_truegdp, Np, flag_country) ;
                end
            end
        end
    end
end