clear;clc;

country = 'Germany';

if strcmp(country, 'Germany')
    foldername = 'C:\Users\Philipp\Documents\Dissertation\sparse nowcasting\estim\PH_GER'; 
    fileprefix = '\PH_GER';
    Nvintages = 157;
elseif strcmp(country, 'United States')
    foldername = 'C:\Users\Philipp\Documents\Dissertation\sparse nowcasting\estim\PH_US';
    fileprefix = '\PH_US' ; 
    Nvintages = 229 ; 
end

counter = 0 ;
priorspec =  1:5 ; % Normal prior, PMNM prior 
rec_rolling = {'rec', 'rolling'} ;
surveys = {'level', 'diff'} ;
Nrs = 1:10 ; % # of factors
Nps = [1 3] ; % lag length in factor VAR
 
for v = 1 : Nvintages
    for prior = priorspec 
        for r = Nrs 
            for p = Nps
                for i_recroll = 1 : length(rec_rolling) 
                    flag_rolling = rec_rolling{i_recroll} ;
                    for i_survey = 1 : length(surveys) 
                        flag_survey = surveys{i_survey} ; 
                        filename = [foldername fileprefix '_v' num2str(v) '_prior' num2str(prior) '_Nr' num2str(r) '_Np' num2str(p) '_' flag_rolling '_' flag_survey  '.mat'] ;  
                        if exist(filename, 'file') == 2
                            counter = counter + 1 ;
                        else
                            disp(filename)
                        end
                    end
                end
            end
        end
    end
end