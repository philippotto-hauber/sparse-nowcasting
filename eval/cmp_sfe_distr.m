clear; clc; close all; 

dir_out = 'C:\Users\Philipp\Documents\Dissertation\sparse nowcasting\eval\sfe_distr\';
if exist(dir_out, 'dir') ~= 7;mkdir(dir_out); end  
dir_in = 'C:\Users\Philipp\Documents\Dissertation\sparse nowcasting\documentation\results_eval mat files\';
dir_gdp = 'C:\Users\Philipp\Documents\Dissertation\sparse nowcasting\data\out\';

Ncountries = {'GER', 'US'};
models = {'BAR' 'NIG', 'PMNM', 'MG', 'HS'};
Nprior = [' ', 1, 3, 2, 4];

Ndraws = 1000; 
Np = 1; 
flag_survey = 'level';
flag_sample = 'rec';
flag_gdp = 'first';

Nhs = 2;
Nrs = [2, 5, 8]; 
for i_country = 1 : length(Ncountries)
    flag_country = Ncountries{i_country}; 
    
    % load true gdp
    load([dir_gdp 'truegdp' flag_country '.mat'])
    Nquarters = length(truegdp_strct.first);
        
    % load results eval mat file squared forecast errors
    load([dir_in 'results_eval_' flag_country '_' flag_sample '_' flag_survey '_Np' num2str(Np) '_' flag_gdp '.mat'])
    
    for Nr = Nrs
        for Nh = Nhs
            sfe_vec = NaN(Nquarters * Ndraws, length(models)); 


            for m = 1 : length(models)

                if strcmp(models{m}, 'BAR')
                    dens = results_eval.benchmark_BAR.horizon.dens;
                else
                    tmp1 = results_eval.priors(Nprior(m)).R(Nr).horizon(Nh).dens;
                    dens = [];
                    for n = 1 : Nquarters
                        dens = [dens; tmp1{n}]; 
                    end        
                end

                sfe = (dens - truegdp_strct.first') .^ 2;
                sfe_vec(:, m) = sfe(:); 
            end
            filename = [dir_out 'sfedistr_' flag_country '_' num2str(Nr) '_' num2str(Nh) '.csv'];
            fid = fopen(filename, 'w'); 
            fprintf(fid, ['%s', repmat(', %s', 1, length(models)-1)], models{:});
            fprintf(fid, '\n');
            for n_row = 1:size(sfe_vec, 1)            
                fprintf(fid, ['%.4f', repmat(', %.4f', 1, length(models)-1)], sfe_vec(n_row, :));   
                fprintf(fid, '\n');
            end
            fclose(fid);
        end
    end
end