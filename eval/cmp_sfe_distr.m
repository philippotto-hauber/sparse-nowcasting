clear; clc; close all; 
load('C:\Users\Philipp\Documents\Dissertation\sparse nowcasting\documentation\results_eval mat files\results_eval_GER_rec_level_Np1_first.mat')
load('C:\Users\Philipp\Documents\Dissertation\sparse nowcasting\data\out\truegdpGER.mat')

models = {'NIG', 'PMNM', 'MG', 'HS'};
Nprior = [1, 3, 2, 4];
Nquarters = 52; 
Ndraws = 1000; 

Nhs = 1:4;
Nrs = [1, 2, 5, 8]; 
counter = 1; 
for Nr = Nrs
    for Nh = Nhs
        sfe_vec = NaN(Nquarters * Ndraws, length(models)); 


        for m = 1 : length(models)

            if strcmp(models{m}, 'BAR')
                dens = results_eval.benchmark_BAR.horizon.dens;
            else
                tmp1 = results_eval.priors(Nprior(m)).R(Nr).horizon(Nh).dens;
                dens = [];
                for n = 1 : 52
                    dens = [dens; tmp1{n}]; 
                end        
            end

            sfe = (dens - truegdp_strct.first') .^ 2;
            sfe_vec(:, m) = sfe(:); 
        end

        subplot(length(Nrs), length(Nhs), counter)
        boxplot(sfe_vec);xticklabels(models);title(['Nr = ' num2str(Nr) ', Nh = ' num2str(Nh)])
        counter = counter + 1; 
    end
end