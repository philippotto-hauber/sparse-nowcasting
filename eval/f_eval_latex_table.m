clear; close all; clc; 
% flags 
flag_survey = 'level';
flag_sample = 'rec';
flag_truegdp = 'first';
Np = 1;
flag_country = 'GER';

% dir in & out
dir_load = ['C:\Users\Philipp\Documents\Dissertation\sparse nowcasting\eval\' flag_country '\' flag_survey ' ' flag_sample '\Np = ' num2str(Np) '\' flag_truegdp '\'] ;
dir_save = 'C:\Users\Philipp\Documents\Dissertation\sparse nowcasting\eval\latex_tables\' ; 
if exist(dir_save, 'dir') ~= 7;mkdir(dir_save); end  

% evaloptions (overwrite and amend)
evaloptions = load_evaloptions(flag_country);

evaloptions.Nhs = [2, 4] ; 

evaloptions.indexstarts = evaloptions.indexstarts([1 3]);
evaloptions.indexends = evaloptions.indexends([1 3]);
evaloptions.names_subsamples = {'full sample', 'post-crisis sample'};

evaloptions.metrics = {'sfe', 'logscore', 'crps'}; 
evaloptions.names_metrics = {'RMSFE', 'logS', 'CRPS'}; 

evaloptions.names_priors = {'NIG', 'MG', 'PMNM', 'HS', 'Nd'};

% load eval structure
load([dir_load 'results_eval.mat'])


% table info
Ncol_prior = 1; 
Ncol_model_benchmark = 1;
Ncols = length(evaloptions.names_subsamples) * length(evaloptions.Nhs) * length(evaloptions.metrics) + Ncol_prior + Ncol_model_benchmark;


% open file
filename = 'test.tex';
fid = fopen([dir_save filename], 'w'); 

% upper body
% caption
% label

% benchmark
vals_bar = getvals(results_eval, [], [], 0, 1, evaloptions);
fprintf(fid, ['%s & & %.2f' repmat('& %.2f ', 1, length(vals_bar)-1), ' %.2f'], 'B-AR', vals_bar);
fprintf(fid, '\\\\\n');

% loop over R's
for r = evaloptions.Nrs
    for p = evaloptions.Npriorspecs
        vals = getvals(results_eval, p, r, 0, 0, evaloptions);
        if p == evaloptions.Npriorspecs(1) 
            % print R=r in first column            
            fprintf(fid, ['%s & %s & %.2f ' repmat('& %.2f ', 1, length(vals_bar)-1), ' %.2f'], ['R=' num2str(r)], evaloptions.names_priors{p}, vals ./ vals_bar);
            fprintf(fid, '\\\\\n');
        else
            fprintf(fid, [' & %s & %.2f ' repmat('& %.2f ', 1, length(vals_bar)-1), ' %.2f'], evaloptions.names_priors{p}, vals ./ vals_bar);
            fprintf(fid, '\\\\\n');
        end
    end
end

% pool (only loop over priors!)
for p = evaloptions.Npriorspecs
    vals = getvals(results_eval, p, [], 1, 0, evaloptions);
    if p == evaloptions.Npriorspecs(1) 
        % print R=r in first column            
        fprintf(fid, ['%s & %s & %.2f ' repmat('& %.2f ', 1, length(vals_bar)-1), ' %.2f'], 'pool', evaloptions.names_priors{p}, vals ./ vals_bar);
        fprintf(fid, '\\\\\n');
    else
        fprintf(fid, [' & %s & %.2f ' repmat('& %.2f ', 1, length(vals_bar)-1), ' %.2f'], evaloptions.names_priors{p}, vals ./ vals_bar);
        fprintf(fid, '\\\\\n');
    end
end

% append lower body


% close file
fclose(fid);

function vals = getvals(results_eval, p, r, flag_pool, flag_benchmark, evaloptions)
% function that loops over samples, forecast metrics and horizons and
% returns the values
    counter = 0; 
    for s = 1 : length(evaloptions.indexstarts)
        indexstart = evaloptions.indexstarts(s); 
        indexend = evaloptions.indexends(s); 
        for ind_m = 1 : length(evaloptions.metrics)
            m = evaloptions.metrics{ind_m};
            for h = evaloptions.Nhs
                counter = counter + 1;
                if flag_benchmark == 1
                    tmp = results_eval.benchmark_BAR.horizon(h).(m)(indexstart:indexend);
                elseif flag_pool == 1
                    tmp = results_eval.priors(p).pool.horizon(h).(m)(indexstart:indexend);
                else
                    tmp = results_eval.priors(p).R(r).horizon(h).(m)(indexstart:indexend);
                end
                if strcmp(m, 'sfe')
                    vals(counter) = sqrt(mean(tmp));
                else
                    vals(counter) = mean(tmp);
                end
            end
        end
    end
end



                
                
                
            











