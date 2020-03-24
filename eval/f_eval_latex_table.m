clear; close all; clc; 
% flags 
flag_survey = 'level';
flag_sample = 'rec';
flag_truegdp = 'first';
Np = 1;
flag_country = 'GER';

if strcmp(flag_country, 'GER')
    name_country = 'Germany';
elseif strcmp(flag_country, 'US')
    name_country = 'United States';
end

% dir in & out
dir_load = ['C:\Users\Philipp\Documents\Dissertation\sparse nowcasting\eval\' flag_country '\' flag_survey ' ' flag_sample '\Np = ' num2str(Np) '\' flag_truegdp '\'] ;
dir_save = 'C:\Users\Philipp\Documents\Dissertation\sparse nowcasting\eval\latex_tables\' ; 
if exist(dir_save, 'dir') ~= 7;mkdir(dir_save); end  

% font size
fontsize = 'tiny'; 

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
Ncols_per_sample = length(evaloptions.names_subsamples) * length(evaloptions.Nhs) * length(evaloptions.metrics);


% open file
filename = ['table_' flag_country '_' flag_truegdp '_' flag_survey '_' flag_sample '_Np' num2str(Np) '.tex'];
fid = fopen([dir_save filename], 'w'); 

% upper body
str_label = 'label1';
str_caption = 'captionA';
fprintf(fid, ['\\begin{threeparttable}[p]\n\\caption{' str_caption '}\n\\label{' str_label '}\n\\' fontsize '\n']);
fprintf(fid, '\\begin{tabular}{ c l ');
fprintf(fid, [repmat('c ', 1, Ncols_per_sample) '}\n\\toprule\n']);

counter = 1;
str_tmp = ['\\multicolumn{' num2str(Ncols_per_sample/length(evaloptions.names_subsamples)) '}{c}{' evaloptions.names_subsamples{counter} '}'];
while counter < length(evaloptions.names_subsamples)
    counter = counter + 1;
    str_tmp = [str_tmp [' & \\multicolumn{' num2str(Ncols_per_sample/length(evaloptions.names_subsamples)) '}{c}{' evaloptions.names_subsamples{counter} '}']];
end
fprintf(fid, [' & & ' str_tmp '\\\\\n']);

counter = 1;
str_tmp = ['\\multicolumn{' num2str(length(evaloptions.Nhs)) '}{c}{' evaloptions.names_metrics{counter} '}'];
while counter < length(evaloptions.names_metrics)
    counter = counter + 1;
    str_tmp = [str_tmp [' & \\multicolumn{' num2str(length(evaloptions.Nhs)) '}{c}{' evaloptions.names_metrics{counter} '}']];
end
fprintf(fid, [' & & ' [repmat(str_tmp, 1, length(evaloptions.names_subsamples)-1) ' & ' str_tmp] '\\\\\n']);

counter = 1;
str_tmp = results_eval.priors(1).R(1).horizon(evaloptions.Nhs(counter)).name;
while counter < length(evaloptions.Nhs)
    counter = counter + 1;
    str_tmp = [str_tmp [' & ' results_eval.priors(1).R(1).horizon(evaloptions.Nhs(counter)).name]];
end
fprintf(fid, [' & & ' repmat([str_tmp ' & '], 1, length(evaloptions.metrics) * length(evaloptions.names_subsamples)-1) str_tmp '\\\\\n']);

fprintf(fid, '\\midrule\n');

% benchmark
vals_bar = getvals(results_eval, [], [], 0, 1, evaloptions);
addline(fid, vals_bar, 'B-AR', [])
fprintf(fid, '\\hspace{0.1cm}\\\\\n');

% loop over R's
for r = evaloptions.Nrs
    for p = evaloptions.Npriorspecs
        vals = getvals(results_eval, p, r, 0, 0, evaloptions);
        if p == evaloptions.Npriorspecs(1) 
            % print R=r in first column      
            addline(fid, vals ./ vals_bar, ['R=' num2str(r)], evaloptions.names_priors{p})
        else
            addline(fid, vals ./ vals_bar, '', evaloptions.names_priors{p})
        end
    end
    fprintf(fid, '\\hspace{0.1cm}\\\\\n');
end

% pool (only loop over priors!)
for p = evaloptions.Npriorspecs
    vals = getvals(results_eval, p, [], 1, 0, evaloptions);
    if p == evaloptions.Npriorspecs(1) 
        % print R=r in first column      
        addline(fid, vals ./ vals_bar, 'pool', evaloptions.names_priors{p})
    else
        addline(fid, vals ./ vals_bar, '', evaloptions.names_priors{p})
    end
end

% append lower body
fprintf(fid, '\\bottomrule\n\\end{tabular}\n\\begin{tablenotes}\n\\' fontsize '\n\\item ');

% notes under tables (split across multiple fprintf's for convenience
fprintf(fid, 'RMSFE, logS and CRPS for the models are relativ to the B-AR benchmark (see text for details)');
fprintf(fid, 'The log score is negatively orientated so that a value in the table below 1 corresponds to a better performance than the benchmark.');
fprintf(fid, 'Forecast horizon h is in months. The full sample period is 2000Q1-2018Q4, the post-crisis sample starts in 2010Q1 and ends in 2018Q4.'); 
fprintf(fid, '\n\\end{tablenotes}\n\\end{threeparttable}\n');

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

function addline(fid, vals, r_name, p_name)

% add leading columns and values
fprintf(fid, ['%s & %s & %.2f' repmat(' & %.2f ', 1, length(vals)-1), ' %.2f'], r_name, p_name, vals);

% add backslashes and new line
fprintf(fid, '\\\\\n');

end




                
                
                
            











