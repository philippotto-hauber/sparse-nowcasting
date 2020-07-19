function f_eval_latex_table(flag_survey, flag_sample, flag_truegdp, Np, flag_country)

    % - dir in & out ---------------------------------------------------- %
    % ------------------------------------------------------------------- %
    dir_save = 'latex_tables/' ; 
    if exist(dir_save, 'dir') ~= 7;mkdir(dir_save); end  

    % - load eval structure and options --------------------------------- %
    % ------------------------------------------------------------------- %
    load(['results_eval mat files/results_eval_' flag_country '_' flag_sample '_' flag_survey '_Np' num2str(Np) '_' flag_truegdp '.mat'])
    evaloptions = load_evaloptions(flag_country);

    % - USER INPUT ------------------------------------------------------ %
    % ------------------------------------------------------------------- %
    fontsize = 'scriptsize'; % 
    linestep = '0.0'; % step size between model blocks

    evaloptions.Nhs = [2] ; % forecast horizons
    evaloptions.Nrs = [1, 2, 5, 8] ; % # of factors to include in the table (pool always included and note that it is calculated over Nr = 1:10)

    % subsamples
    evaloptions.indexstarts = evaloptions.indexstarts([1 3]);
    evaloptions.indexends = evaloptions.indexends([1 3]);
    evaloptions.names_subsamples = {'full sample', 'post-crisis sample'};

    % forecast evaluation metrics
    evaloptions.metrics = {'sfe', 'logscore', 'crps'}; 
    evaloptions.names_metrics = {'RMSFE', 'logS', 'CRPS'}; 
    
    % prior names
    evaloptions.names_priors = {'NIG', 'MG', 'PMNM', 'HS', 'Nd'};

    % label
    str_label = ['table:' flag_country '_' flag_truegdp '_' flag_survey '_' flag_sample '_Np' num2str(Np)];
    
    % caption
    if strcmp(flag_country, 'GER'); name_country = 'Germany'; elseif strcmp(flag_country, 'US'); name_country = 'United States'; end
    str_caption = [name_country ' (' flag_truegdp ', ' flag_survey ', ' flag_sample ')']; 

    % notes to table
    if strcmp(flag_country, 'GER')
    str_notes = ['RMSFE is the root mean squared forecast error, logS and CRPS are the average log score and continuous ranked probability score. All entries for the factor models are relative to the B-AR benchmark (see text for details)' ...
                 ' and negatively orientated so that a value in the table below 1 corresponds to a better performance than the benchmark.' ...
                 ' The forecast horizon h is in months. The full sample period is 2006Q1-2018Q4, the post-crisis sample starts in 2010Q1 and ends in 2018Q4.'];
    elseif strcmp(flag_country, 'US')
    str_notes = ['RMSFE is the root mean squared forecast error, logS and CRPS are the average log score and continuous ranked probability score. All entries for the factor models are relative to the B-AR benchmark (see text for details)' ...
                 ' and negatively orientated so that a value in the table below 1 corresponds to a better performance than the benchmark.' ...
                 ' The forecast horizon h is in months. The full sample period is 2000Q1-2018Q4, the post-crisis sample starts in 2010Q1 and ends in 2018Q4.'];    
    end
    
    % - open file ------------------------------------------------------- %
    % ------------------------------------------------------------------- %
    filename = ['table_' flag_country '_' flag_truegdp '_' flag_survey '_' flag_sample '_Np' num2str(Np) '.tex'];
    fid = fopen([dir_save filename], 'w'); 

    % - upper body ------------------------------------------------------ %
    % ------------------------------------------------------------------- %
    
    Ncols_per_sample = length(evaloptions.names_subsamples) * length(evaloptions.Nhs) * length(evaloptions.metrics);

    fprintf(fid, ['\\begin{threeparttable}[p]\n\\caption{' str_caption '}\n\\label{' str_label '}\n\\' fontsize '\n']);
    fprintf(fid, '\\begin{tabularx}{0.7\\textwidth}{ c l');

    tmp1 = repmat(' Y ', 1, length(evaloptions.Nhs));
    counter = 1;
    tmp2 = tmp1; 
    while counter < length(evaloptions.metrics)
        tmp2 = [tmp2 ' ' tmp1];
        counter = counter + 1;
    end
    
    counter = 1; 
    tmp3 = tmp2; 
    while counter < length(evaloptions.names_subsamples)        
        tmp3 = [tmp3 '  ' tmp2];
        counter = counter + 1;
    end
    
    fprintf(fid,[tmp3 '}\n\\toprule\n']);
    
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
    
    counter = 1;
    fprintf(fid, [' & & ' str_tmp]); 
    while counter < length(evaloptions.names_subsamples)
        fprintf(fid, [' & ' str_tmp]);
        counter = counter + 1;
    end
    fprintf(fid, '\\\\\n');

    counter = 1;
    str_tmp = results_eval.priors(1).R(1).horizon(evaloptions.Nhs(counter)).name;
    while counter < length(evaloptions.Nhs)
        counter = counter + 1;
        str_tmp = [str_tmp [' & ' results_eval.priors(1).R(1).horizon(evaloptions.Nhs(counter)).name]];
    end
    fprintf(fid, [' & & ' repmat([str_tmp ' & '], 1, length(evaloptions.metrics) * length(evaloptions.names_subsamples)-1) str_tmp '\\\\\n']);

    fprintf(fid, '\\midrule\n');
    
    % - benchmark ------------------------------------------------------- %
    % ------------------------------------------------------------------- %
    
    vals_bar = getvals(results_eval, [], [], 0, 1, evaloptions);
    addline(fid, vals_bar, 'B-AR', [])
    fprintf(fid, ['\\vspace{' linestep 'cm}\\\\\n']);

    % - models for different R ------------------------------------------ %
    % ------------------------------------------------------------------- %
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
        fprintf(fid, ['\\vspace{' linestep 'cm}\\\\\n']);
    end

    % - pool (only loop over priors!) ----------------------------------- %
    % ------------------------------------------------------------------- %
    for p = evaloptions.Npriorspecs
        vals = getvals(results_eval, p, [], 1, 0, evaloptions);
        if p == evaloptions.Npriorspecs(1) 
            % print R=r in first column      
            addline(fid, vals ./ vals_bar, 'pool', evaloptions.names_priors{p})
        else
            addline(fid, vals ./ vals_bar, '', evaloptions.names_priors{p})
        end
    end

    % - append lower body ----------------------------------------------- %
    % ------------------------------------------------------------------- %
    fprintf(fid, ['\\bottomrule\n\\end{tabularx}\n\\begin{tablenotes}\n\\' fontsize '\n\\item ' str_notes]);
    fprintf(fid, '\n\\end{tablenotes}\n\\end{threeparttable}\n');

    % - close file ------------------------------------------------------ %
    % ------------------------------------------------------------------- %
    fclose(fid);
end

% - FUNCTIONS ----------------------------------------------------------- %
% ----------------------------------------------------------------------- %

function vals = getvals(results_eval, p, r, flag_pool, flag_benchmark, evaloptions)
% function that loops over samples, forecast metrics and horizons and
% returns the values of RMSFE, log score and CRPS
    counter = 0; 
    for s = 1 : length(evaloptions.indexstarts)
        indexstart = evaloptions.indexstarts(s); 
        indexend = evaloptions.indexends(s); 
        for ind_m = 1 : length(evaloptions.metrics)
            m = evaloptions.metrics{ind_m};
            for h = evaloptions.Nhs
                counter = counter + 1;
                if flag_benchmark == 1
                    tmp = results_eval.benchmark_BAR.horizon(h).(m)(indexstart:indexend, :);
                elseif flag_pool == 1
                    tmp = results_eval.priors(p).pool.horizon(h).(m)(indexstart:indexend, :);
                else
                    tmp = results_eval.priors(p).R(r).horizon(h).(m)(indexstart:indexend, :);
                end
                if strcmp(m, 'sfe')
                    vals(counter) = mean(sqrt(mean(tmp, 1))); % RMSFE is calculated for each draw, then averaged!
                else
                    vals(counter) = mean(tmp);
                end
            end
        end
    end
end

function addline(fid, vals, r_name, p_name)

% r_name in bold!
if not(isempty(r_name))
    r_name = ['\textbf{' r_name '}'];
end

% add leading columns and values
fprintf(fid, ['%s & %s & %.2f' repmat(' & %.2f ', 1, length(vals)-1), ' %.2f'], r_name, p_name, vals);

% add backslashes and new line
fprintf(fid, '\\\\\n');

end




                
                
                
            











