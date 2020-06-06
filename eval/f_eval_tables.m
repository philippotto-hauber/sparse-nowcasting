function f_eval_tables(flag_survey, flag_sample, flag_truegdp, Np, flag_country) 

% - options ---
% ----------------------------
evaloptions = load_evaloptions(flag_country) ; 

% - directories ---
% ----------------------------
dir_save = [pwd '/' flag_country '/' flag_survey ' ' flag_sample '/Np = ' num2str(Np) '/' flag_truegdp '/tables/'] ; 
if exist(dir_save, 'dir') ~= 7;mkdir(dir_save); end  

% - load forecast structure
% -----------------
load(['results_eval mat files/results_eval_' flag_country '_' flag_sample '_' flag_survey '_Np' num2str(Np) '_' flag_truegdp '.mat'])

% - loop over subsamples
% -----------------
datatable_all = [] ;

for i = 1 : length(evaloptions.indexstarts)

    % get indicies
    indexstart = evaloptions.indexstarts(i) ; 
    indexend = evaloptions.indexends(i) ; 
    
    % clear matrices
    rmsfe_all = [] ; 
    logscore_all = [] ; 
    crps_all = [] ; 
    
    % loop over number of factors
    for index_r = 1:length(evaloptions.Nrs) 
        
        r = evaloptions.Nrs(index_r) ; 
        
        % - empty matrix to store squared forecast errors, log score and crps
        % -----------------------------------------
        
        rmsfe = [] ; 
        logscore = [] ; 
        crps = [] ; 

        for h = 1:evaloptions.Nhs 
            
            temp_rmsfe = [ sqrt(mean(results_eval.benchmark_BAR.horizon(h).sfe(indexstart:indexend))) ] ; 
            temp_logscore = [  mean(results_eval.benchmark_BAR.horizon(h).logscore(indexstart:indexend))] ;
            temp_crps = [ mean(results_eval.benchmark_BAR.horizon(h).crps(indexstart:indexend))] ; 

            for p = evaloptions.Npriorspecs   

                    temp_rmsfe = [temp_rmsfe; sqrt(mean(results_eval.priors(p).R(r).horizon(h).sfe(indexstart:indexend)))] ;
                    temp_logscore = [temp_logscore; mean(results_eval.priors(p).R(r).horizon(h).logscore(indexstart:indexend))] ;
                    temp_crps = [temp_crps; mean(results_eval.priors(p).R(r).horizon(h).crps(indexstart:indexend))] ;

            end
            rmsfe = [rmsfe temp_rmsfe] ; 
            logscore = [logscore temp_logscore] ; 
            crps = [crps temp_crps] ; 
        end
        
        % - compute metrics relative to benchmark
        % --------------------------------------------
        
        relrmsfe = [rmsfe(1,:); rmsfe(2:end,:)./rmsfe(1,:)] ;
        rellogscore = [logscore(1,:); logscore(2:end,:)./logscore(1,:)] ;
        relcrps = [crps(1,:); crps(2:end,:)./crps(1,:)] ;
        
        % - append current r to master arrays
        % --------------------------------------------        
        
        if r == 1            
            rmsfe_all = [rmsfe_all ; relrmsfe ] ; 
            logscore_all = [logscore_all ; rellogscore ] ; 
            crps_all = [crps_all ; relcrps ] ; 
        else
            rmsfe_all = [rmsfe_all ; relrmsfe(end-length(evaloptions.Npriorspecs)+1:end,:) ] ; 
            logscore_all = [logscore_all ; rellogscore(end-length(evaloptions.Npriorspecs)+1:end,:)  ] ; 
            crps_all = [crps_all ; relcrps(end-length(evaloptions.Npriorspecs)+1:end,:)  ] ; 
        end

    end
    
    % - equal weight pool
    % -----------------------------------------
    
    rmsfe_pool = [] ; 
    logscore_pool = [] ; 
    crps_pool = [] ; 

    for h = 1:evaloptions.Nhs 

        temp_rmsfe = [ sqrt(mean(results_eval.benchmark_BAR.horizon(h).sfe(indexstart:indexend))) ] ; 
        temp_logscore = [ mean(results_eval.benchmark_BAR.horizon(h).logscore(indexstart:indexend))] ;
        temp_crps = [ mean(results_eval.benchmark_BAR.horizon(h).crps(indexstart:indexend))] ; 


        for p = evaloptions.Npriorspecs   

                temp_rmsfe = [temp_rmsfe; sqrt(mean(results_eval.priors(p).pool.horizon(h).sfe(indexstart:indexend)))] ;
                temp_logscore = [temp_logscore; mean(results_eval.priors(p).pool.horizon(h).logscore(indexstart:indexend))] ;
                temp_crps = [temp_crps; mean(results_eval.priors(p).pool.horizon(h).crps(indexstart:indexend))] ;

        end
        rmsfe_pool = [rmsfe_pool temp_rmsfe] ; 
        logscore_pool = [logscore_pool temp_logscore] ; 
        crps_pool = [crps_pool temp_crps] ; 
    end
        
    % - compute metrics relative to benchmark
    % --------------------------------------------

    relrmsfe_pool = [rmsfe_pool(1,:); rmsfe_pool(2:end,:)./rmsfe_pool(1,:)] ;
    rellogscore_pool = [logscore_pool(1,:); logscore_pool(2:end,:)./logscore_pool(1,:)] ;
    relcrps_pool = [crps_pool(1,:); crps_pool(2:end,:)./crps_pool(1,:)] ;
    
    % - append pool to master arrays
    % --------------------------------------------             
         
    rmsfe_all = [rmsfe_all ; relrmsfe_pool(end-length(evaloptions.Npriorspecs)+1:end,:) ] ; 
    logscore_all = [logscore_all ; rellogscore_pool(end-length(evaloptions.Npriorspecs)+1:end,:)  ] ; 
    crps_all = [crps_all ; relcrps_pool(end-length(evaloptions.Npriorspecs)+1:end,:)  ] ; 
        
    % - save big table to xls
    % ------------------------------
    datatable = [rmsfe_all logscore_all crps_all] ; 
    csvwrite([dir_save 'datatable_allR_plus_pool_' results_eval.quarters{indexstart} '_' results_eval.quarters{indexend} '.xlsx'],datatable)    
    
    % - save append current subsample to big big table
    % ----------------------------------------------------
    datatable_all = [ datatable_all datatable ] ; 
end

% - save big big table to xls
% ------------------------------
csvwrite([dir_save 'datatable_all_R_plus_pool_all_samples.xlsx'],datatable_all)
