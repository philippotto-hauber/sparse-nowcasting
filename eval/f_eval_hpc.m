function f_eval_hpc(flag_survey, flag_sample, flag_truegdp, Np, flag_country)

% - options ---
% ----------------------------

evaloptions = load_evaloptions(flag_country) ; 

% - directories ---
% ----------------------------
dir_truegdp = '' ;
%dir_truegdp = 'C:\Users\Philipp\Documents\Dissertation\sparse nowcasting\data\out\';
dir_models = ['../PH_' flag_country '/matfiles/'] ; 
%dir_models = ['../../../Dissertation/sparse nowcasting/documentation/PH_' flag_country '/'];
   
dir_benchmark = ['benchmark_' flag_country '/'] ; 
%dir_benchmark = ['../../../Dissertation/sparse nowcasting/documentation/benchmark_' flag_country '/'] ; 
dir_out = ['results_eval mat files/'] ; 
if exist(dir_out, 'dir') ~= 7;mkdir(dir_out); end  

% - load true gdp mat-file ---
% ----------------------------
load([dir_truegdp 'truegdp' flag_country '.mat'])

% - start looping
% --------------------------

for p = evaloptions.Npriorspecs
    % - store prior names
    % ------------------------
    switch p
        case 1
            priorname = 'Normal-Inverse Gamma' ; 
        case 2
            priorname = 'Multiplicative Gamma' ; 
        case 3
            priorname = 'Point mass Normal mixture' ;
        case 4
            priorname = 'Horseshoe+' ;
        case 5
            priorname = 'Normal-diffuse' ;
    end
    results_eval.priors(p).name = priorname ; 

    for q = 1 : evaloptions.Nquarters   
        if strcmp(flag_truegdp,'first')
            truegdp = evaloptions.multfac*truegdp_strct.first(q) ;
        elseif strcmp(flag_truegdp,'second')
            truegdp = evaloptions.multfac*truegdp_strct.second(q) ;
        elseif strcmp(flag_truegdp,'final')
            truegdp = evaloptions.multfac*truegdp_strct.final(q) ;
        end

        results_eval.quarters{q} = truegdp_strct.quarters{q} ;

        % - get vintage indices
        % ------------------------
        [index_vs, flag_models_vs, flag_BAR_vs] = f_mapping_q_to_v(q,evaloptions.Nhs,flag_country) ; 
        if evaloptions.Nhs ~= length(index_vs)    
            disp('Number of vintages per quarter do not match. Abort execution')
            return
        end

        for h = 1:length(index_vs)
            dens_pool = cell( 1 , length( evaloptions.Nrs ) ) ; 
            for index_r = 1 : length( evaloptions.Nrs ) + 1 % factors + equal weight pool

                % - create subfields
                % ------------------------
                if index_r == length( evaloptions.Nrs ) + 1 
                    results_eval.priors(p).pool.horizon(1).name = 'h=3' ;
                    results_eval.priors(p).pool.horizon(2).name = 'h=2' ;
                    results_eval.priors(p).pool.horizon(3).name = 'h=1' ;
                    results_eval.priors(p).pool.horizon(4).name = 'h=0' ; 
                else
                    r = evaloptions.Nrs( index_r ) ;
                    results_eval.priors(p).R(r).horizon(1).name = 'h=3' ;
                    results_eval.priors(p).R(r).horizon(2).name = 'h=2' ;
                    results_eval.priors(p).R(r).horizon(3).name = 'h=1' ;
                    results_eval.priors(p).R(r).horizon(4).name = 'h=0' ; 
                end
                 if r == 1

                    % ------------------------------------------------------------------------ %
                    % - benchmark: BAR(-1)
                    % ------------------------------------------------------------------------ %

                    % - load results mat-file
                    % -------------------------
                    if strcmp(flag_sample, 'rolling')
                                load([dir_benchmark 'PH_' flag_country '_v' num2str(index_vs(h)) '_roll.mat'])
                    elseif strcmp(flag_sample, 'rec')
                                load([dir_benchmark 'PH_' flag_country '_v' num2str(index_vs(h)) '_' flag_sample '.mat'])
                    end

                    % - select correct row and multiply with 100

                    % -------------------------
                    if flag_BAR_vs(h) == 1 % forecast
                        draws_temp = evaloptions.multfac * draws.forecast( evaloptions.Nthin : evaloptions.Nthin : end ) ;
                    else
                        draws_temp = evaloptions.multfac * draws.nowcast( evaloptions.Nthin : evaloptions.Nthin : end ) ;
                    end

                    % - store density
                    % -------------------------
                    results_eval.benchmark_BAR.horizon(h).dens(q,:) = draws_temp ;
                    
                    % - store squared forecast errors, log score and crps
                    % -------------------------
                    results_eval.benchmark_BAR.horizon(h).sfe(q,:) = f_computesfe(draws_temp, truegdp, evaloptions.computesfe);;
                    results_eval.benchmark_BAR.horizon(h).logscore(q,1) = f_computelogscore(draws_temp, truegdp, evaloptions.computelogscore);
                    results_eval.benchmark_BAR.horizon(h).crps(q,1) = f_computeCRPS(draws_temp', truegdp); % crps ;
                end

                % ------------------------------------------------------------------------ %
                % - models
                % ------------------------------------------------------------------------ %

                if index_r <= length( evaloptions.Nrs )
                    % - load results mat-file
                    % -------------------------
                    load([dir_models 'PH_' flag_country '_v' num2str(index_vs(h)) '_prior' num2str(p) '_Nr' num2str(r) '_Np' num2str(Np) '_' flag_sample '_' flag_survey '.mat'])

                    % - select correct row and multiply with multiplication factor
                    % -------------------------
                    if flag_models_vs(h) == 1 % forecast
                        draws_temp = evaloptions.multfac * draws.forecast( evaloptions.Nthin : evaloptions.Nthin : end ) ;
                    else
                        draws_temp = evaloptions.multfac * draws.nowcast( evaloptions.Nthin : evaloptions.Nthin : end ) ;
                    end     

                    % - store density
                    % -------------------------
                    results_eval.priors(p).R(r).horizon(h).dens{q} = draws_temp ;  

                    % - store for pool!
                    % -------------------------
                    dens_pool{r} = draws_temp ; 

                    % - compute squared forecast errors, log score and
                    % -----------------------------
                    results_eval.priors(p).R(r).horizon(h).sfe(q,:) = f_computesfe(draws_temp, truegdp, evaloptions.computesfe);
                    results_eval.priors(p).R(r).horizon(h).logscore(q,1) = f_computelogscore(draws_temp, truegdp, evaloptions.computelogscore);
                    results_eval.priors(p).R(r).horizon(h).crps(q,1)= f_computeCRPS(draws_temp', truegdp);
                else
                    % equal weight pool
                    results_eval.priors(p).pool.horizon(h).dens{q} = f_pooldens_eqwgts(dens_pool , evaloptions.Nmultpool*evaloptions.Ndraws ) ;
                    results_eval.priors(p).pool.horizon(h).sfe(q,:) = f_computesfe(results_eval.priors(p).pool.horizon(h).dens{q}, truegdp, evaloptions.computesfe);
                    results_eval.priors(p).pool.horizon(h).logscore(q,1) = f_computelogscore(results_eval.priors(p).pool.horizon(h).dens{q}, truegdp, evaloptions.computelogscore);
                    results_eval.priors(p).pool.horizon(h).crps(q,1) = f_computeCRPS(results_eval.priors(p).pool.horizon(h).dens{q},truegdp);
                end                    
            end                        
        end                    
    end
end
    
    % - save results to mat-file
    % ----------------------------- 
    save([dir_out 'results_eval_' flag_country '_' flag_sample '_' flag_survey '_Np' num2str(Np) '_' flag_truegdp '.mat'],'results_eval')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [index_vs, flag_models, flag_BAR] = f_mapping_q_to_v(q,Nh,flag_country)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% This code maps the quarter q to the relevant vintages that produce
%%%% now-/forecasts for that quarter. It returns an index of those 
%%%% vintages - index_vs - as well as a flag of length(index_vs)
%%%% highlighting the relation of that vintage to the nowcasting quarter.
%%%% This is relevant for extracting the correct densities from the 
%%%% saved output mat-files PH12_x_x_x and determining the forecast horizon
%%%% in the naive B-AR model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(flag_country,'US')    
index_vs = (q*3-2) + [0:Nh-1] ; 
flag_models = zeros(1,length(index_vs)) ; 
flag_models(1) = 1 ;
flag_BAR = zeros(1,length(index_vs)) ;
flag_BAR(1) = 1 ;
elseif strcmp(flag_country,'GER')
index_vs = (q*3-2) + [0:Nh-1] ; 
flag_models = zeros(1,length(index_vs)) ; 
flag_models(1) = 1 ;
flag_BAR = zeros(1,length(index_vs)) ;
flag_BAR(1) = 1 ;  % flag_BAR and flag_models should be the same!!!!!
end
end

function denspool = f_pooldens_eqwgts(dens, Nreplics)

denspool = NaN(1,Nreplics);
Nmodels = size( dens , 2 ) ; 
Ndraws = NaN( 1 , Nmodels ) ; 
for i = 1 : Nmodels
    Ndraws( i ) = length( dens{ i } ) ; 
end
minNdraws = min( Ndraws ) ; 

for m = 1 : Nreplics
    indic_dens = unidrnd( Nmodels );
    indic_draws = unidrnd( minNdraws );
    dens_temp = dens{ 1 , indic_dens } ;
    denspool(1,m) = dens_temp( indic_draws );
end
end

function sfe = f_computesfe(draws, truegdp, flag_computesfe)
    if strcmp(flag_computesfe, 'mean')
        sfe = (mean(draws) - truegdp)^2;
    else
        sfe = (draws - truegdp) .^ 2;
    end
end

function logscore = f_computelogscore(draws,truegdp,flag_computelogscore)

if strcmp( flag_computelogscore , 'kde')    
    
    % kernel density estimate following Chovez (2010) ????
    range_kde = max(draws) - min(draws);
    min_kde = min(min(draws), truegdp) - range_kde/10;
    max_kde = max(max(draws), truegdp) + range_kde/10; 
    [~,dens,xs,~]=kde(draws, 2^8, min_kde, max_kde) ; 
    index_truegdp = sum( xs <= truegdp ) ; 
    logscore = log( dens( index_truegdp ) ) ; 
    
elseif strcmp( flag_computelogscore , 'ksdensity')
    
    % built-in MATLAB function with default settings
    mindens = min( min( draws ) , truegdp ) ; 
    maxdens = max( max( draws ) , truegdp ) ; 
    rangedens = maxdens - mindens ; 
    pts = ( mindens - rangedens / 10 ) : rangedens/ 100 : ( maxdens + rangedens / 10 ) ;
    [dens,xs] = ksdensity(draws,pts) ;
    index_truegdp = sum( xs <= truegdp ) ; 
    logscore = log( dens( index_truegdp ) ) ; 
    
elseif strcmp( flag_computelogscore , 'normal approx')
    
    % assume predictive density is Normal
    logscore = log( normpdf( truegdp , mean( draws ) , std( draws ) ) ) ;
    
end

logscore = logscore * (-1) ; 
end

function crps = f_computeCRPS(X,y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This functions computes the continously ranked 
%%% probability score (CRPS) following the empirical
%%% CDF approach outlined in Krueger et al (2017)
%%% "Probabilistic Forecasting and Comparative
%%%  Model Assessement Based on Markov Chain 
%%%  Monte Carlo Output
%%% -----------------------------------------------------
%%% INPUTS: X - mx1 vector of draws from the predictive density
%%%         y - scalar realization
%%% OUTPUT: scalar CRPS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ------------------
% - equation 10
% ------------------

crps1 = sum(abs( X - y )) / length( X );

crps2 = 0 ;
for m = 1 : length(X)
    crps2 = crps2 + sum( abs( X - X(m) ) ) ;
end

crps = crps1 - 1/( 2 * length(X) ^ 2 ) * crps2;
end


