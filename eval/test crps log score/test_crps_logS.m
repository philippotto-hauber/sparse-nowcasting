clear; close all; 

addpath('../../eval/') % kde.m
h = 2;
R = 5;
dir_in = '../../../../../Desktop/results_eval mat files/';
flag_country = 'GER';

load([dir_in 'results_eval_' flag_country '_rec_level_Np1_first.mat'])

ind_q = 12;
dens = results_eval.priors(1).R(R).horizon(h).dens{ind_q};



% back out true gdp from sfe = (med_dens - truegdp)^2
sfe_dens = results_eval.priors(1).R(R).horizon(h).sfe(ind_q);
med_dens = median(dens);
truegdp_backout = med_dens - sqrt(sfe_dens);
load(['C:\Users\Philipp\Documents\Dissertation\sparse nowcasting\data\out\truegdp' flag_country '.mat'])
truegdp = truegdp_strct.first(ind_q);

crps = f_computeCRPS(dens, truegdp);
check_crps = results_eval.priors(1).R(R).horizon(h).crps(ind_q);

logS = f_computelogscore(dens, truegdp, 'kde');
check_logS = results_eval.priors(1).R(R).horizon(h).logscore(ind_q);

% export dens 
writematrix([truegdp, dens], 'data_test_crps_logS.csv')

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