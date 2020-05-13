function PH_GER( Nspec )
if (isdeployed)
    maxNumCompThreads(1);
    Nspec = str2double(Nspec) ;
end

% ----------------------------------------------------------------------- %
% ----- options --------------------------------------------------------- %
% ----------------------------------------------------------------------- %

options.flag_samplemoments = 0 ;
dir_out = '' ; % if not empty, needs to include \ at the end  

Nvintages = 1:157 ; 
Nrs = 1:10 ; 
Npriors = 1:5 ; 
Nps = [1 3] ;
Nmod = 1:4 ;

% ----------------------------------------------------------------------- %
% ----- model spec ------------------------------------------------------ %
% ----------------------------------------------------------------------- %

% calculate combinations of specifications
Nrs_temp = repmat(Nrs, 1, length(Npriors) * length(Nvintages) * length(Nps) * length(Nmod)) ; 
Npriors_temp = repmat(kron(Npriors, ones(1, length(Nrs))), 1 , length(Nvintages) * length(Nps) * length(Nmod)) ;
Nvintages_temp = repmat(kron(Nvintages, ones(1, length(Nrs) * length(Npriors))), 1, length(Nps) * length(Nmod)) ; 
Nps_temp = repmat(kron(Nps, ones(1, length(Nrs) * length(Npriors) * length(Nvintages))), 1, length(Nmod)) ;
Nmod_temp = kron(Nmod, ones(1, length(Nrs) * length(Npriors) * length(Nvintages) * length(Nps))) ;

% back out current one
modspec = Nmod_temp(Nspec) ; 
v = Nvintages_temp(Nspec) ; 
options.Nr = Nrs_temp(Nspec) ; % # of static factors
options.Ns = options.Nr ; % # of static factors
options.Np = Nps_temp(Nspec) ; % # of lags in factor VAR
options.Nj = 0 ; % # of lags in eps
options.priorswitch = Npriors_temp(Nspec) ;

% model flag
if modspec == 1
    % surveys in levels, recursive estimation window
    flag_surveydiff = 'level' ;
    flag_rolling = 'rec' ; 
elseif modspec == 2 
    % surveys in 1st diffs, recursive estimation window
    flag_surveydiff = 'diff' ;
    flag_rolling = 'rec' ;
elseif modspec == 3
    % surveys in levels, rolling estimation window
    flag_surveydiff = 'level' ;
    flag_rolling = 'rolling' ;
elseif modspec == 4
    % surveys in 1st diffs, rolling estimation window
    flag_surveydiff = 'diff' ;
    flag_rolling = 'rolling' ;
end

% ----------------------------------------------------------------------- %
% ----- load data ------------------------------------------------------- %
% ----------------------------------------------------------------------- %

load('datasetsGER.mat')
options.vintagedate = datasets.vintage(v).vintage ; 
options.modspec = modspec ;

if strcmp(flag_surveydiff, 'diff')
    X = [ datasets.vintage(v).data_ifo_diff , ...
          datasets.vintage(v).data_bubartd( 2 : end , : ) , ...
          datasets.vintage(v).data_financial( 2 : end , :) ]' ; 
    yQ = datasets.vintage(v).data_gdp( 2 : end , 1 )' ;
elseif strcmp(flag_surveydiff, 'level')
    X = [ datasets.vintage(v).data_ifo_level , ...
          datasets.vintage(v).data_bubartd , ...
          datasets.vintage(v).data_financial(:,:)]' ; 
    yQ = datasets.vintage(v).data_gdp' ;
else 
    disp('incorrectly specified flag_surveydiff. Abort!');
    abort;
end

index_use = sum(isnan(X),2)<(0.5 * size(X,2) ) ; 

X = X(index_use,:) ;  

% remove outliers from X
[X_xoutl,n_outl]=remove_outliers(X') ;
X = X_xoutl' ; 

% standardize
Xstand = ( X - nanmean( X , 2 ) ) ./ nanstd( X , [] , 2 ) ;
X = Xstand ; 
mean_gdp = nanmean( yQ ) ; 
std_gdp = nanstd( yQ ) ; 
yQ = ( yQ - mean_gdp ) / std_gdp ; 


% ------------------
% - recursive or rolling estimation window ?

if strcmp(flag_rolling, 'rolling') 
    yQ = yQ(1,end-120:end) ;
    X = X(:,end-120:end) ;
elseif strcmp(flag_rolling, 'rec')
    % do nothing
else
    disp('incorrectly specified flag_rolling. Abort!');
    abort;
end

% --------------
% options
options.Nm = size(X,1);
options.Nq = size(yQ,1);
options.Nn = options.Nm + options.Nq ; 
options.Nt = size(X,2);
vintagemonth = month(datasets.vintage(v).vintage,'yyyy-mm') ; 
if ismember(vintagemonth,[3 6 9 12])
    options.Nh = 3 ; 
elseif ismember(vintagemonth,[1 4 7 10])
    options.Nh = 2 ; 
elseif ismember(vintagemonth,[2 5 8 11])
    options.Nh = 1 ;
end

% ----------------------------------------------------------------------- %
% ----- Gibbs Sampler --------------------------------------------------- %
% ----------------------------------------------------------------------- %

% options
options.Nburnin = 5000 ; % # of burn-ins
options.Nreplic = 10000 ; % # of replics
options.Nthin = 10 ; % store each options.thinning-th draw
options.Ndisplay = 1000 ;  % display each options.display-th iteration

% priors
priors = loadpriors(options,options.priorswitch);

% call GibbsSampler.m
draws = GibbsSampler( X, yQ, priors , options );

% restandardize nowcast (and forecast)
if options.flag_samplemoments == 0 
    draws.nowcast = std_gdp * draws.nowcast + mean_gdp ; 
    if options.Nh == 3 
        draws.forecast = std_gdp * draws.forecast + mean_gdp ; 
    end
elseif options.flag_samplemoments == 1
    draws.nowcast_mean = std_gdp * draws.nowcast_mean + mean_gdp ; 
    draws.nowcast_var = std_gdp^2 * draws.nowcast_var ; 
    if options.Nh == 3 
        draws.forecast_mean = std_gdp * draws.forecast_mean + mean_gdp ; 
        draws.forecast_var = std_gdp^2 * draws.forecast_var ; 
    end
end

% ----------------------------------------------------------------------- %
% ----- save results to mat-file ---------------------------------------- %
% ----------------------------------------------------------------------- %
save([dir_out 'PH_GER_v' num2str(Nvintages_temp(Nspec))  '_prior' num2str(Npriors_temp(Nspec)) '_Nr' num2str(Nrs_temp(Nspec)) '_Np' num2str(Nps_temp(Nspec)) '_' flag_rolling '_' flag_surveydiff  '.mat'], 'draws');






