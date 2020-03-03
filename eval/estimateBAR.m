clear;close all;clc;

% -------------------------------------------- %
% - preliminaries

dir_local = 'C:\Users\Philipp\Documents\Dissertation\sparse nowcasting\' ; 
list_countries = {'GER', 'US'} ; 

% -------------------------------------------- %
% - loop over countries

for i = 1 : length(list_countries)
    flag_country = list_countries{i} ; 

    if strcmp(flag_country, 'GER')
        Nvintages = 157 ; 
        filename_save = 'eval\GER\benchmark\' ; 
        load([dir_local 'data\out\datasetsGER.mat'])
    elseif strcmp(flag_country, 'US')
        Nvintages = 229 ; 
        filename_save = '\eval\US\benchmark\' ; 
        load([dir_local 'data\out\datasetsUS.mat'])
    end


    % - B-AR(1) stuff
    % --------------------------
    pi1 = 3;
    pi2 = 2;
    Np = 1;
    yQ = 0 ; 
    for flag_rolling = [0 1]
        if flag_rolling == 0 
            name_flag_rolling = 'rec' ;
        elseif flag_rolling == 1 
            name_flag_rolling = 'roll' ;
        end

        for v = 1 : Nvintages

            % - update yQold
            % -----------------------------
            yQold = yQ ; 

            % - get vintage month & data
            % -----------------------------
            if strcmp(flag_country, 'GER')
                vintage_month = month(datasets.vintage(v).vintage) ;         

                y = datasets.vintage(v).data_gdp ;    
                meangdp = nanmean( y ) ; 
                stdgdp = nanstd( y ) ; 
                yQ = y(~isnan(y)) ; 
                yQ_stand = yQ ; 

            elseif strcmp(flag_country, 'US')
                vintage_month = month(datasets(v).vintage) ; 

                y = datasets(v).data_gdp ;     
                meangdp = datasets(v).meangdp ; 
                stdgdp  = datasets(v).stdgdp ;       
                yQ_stand = y(~isnan(y)) ; 
            end

            % - check for rolling estimation window
            if flag_rolling == 1
                yQ_stand = yQ_stand(end-40+1:end) ; % 10-year rolling window
            end

            % determine forecast horizon
            % ----------------------------- 
            if ismember(vintage_month,[1 4 7 10])
                flag_backcast = 1 ; 
            else
                flag_backcast = 0 ; 
            end

            if ismember(vintage_month,[3 6 9 12])
                foreH = 2 + flag_backcast ;
            else
                foreH = 1 + flag_backcast ; 
            end   

            % - estimate B-AR(1)
            % ------------------------------------     

            dens = f_BayARpreddens(yQ_stand,foreH,Np,pi1,pi2,10,5,10);
            draws.nowcast = stdgdp * dens(1 + flag_backcast , :  ) + meangdp ; 
            if foreH == 2 + flag_backcast 
                draws.forecast = stdgdp * dens(2 + flag_backcast , : )  + meangdp ; 
            end

            save([ dir_local filename_save '\PH_' flag_country '_v' num2str(v) '_' name_flag_rolling '.mat'],'draws')
        end
    end
end

function store_dens = f_BayARpreddens(data,foreH,P,pi1,pi2,nreplic,nburnin,nthin)

	% calculate prior
	beta0 = zeros(P,1);
	pp = 1:P;
	Pbeta0 =  diag(1./(pi1./(pp.^pi2))); 
	upsilon_sig2 = 1;
	s_sig2 = 0.01;

	% empty matrices to store pred
	store_dens = NaN(foreH,nreplic/nthin);
	store_beta = NaN(P,nreplic/nthin);

	% starting values => OLS
	y = data(P+1:end,1);
	X = [];
	for p=1:P
		X = [X data(P+1-p:end-p,1)];
	end

	beta = (X'*X)\X'*y;
	sig2 = var(y - X*beta);
	tic
	for r = 1:nreplic+nburnin
		
		if mod(r,1000)==0
			disp('Number of iterations');
			disp(r);toc;
			disp('-------------------------------------------')
			disp('-------------------------------------------')
		end
		
		% conditional on sig2, sample beta
		% ========================================
		Vbeta = (Pbeta0 + X'*X/sig2)\speye(P);
		betahat = Vbeta*(Pbeta0*beta0 + X'*y/sig2);
		beta = betahat + chol(Vbeta,'lower')*randn(P,1);
		
		% conditional on beta, sample sig2
		% ========================================
		eps = y - X*beta;
		sig2 = 1/gamrnd(upsilon_sig2 + length(y)/2,1/(s_sig2 + eps'*eps/2));
		
		% conditional on params, sample forecasts
		% ========================================
		yfore = [y;NaN(foreH,1)];
		for h = 1:foreH
			Xfore = [];
			for p=1:P
				Xfore = [Xfore yfore(length(y)+h-p,1)];
			end
			yfore(length(y) + h,1) = Xfore * beta + sqrt(sig2) * randn;
		end

		% store draws
		% ========================================

		if r>nburnin && mod(r-nburnin,nthin) == 0
			saveindex = (r-nburnin)/nthin;
			store_dens(:,saveindex) = yfore(end-foreH + 1:end,1);
			store_beta(:,saveindex) = beta;
		end
	end
end










        

    

