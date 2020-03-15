function evaloptions = load_evaloptions(flag_country) 

evaloptions.Nrs = 1:10 ; % # of factors
evaloptions.Npriorspecs = 5 ; % 1: Normal, 2: MG, 3: PMNM, 4: HS++
evaloptions.Nthin = 1 ; % # use every Nthin-th draw
evaloptions.Ndraws = 1000 ; % # of draws
evaloptions.Nhs = 4 ; % # of horizons, i.e. nowcasts per quarter
evaloptions.Nmultpool = 1 ; % factor by which Ndraws is multiplied when pooling (to get "smoother" pools!)
if strcmp(flag_country, 'GER')
    evaloptions.Nquarters = 52; % # of quarters
    evaloptions.multfac = 1 ; 
    evaloptions.indexstarts = [1 1 17] ; % 2006Q1 2006Q1 2010Q1
    evaloptions.indexends =  [52 16 52] ; % 2018Q4 2009Q4 2018Q4
elseif strcmp(flag_country, 'US')
    evaloptions.Nquarters = 76 ; 
    evaloptions.multfac = 100 ;
    evaloptions.indexstarts = [1 1 41] ; % 2000Q1 2000Q1 2010Q1
    evaloptions.indexends =  [76 40 76] ; % 2018Q4 2009Q4 2018Q4
end

evaloptions.computelogscore = 'ksdensity' ; 

end
