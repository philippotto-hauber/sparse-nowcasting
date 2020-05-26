function specnum = f_get_specnum(v, Nprior, Nr, Np, flag_survey, flag_recroll, flag_country) 
% this function returns the unique spec number relating to the given inputs
% of spec variablesthat is used to call the compiled functions on the 
% server. Useful for finding the spec number of a given spec that needs to
% be executed locally. 

if strcmp(flag_country, 'US')
    Nvintages = 1:219;
elseif strcmp(flag_country, 'GER')
    Nvintages = 1:157; 
else
    error('Incorrect country flag (options:GER/US)')
end

Nrs = 1:10 ; 
Npriors = 1:5 ; 
Nps = [1 3] ;
Nmod = 1:4 ;

% calculate combinations of specifications
Nrs_temp = repmat(Nrs, 1, length(Npriors) * length(Nvintages) * length(Nps) * length(Nmod)) ; 
Npriors_temp = repmat(kron(Npriors, ones(1, length(Nrs))), 1 , length(Nvintages) * length(Nps) * length(Nmod)) ;
Nvintages_temp = repmat(kron(Nvintages, ones(1, length(Nrs) * length(Npriors))), 1, length(Nps) * length(Nmod)) ; 
Nps_temp = repmat(kron(Nps, ones(1, length(Nrs) * length(Npriors) * length(Nvintages))), 1, length(Nmod)) ;
Nmod_temp = kron(Nmod, ones(1, length(Nrs) * length(Npriors) * length(Nvintages) * length(Nps))) ;

% map flags to Nmod
if strcmp(flag_survey, 'level') && strcmp(flag_recroll, 'rec')
    Nmod = 1;
elseif strcmp(flag_survey, 'diff') && strcmp(flag_recroll, 'rec')
    Nmod = 2 ; 
elseif strcmp(flag_survey, 'level') && strcmp(flag_recroll, 'rolling')
    Nmod = 3;
elseif strcmp(flag_survey, 'diff') && strcmp(flag_recroll, 'rolling')
    Nmod = 4 ; 
else
    error('incorrect flag_survey or flag_recroll')
end
    
% find spec number
specnum = find( Nrs_temp == Nr & Nvintages_temp == v & Npriors_temp == Nprior & Nps_temp == Np & Nmod_temp == Nmod) ;



