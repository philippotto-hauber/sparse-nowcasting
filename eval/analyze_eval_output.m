clear; close all; 
h = 4;

R = 5;
dir_in = '../../../../Desktop/results_eval mat files/';
flag_country = 'GER';

load([dir_in 'results_eval_' flag_country '_rec_level_Np1_first.mat'])

bar = results_eval.benchmark_BAR.horizon(h).sfe;
nig = results_eval.priors(1).R(R).horizon(h).sfe;
nd = results_eval.priors(5).R(R).horizon(h).sfe;
pmnm = results_eval.priors(3).R(R).horizon(h).sfe;
hs = results_eval.priors(4).R(R).horizon(h).sfe;



disp([sqrt(mean(nd)) / sqrt(mean(bar));...
sqrt(mean(nig)) / sqrt(mean(bar));...
sqrt(mean(pmnm)) / sqrt(mean(bar));...
sqrt(mean(hs)) / sqrt(mean(bar))])

if strcmp(flag_country, 'GER')
ind_postcrisis = 17:52;
elseif strcmp(flag_country, 'US')
    ind_postcrisis = 41:76;
end

disp([sqrt(mean(nd(ind_postcrisis))) / sqrt(mean(bar(ind_postcrisis)));...
sqrt(mean(nig(ind_postcrisis))) / sqrt(mean(bar(ind_postcrisis)));...
sqrt(mean(pmnm(ind_postcrisis))) / sqrt(mean(bar(ind_postcrisis)));...
sqrt(mean(hs(ind_postcrisis))) / sqrt(mean(bar(ind_postcrisis)))])