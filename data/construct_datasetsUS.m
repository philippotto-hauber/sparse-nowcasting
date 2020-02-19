clear;close all;clc

addpath C:\Users\Philipp\Documents\IRIS_Tbx; irisstartup

% ------------------------------- %
% - input
% ------------------------------- %

filename_in = 'C:\Users\Philipp\Documents\Dissertation\sparse nowcasting\data\in' ; 
filename_out = 'C:\Users\Philipp\Documents\Dissertation\sparse nowcasting\data\out' ; 

first_vintage = '1999-12' ; 
start_year = 2000 ;
end_year = 2018 ;
n_year = end_year - start_year + 1 ; 

% ------------------------------- %
% - vintage names
% ------------------------------- %

% vintages
counter = 2 ; 
vintages = cell(1,n_year*12+1) ; 
vintages{1,1} = first_vintage ; 
for y = start_year:end_year
    for m = 1:12
        if m<10
            temp = [num2str(y) '-0' num2str(m)] ;
        else
            temp = [num2str(y) '-' num2str(m)] ; 
        end
        vintages{1,counter} = temp ; 
        counter = counter + 1 ; 
    end
end

% ------------------------------- %
% - loop to construct vintages
% ------------------------------- %

samplestart_year = 1985 ;
samplestart_month = 1 ; 
samplestart = samplestart_year + samplestart_month / 12 ; 

for v = 1 : length(vintages)
   
    vintageyear = year(vintages{v},'yyyy-mm') ; 
    vintagemonth = month(vintages{v},'yyyy-mm') ;     
    vintageday = eomday(vintageyear, vintagemonth) ; % only needed for gdp vintages!
    
    % ------------------------------------- % 
    % - FRED-MD data ---------------------- %
    % ------------------------------------- %
    csv_in = [filename_in '\FRED_MD\' vintages{v} '.csv'] ;   
    [data_fredmd, dates_fredmd] = f_load_FRED_MD_data(csv_in,samplestart,vintageyear,vintagemonth) ; 
    
    % ------------------------------------- % 
    % - PHILLY FED BOS -------------------- %
    % ------------------------------------- %
    
    % raw, unadjusted data (from Philly Fed webpage)
    [ data_bos_nsa_level , data_bos_nsa_diff ] = f_load_BOS_data(filename_in, vintageyear,vintagemonth,samplestart) ;
    
    % seasonally adjust data (=> X13-ARIMA using IRIS tool box!)
    ts_temp = tseries( datrange( mm(1985,1) , mm(vintageyear,vintagemonth) ), data_bos_nsa_level ) ;
    x12_temp =  x12( ts_temp ) ;

    data_bos_sa_level = x12_temp.Data ; 
    
    data_bos_sa_diff = data_bos_sa_level( 2 : end , : ) - data_bos_sa_level( 1 : end - 1 , : ) ; 
    
    % remove first (two) rows that were lost due to transformations of FRED-MD data
    data_bos_sa_level = data_bos_sa_level( 3 : end , : ) ; 
    data_bos_nsa_level = data_bos_nsa_level( 3 : end , : ) ; 
    data_bos_sa_diff  = data_bos_sa_diff( 2 : end , : ) ; 
    data_bos_nsa_diff = data_bos_nsa_diff( 2 : end , : ) ; 
    
    % ------------------------------------- % 
    % - gross domestic product ------------ %
    % ------------------------------------- %
   
    [data_gdp, dates_gdp, meangdp, stdgdp, truegdp] = f_load_gdp_data(filename_in, vintageyear,vintagemonth,vintageday,samplestart) ;
    
    % remove first (two) rows that were lost due to transformations of FRED-MD data
    data_gdp = data_gdp( 3 : end , 1 ) ; 
    
    
    % ------------------------------------- % 
    % - store data set in structure ------- %
    % ------------------------------------- %
    
    datasets(v).vintage =  [num2str(vintageyear) '-' num2str(vintagemonth) '-' num2str(vintageday)] ; 
    datasets(v).meangdp =  meangdp ;
    datasets(v).stdgdp =  stdgdp ;
    datasets(v).data_fredmd =  data_fredmd ;
    datasets(v).data_bos_nsa_level =  data_bos_nsa_level ;
    datasets(v).data_bos_nsa_diff =  data_bos_nsa_diff ;
    datasets(v).data_bos_sa_level =  data_bos_sa_level ;
    datasets(v).data_bos_sa_diff =  data_bos_sa_diff ;
    datasets(v).data_gdp =  data_gdp ;
    datasets(v).dates = dates_fredmd ;   
    datasets(v).truegdp = truegdp ; 
end

save([filename_out '\datasetsUS.mat'],'datasets')

%---------------------------------------------------------%
%- FUNCTIONS ---------------------------------------------%
%---------------------------------------------------------%

function [ yt_stand , dates_num ] = f_load_FRED_MD_data(csv_in,samplestart,vintageyear,vintagemonth)
% ------------------------------------- % 
% - Parts I and II from FRED codes 
% - https://research.stlouisfed.org/econ/mccracken/fred-databases/
% ------------------------------------- %

% Load data from CSV file
dum=importdata(csv_in,',');

% Variable names
series=dum.textdata(1,2:end);

% Transformation numbers
tcode=dum.data(1,:);

% Raw data
rawdata=dum.data(2:end,:);

% Month/year of final observation
final_datevec=datevec(dum.textdata(end,1));
final_month=final_datevec(2);
final_year=final_datevec(1);

% Dates (monthly) are of the form YEAR+MONTH/12
% e.g. March 1970 is represented as 1970+3/12
% Dates go from 1959:01 to final_year:final_month (see above)
dates_num = (1959+1/12:1/12:final_year+final_month/12)';
index_start = find( abs( dates_num - samplestart ) < 1e-10 ) ; 
rawdata = rawdata( index_start : end , : ) ;
dates_num = dates_num( index_start : end , : ) ;

% Transform raw data to be stationary using auxiliary function
% prepare_missing()
%ytold = prepare_missing( rawdata , tcode );
% change trafo code from 4 to 5 for housing variables 
housingvars = { 'HOUST' , 'HOUSTNE' , 'HOUSTMW' , 'HOUSTS' ,'HOUSTW' , 'PERMIT' , 'PERMITNE' , 'PERMITMW' , 'PERMITS' ,'PERMITW' } ;
for i = 1 : length( housingvars )  
    varindex = find( strcmp(series, housingvars{ i } ) == 1 ) ; 
    tcode( varindex ) = 5 ; 
end

yt = prepare_missing( rawdata , tcode );

% remove outliers => done before Gibbs Sampler is called!
%[yt , ~ ] = remove_outliers( yt ) ;

% Reduce sample to usable dates: remove first two months because some
% series have been first differenced
yt=yt( 3 : end , : ) ;
dates_num=dates_num( 3 : end , : ) ;

% add rows of NaN at the end of sample to account for difference between
% sample end and vintage month
diff_month = round(( ( vintageyear + vintagemonth/12 ) - dates_num(end) ) * 12) ;
yt = [yt; NaN( diff_month , size(yt,2) ) ] ; 
dates_num = [dates_num; dates_num(end) + (1:diff_month)'./12 ] ;

% standardize yt
yt_stand = NaN(size(yt)) ; 
for i = 1 : size(yt,2) 
    yt_stand(:,i) = ( yt(:,i) - nanmean(yt(:,i)) ) / nanstd(yt(:,i)) ; 
end
end

function yt   = prepare_missing(rawdata,tcode)
% =========================================================================
% DESCRIPTION: 
% This function transforms raw data based on each series' transformation
% code.
%
% -------------------------------------------------------------------------
% INPUT:
%           rawdata     = raw data 
%           tcode       = transformation codes for each series
%
% OUTPUT: 
%           yt          = transformed data
%
% -------------------------------------------------------------------------
% SUBFUNCTION:
%           transxf:    transforms a single series as specified by a 
%                       given transfromation code
%
% =========================================================================
% APPLY TRANSFORMATION:
% Initialize output variable
yt        = [];                                     

% Number of series kept
N = size(rawdata,2);  

% Perform transformation using subfunction transxf (see below for details)
for i = 1:N
    dum = transxf(rawdata(:,i),tcode(i));
    yt    = [yt, dum];
end

end

function y=transxf(x,tcode)
% =========================================================================
% DESCRIPTION:
% This function transforms a single series (in a column vector)as specified
% by a given transfromation code.
%
% -------------------------------------------------------------------------
% INPUT:
%           x       = series (in a column vector) to be transformed
%           tcode   = transformation code (1-7)
%
% OUTPUT:   
%           y       = transformed series (as a column vector)
%
% =========================================================================
% SETUP:
% Number of observations (including missing values)
n=size(x,1);
 
% Value close to zero 
small=1e-6;

% Allocate output variable
y=NaN*ones(n,1);

% =========================================================================
% TRANSFORMATION: 
% Determine case 1-7 by transformation code
switch(tcode);
    
  case 1, % Level (i.e. no transformation): x(t)
    y=x;

  case 2, % First difference: x(t)-x(t-1)
    y(2:n)=x(2:n,1)-x(1:n-1,1);
  
  case 3, % Second difference: (x(t)-x(t-1))-(x(t-1)-x(t-2))
    y(3:n)=x(3:n)-2*x(2:n-1)+x(1:n-2);

  case 4, % Natural log: ln(x)
    if min(x) < small; 
        y=NaN; 
    else
        y=log(x);
    end;
  
  case 5, % First difference of natural log: ln(x)-ln(x-1)
    if min(x) > small;
        x=log(x);
        y(2:n)=x(2:n)-x(1:n-1);
    end;
  
  case 6, % Second difference of natural log: (ln(x)-ln(x-1))-(ln(x-1)-ln(x-2))
    if min(x) > small;
        x=log(x);
        y(3:n)=x(3:n)-2*x(2:n-1)+x(1:n-2);
    end;
  
  case 7, % First difference of percent change: (x(t)/x(t-1)-1)-(x(t-1)/x(t-2)-1)
    y1(2:n)=(x(2:n)-x(1:n-1))./x(1:n-1);
    y(3:n)=y1(3:n)-y1(2:n-1);
end
end

function [Y,n]=remove_outliers(X)
% =========================================================================
% DESCRIPTION:
% This function takes a set of series aligned in the columns of a matrix
% and replaces outliers with the value NaN.
%
% -------------------------------------------------------------------------
% INPUT:
%           X   = dataset (one series per column)
% 
% OUTPUT:
%           Y   = dataset with outliers replaced with NaN 
%           n   = number of outliers found in each series
%
% -------------------------------------------------------------------------
% NOTES:
%           1) Outlier definition: a data point x of a series X(:,i) is
%           considered an outlier if abs(x-median)>10*interquartile_range.
%
%           2) This function ignores values of NaN and thus is capable of
%           replacing outliers for series that have missing values.
%
% =========================================================================
% FUNCTION:

% Calcualte median of each series
median_X=nanmedian(X,1);

% Repeat median of each series over all data points in the series
median_X_mat=repmat(median_X,size(X,1),1);

% Calculate quartiles 
Q=prctile(X,[25, 50, 75],1);

% Calculate interquartile range (IQR) of each series
IQR=Q(3,:)-Q(1,:);

% Repeat IQR of each series over all data points in the series
IQR_mat=repmat(IQR,size(X,1),1);

% Determine outliers 
Z=abs(X-median_X_mat);
outlier=Z>(10*IQR_mat);

% Replace outliers with NaN
Y=X;
Y(outlier)=NaN;

% Count number of outliers
n=sum(outlier,1);
end

function [ data_bos_level , data_bos_diff ] = f_load_BOS_data(filename_in, vintageyear,vintagemonth,samplestart)

% load data
[num,txt,~] = xlsread([filename_in '\bos_historyx.xls']) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------------------------------------------------------------------- %
% -- MANUAL INPUT - NEEDS TO BE ADJUSTED IF XLS FILE CHANGES -------- %
firstyear = 1968 ; firstmonth = 5 ; 
finalyear = 2019 ; finalmonth = 1 ;
% ------------------------------------------------------------------- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dates_bos = (firstyear + firstmonth/12 : 1/12 : finalyear + finalmonth/12)' ; 
indexstart = find( abs((dates_bos - samplestart)) < 1e-10 ) ; % why does simple find not work ????
indexend = find( abs(dates_bos - (vintageyear + vintagemonth/12)) < 1e-10 ) ;
namevars = {'gacdfna','gafdfna','nocdfna','nofdfna','shcdfna','shfdfna',...
            'uocdfna','uofdfna','dtcdfna','dtfdfna','ivcdfna','ivfdfna',...
            'ppcdfna','ppfdfna','prcdfna','prfdfna','necdfna','nefdfna',...
            'awcdfna','awfdfna'} ;  
[~, ~, indexvars] = intersect( namevars , txt(1,2:end) , 'stable' ) ; 
data_bos_level = num( indexstart : indexend , indexvars ) ; 
data_bos_diff = diff( num( indexstart : indexend , indexvars ) ) ; 
end

function [data_gdp_stand, dates_gdp, meangdp, stdgdp, truegdp] = f_load_gdp_data(filename_in, vintageyear,vintagemonth,vintageday,samplestart)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ----------------------------------------------------------------------- %
% -- MANUAL INPUT - NEEDS TO BE ADJUSTED IF XLS FILE CHANGES ------------ %
firstyear = 1980 ; firstmonth = 1 ; % firstmonth = ceil(firstquarter/3)
finalyear = 2018 ; finalmonth = 9 ; % finalmonth = finalquarter * 3
% ----------------------------------------------------------------------- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% load data
[num,txt,~] = xlsread([filename_in, '\vintagesGDP.xls'],'Vintages Starting 1999-01-29') ;

% get vintage dates
for n = 1 : size(txt,2) - 1 
    temp = txt{1,n+1} ; 
    numvintagedates(n,1) = datenum(temp(end-7:end),'yyyymmdd') ;
end
numcurrentvintage = datenum([num2str(vintageyear) '-' num2str(vintagemonth) '-' num2str(vintageday)],'yyyy-mm-dd') ; 

% extract gdp data   
indexvintage = sum(numvintagedates <= numcurrentvintage) ;
gdptemp = num(:,indexvintage) ;

%  convert to monthly frequency
gdptemp_m = kron( gdptemp , [NaN;NaN;1] ) ; 
datesgdp_m = (firstyear + firstmonth/12 : 1/12 : finalyear + finalmonth/12)' ;

% log diff
gdptemp_m = log(gdptemp_m(4:end,1)) - log(gdptemp_m(1:end-3,1)) ; 
datesgdp_m = datesgdp_m(4:end,:) ; 

% adjust sample
indexstart = find(datesgdp_m == samplestart) ; 
indexend = find(datesgdp_m == vintageyear + vintagemonth/12) ; 
data_gdp = gdptemp_m(indexstart:indexend,1) ; 
dates_gdp = datesgdp_m(indexstart:indexend,:) ; 

% standardize data_gdp
meangdp = nanmean(data_gdp) ; 
stdgdp = nanstd(data_gdp) ; 
data_gdp_stand = ( data_gdp - nanmean(data_gdp) ) / nanstd(data_gdp) ; 

% get "true" gdp according to first estimate and last available vintage
temp = num(:,indexvintage) ;
indexrow_nowcast = length(temp) - sum(isnan(temp)) + 1 ; 
indexcol_nowcast = sum(isnan(num(indexrow_nowcast,:))) + 1 ; 
truegdp.first.nowcast = log(num(indexrow_nowcast,indexcol_nowcast)) - log(num(indexrow_nowcast-1,indexcol_nowcast)) ; 
truegdp.second.nowcast = log(num(indexrow_nowcast,indexcol_nowcast+1)) - log(num(indexrow_nowcast-1,indexcol_nowcast+1)) ; 
truegdp.final.nowcast = log(num(indexrow_nowcast,end)) - log(num(indexrow_nowcast - 1,end)) ; 

if ismember(vintagemonth,[12 9 6 3])
    indexrow_forecast = indexrow_nowcast + 1 ; 
    indexcol_forecast = sum(isnan(num(indexrow_forecast,:))) + 1 ; 
    truegdp.first.forecast = log(num(indexrow_forecast,indexcol_forecast)) - log(num(indexrow_forecast-1,indexcol_forecast)) ; 
    truegdp.second.forecast = log(num(indexrow_forecast,indexcol_forecast+1)) - log(num(indexrow_forecast-1,indexcol_forecast+1)) ; 
    truegdp.final.forecast = log(num(indexrow_forecast,end)) - log(num(indexrow_forecast - 1,end)) ; 
end
end






