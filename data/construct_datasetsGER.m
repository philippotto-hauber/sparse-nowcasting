clear ; close all ; clc ; 
% ----------------------------------------------------------------------- %

% ------------------------------- %
% - input
% ------------------------------- %

filename_in = 'C:\Users\Philipp\Documents\Dissertation\sparse nowcasting\data\in' ; 
filename_out = 'C:\Users\Philipp\Documents\Dissertation\sparse nowcasting\data\out' ; 

first_vintage = '2005-12' ; 
start_year = 2006 ;
end_year = 2018 ;
n_years = 2018 - 2006 + 1 ; 

% ------------------------------- %
% - get vintages
% ------------------------------- %
vintages = cell(1,n_years*12 + 1) ; 
vintages{1,1} = first_vintage ; 
counter = 2 ; 
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
Nvintages = length( vintages ) ; 

sample_start = 1992 + 1/12 ; 

% load raw vintages
load([filename_in '\datasets_orig.mat'])
datasets_temp = datasets ; clear datasets ; 

% load IRIS
addpath C:\Users\Philipp\Documents\IRIS_Tbx; irisstartup
%addpath C:\Users\Hauber\Documents\IRIS_Tbx; irisstartup

for v = 1 : Nvintages
    
    % ----------------------------------------------------------------- % 
    % - handle dates
    % ----------------------------------------------------------------- %    
    dates_num = sample_start : 1/ 12 : ( ( year( vintages{ v } ) + month( vintages{ v } ) / 12 ) ) ;
    dates_master = f_convertdates( dates_num' ) ; 
   
    % ----------------------------------------------------------------- % 
    % - Bundesbank real-time data (production, orders, turnover etc.)
    % ----------------------------------------------------------------- %    
    
    % monthly vars
    % ------------------
    index_monthly = strcmp( datasets_temp.vintage(v).data_BuBaRTD.type , 'm' ) ; 
    index_samplestart = find( abs( datasets_temp.vintage(v).data_BuBaRTD.dates - sample_start ) < 1e-5 ) ;
    data_bubartd = datasets_temp.vintage(v).data_BuBaRTD.data( index_samplestart : end , index_monthly ) ; 
    
    % gross domestic product
    % ------------------
    index_gdp = strcmp( datasets_temp.vintage(v).data_BuBaRTD.names , 'gross domestic product' ) ;
    data_gdp = datasets_temp.vintage(v).data_BuBaRTD.data( index_samplestart : end , index_gdp ) ; 
    
    % ----------------------------------------------------------------- % 
    % - financial market data (exchange rates, interest rates etc.)
    % ----------------------------------------------------------------- % 
    
    index_samplestart = find( abs( datasets_temp.vintage(v).data_financial.dates - sample_start ) < 1e-5 ) ;
    data_financial = datasets_temp.vintage(v).data_financial.data( index_samplestart : end , : ) ; 
    
    % ----------------------------------------------------------------- % 
    % - ifo (headline, manufacturing, retail and wholesale)
    % ----------------------------------------------------------------- % 
    index_samplestart = find( abs( datasets_temp.vintage(v).data_ifo.dates - sample_start ) < 1e-5 ) ;
    
    ts_temp = tseries( datrange( mm( floor( sample_start ) , round( ( sample_start - floor( sample_start ) ) * 12 ) ) , mm(year(vintages{v}),month(vintages{v})) ), datasets_temp.vintage(v).data_ifo.rawdata( index_samplestart : end , 1 : 37 ) ) ; % exclude services and quarterly series on capacity utilisation
    x12_temp =  x12( ts_temp ) ;

    data_ifo_level = x12_temp.Data ; 
    
    data_ifo_diff = data_ifo_level( 2 : end , : ) - data_ifo_level( 1 : end - 1 , : ) ; 
    
    % ------------------------------------------------------ % 
    % - adjust number of observations
    % ------------------------------------------------------ % 
    
     maxobs = max( [ size(data_bubartd , 1 ) , ...
                size(data_financial , 1 ) , ...
                size(data_ifo_level , 1 )  , ...
                size(data_gdp , 1 )  ] ) ; 
    
    datasets.vintage(v).data_ifo_level = [ data_ifo_level ; NaN( maxobs - size( data_ifo_level , 1 ) , size( data_ifo_level , 2 ) ) ] ;          
    datasets.vintage(v).data_ifo_diff = diff( data_ifo_level ) ;  
    datasets.vintage(v).data_financial = [ data_financial ; NaN( maxobs - size( data_financial , 1 ) , size( data_financial , 2 ) ) ] ;          
    datasets.vintage(v).data_bubartd = [ data_bubartd ; NaN( maxobs - size( data_bubartd , 1 ) , size( data_bubartd , 2 ) ) ] ;        
    datasets.vintage(v).data_gdp = [ data_gdp ; NaN( maxobs - size( data_gdp , 1 ) , size( data_gdp , 2 ) ) ] ; 
    
    % names, dates, ....
    datasets.vintage( v ).vintage = vintages{ v } ; 
    datasets.vintage( v ).dates = dates_master( 1 : maxobs ) ;
end
save([filename_out '\datasetsGER.mat'],'datasets')

% ----------------------------------------------------------------------- %
% -- FUNCTIONS ---------------------------------------------------------- %
% ----------------------------------------------------------------------- %

function dates_str = f_convertdates( dates_num ) 

dates_str = cell( length( dates_num ) , 1 ) ; 
for t = 1 : length( dates_num ) 
    date_year = floor( dates_num( t ) ) ;
    date_month = ( dates_num( t ) - floor( dates_num( t ) ) ) * 12 ; 
    if date_month == 0
        temp = [num2str(date_year - 1) '-' num2str(12)]  ;  
    else
        temp = [num2str(date_year) '-' num2str(date_month)] ; 
    end
    dates_str{ t , 1 } = temp ; 
end

end

