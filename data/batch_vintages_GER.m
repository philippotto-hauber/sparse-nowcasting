clear; close all; clc
% needs to be run from C:\Users\Philipp\Google Drive\Prognose\Echtzeitdatensatz

cd 'C:\Users\Philipp\Google Drive\Prognose\Echtzeitdatensatz'
filename_out = 'C:\Users\Philipp\Documents\Dissertation\sparse nowcasting\data\in'

[ vintages , Nvintages ] = f_getvintages( 2005 , 12 , 2018 , 12 ) ;
        
for v = 1 : length(vintages)
    v
    datasets.vintage(v).data_ifo = f_load_ifo(vintages{v}) ;
    datasets.vintage(v).data_ESIBCI = f_load_ESIBCI(vintages{v}) ;
    datasets.vintage(v).data_BuBaRTD = f_load_BuBaRTD(vintages{v}) ;
    datasets.vintage(v).data_financial = f_load_financial(vintages{v}) ;
    datasets.vintage(v).vintagedate = vintages{v} ; 
end

% save datasets to mat
save([filename_out '\datasets_orig'],'datasets');

    
