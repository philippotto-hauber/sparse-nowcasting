clear;close all;clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This code extracts the first, second and final estimates of quarterly
%%% GDP growth, storing the results in a mat-file called truegdp.mat 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ------------------------------------
% ----- set years and quarters -------
% ------------------------------------

years = 2006:2018 ; 
quarters = 1:4 ; 

% --------------------------------------------
% ----- dir to store structure in  -----------
% --------------------------------------------

dir_in = 'C:\Users\Philipp\Documents\Dissertation\sparse nowcasting\data\in' ;
dir_out = 'C:\Users\Philipp\Documents\Dissertation\sparse nowcasting\data\out' ;


% ---------------------------
% ----- load GDP vintages ---
% ---------------------------

[num,txt,~] = xlsread([dir_in '\Vintages_GDP_GER.xls']) ;
dates = txt((size(txt,1)-size(num,1)-1):end,1) ;
counter_q = 1 ; 
for y = years 
    for q = quarters
       truegdp_strct.quarters{counter_q} = [num2str(y) 'Q' num2str(q)] ; 
       index_row = find(strcmp(dates,datestr(datenum([y,q*3-2,1]),'yyyy-mm'))==1) ; 
       [index_col_first, index_col_second, index_col_final] = f_find_indexcols(num,index_row) ; 

       truegdp_strct.first(counter_q) = 100 * ( log(num(index_row,index_col_first)) - log(num(index_row-1,index_col_first ) ) ) ; 
       truegdp_strct.second(counter_q) = 100 * ( log(num(index_row,index_col_second)) - log(num(index_row-1,index_col_second ) ) ) ;
       truegdp_strct.final(counter_q) = 100 * ( log(num(index_row,index_col_final)) - log(num(index_row-1,index_col_final ) ) ) ;
       
       counter_q = counter_q + 1 ; 
    end
end

save([dir_out '\truegdpGER.mat'],'truegdp_strct')

function [index_col_first, index_col_second, index_col_final] = f_find_indexcols(num,index_row)

	% index_col_first
	index_col_first = 1 ;
	while isnan(num(index_row,index_col_first))
		index_col_first = index_col_first + 1 ;
	end

	% index_col_second and index_col_final
	index_col_second = index_col_first + 1 ; 
	index_col_final = size(num,2) ; 
end

        