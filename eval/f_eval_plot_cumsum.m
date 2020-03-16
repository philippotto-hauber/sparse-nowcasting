function f_eval_plot_cumsum(flag_survey, flag_sample, flag_truegdp, Np, flag_country)

    % - options ---
    % ----------------------------
    evaloptions = load_evaloptions(flag_country) ; 

    % - directories ---
    % ----------------------------
    dir_load = ['C:\Users\Philipp\Documents\Dissertation\sparse nowcasting\eval\' flag_country '\' flag_survey ' ' flag_sample '\Np = ' num2str(Np) '\' flag_truegdp '\'] ;
    dir_save = [dir_load 'graphs\'] ; 
    if exist(dir_save, 'dir') ~= 7;mkdir(dir_save); end  

    % - load forecast structure
    % -----------------
    load([dir_load 'results_eval.mat'])
    
    % - filenames and 
    % --------------------------
    filename_ps = [dir_save 'graphs_cumsum.ps'] ; 
    filename_pdf = [dir_save 'graphs_cumsum.pdf'] ; 
    
    if isfile(filename_pdf); delete(filename_pdf); end
    
    % - graph labels
    % --------------------------
     xtickslabelname = evaloptions.quarters_graph' ; 
     
    % - metrics
    % --------------------------
    metrics = { 'sfe' , 'logscore' , 'crps' }  ;
    Nmetrics = length( metrics ) ;  
    
    % - start looping
    % --------------------------
    
    for index_r = 1 : length( evaloptions.Nrs ) + 1
        if index_r <= length( evaloptions.Nrs )
            r = Nrs( index_r ) ; 
        end

        for m = 1 : Nmetrics

            metric = metrics{ m } ;
            if strcmp( metric , 'sfe' )
                ylabelname = 'cumul. squared forecast error' ; 
            elseif strcmp( metric , 'logscore' )
                ylabelname = 'cumul. log score' ; 
            elseif strcmp( metric , 'crps' )
                ylabelname = 'cumul. CRPS' ; 
            end

            for i = 1 : length(evaloptions.indexstarts)

                figure
                fig = gcf;
                fig.PaperOrientation = 'portrait';

                indexstart = evaloptions.indexstarts(i) ; 
                indexend = evaloptions.indexends(i) ; 

                xtickslabelname = evaloptions.quarters_graph(indexstart:indexend)' ; 

                for h = 1:evaloptions.Nhs
                    if index_r == ( length( evaloptions.Nrs ) + 1 )
                        titlename = [name_country_long ', pool, ' results_eval.priors(1).pool.horizon(h).name ' (in months)'] ;
                    else
                        titlename = [name_country_long ', R=' num2str(r) ', ' results_eval.priors(1).R(1).horizon(h).name ' (in months)'] ;
                    end
                    subplot(Nhs,1,h)

                    temp_metric_benchmark = results_eval.benchmark_BAR.horizon(h).(metric)(indexstart:indexend) ; 

                    for p = evaloptions.Npriorspecs
                        if index_r == ( length( evaloptions.Nrs ) + 1 )
                            temp_metric = [temp_metric results_eval.priors(p).pool.horizon(h).(metric)(indexstart:indexend)] ; 
                        else
                            temp_metric = [temp_metric results_eval.priors(p).R(r).horizon(h).(metric)(indexstart:indexend)] ; 
                        end
                    end

                    % - plot data
                    % -----------------
                    f_plot_cumul( temp_metric - temp_metric_benchmark , titlename , ylabelname , xtickslabelname , metric )
                end
                
                % - save plot
                % -----------------        
                print(filename_ps,'-dpsc','-fillpage','-append') ;   
                close
            end
        end
    end

    % -convert ps figure to pdf
    % ------------------------------ 
    %
    ps2pdf('psfile', filename_ps , 'pdffile', filename_pdf , ...
       'gspapersize', 'a4', ...
       'gscommand', 'C:\Program Files\gs9.26\bin\gswin64c.exe', ...
       'gsfontpath', 'C:\Program Files\gs9.26\Resource\Font', ...
       'gslibpath', 'C:\Program Files\gs9.26\lib')
   delete(filename_ps) 
end

function f_plot_cumul(data,titlename,ylabelname,xtickslabelname,metric)

    % cumulate data
    data_cumsum = cumsum(data) ; 

    p1 = plot(data_cumsum(:,1),'Color',[0 0 0],'LineStyle',':','LineWidth',1.5);
    hold on
    p2 = plot(data_cumsum(:,2),'Color',[1 .5 0],'LineStyle','-','LineWidth',1.5);
    p3 = plot(data_cumsum(:,3),'Color',[0.6 0 0.6],'LineStyle','-','LineWidth',1.5);
    p4 = plot(data_cumsum(:,4),'Color',[0 0 1],'LineStyle','-','LineWidth',1.5);
    set(gca, 'FontSize', 7)

    % flip y-axis for squared forecast errors and CRPS
    if strcmp( metric , 'sfe' ) || strcmp( metric , 'crps' )
        set(gca, 'YDir','reverse')
    end

    box off
    title(titlename,'Fontsize',10,'FontWeight','normal')
    lgd = legend([p1 p2 p3 p4 p5 p6], ...
           'NIG', ...
           'MG',...
           'PMNM',...
           'HS+',...
           'Nd',...
           'Location','best') ;
    lgd.FontSize = 7;  
    ylabel(ylabelname,'Fontsize',7)
    xticks(1:size(data,1)) 
    xticklabels(xtickslabelname)
end

