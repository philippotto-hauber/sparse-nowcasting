function f_eval_plot_pointnowcasts(flag_survey, flag_sample, flag_truegdp, Np, flag_country) 

    % - options ---
    % ----------------------------
    evaloptions = load_evaloptions(flag_country) ; 

    % - directories ---
    % ----------------------------
    dir_save = [pwd '/' flag_country '/' flag_survey ' ' flag_sample '/Np = ' num2str(Np) '/' flag_truegdp '/graphs/'] ; 
    if exist(dir_save, 'dir') ~= 7;mkdir(dir_save); end  
    dir_truegdp = '' ;

    % - load forecast structure
    % -----------------
    load(['results_eval mat files/results_eval_' flag_country '_' flag_sample '_' flag_survey '_Np' num2str(Np) '_' flag_truegdp '.mat'])

    % - load true gdp mat-file ---
    % ----------------------------
    load([dir_truegdp 'truegdp' flag_country '.mat'])
    
    % - filenames and 
    % --------------------------
    filename_ps = [dir_save 'graphs_pointnowcasts.ps'] ; 
    filename_pdf = [dir_save 'graphs_pointnowcasts.pdf'] ; 
    
    if isfile(filename_pdf); delete(filename_pdf); end
    
    % - graph labels
    % --------------------------
    ylabelname = 'percent' ; 
    xtickslabelname = evaloptions.quarters_graph' ; 
    
    % - start looping
    % --------------------------

    for index_r = 1 : length( evaloptions.Nrs ) + 1
        if index_r <= length(evaloptions.Nrs)
            r = index_r ; 
        end
        
        for i = 1 : length(evaloptions.indexstarts)

            indexstart = evaloptions.indexstarts(i) ; 
            indexend = evaloptions.indexends(i) ; 

            if strcmp( flag_truegdp , 'first')
                truegdp = evaloptions.multfac*truegdp_strct.first( indexstart : indexend ) ; 
            elseif strcmp( flag_truegdp , 'second')
                truegdp = evaloptions.multfac*truegdp_strct.second( indexstart : indexend ) ; 
            elseif strcmp( flag_truegdp , 'final')
                truegdp = evaloptions.multfac*truegdp_strct.final( indexstart : indexend ) ;
            end

            dens_medians = NaN( 6 , length( truegdp ) ) ; 
            dens_upper = NaN( 6 , length( truegdp ) ) ; 
            dens_lower = NaN( 6 , length( truegdp ) ) ; 

            for h = 1:evaloptions.Nhs

                figure
                fig = gcf;
                fig.PaperOrientation = 'landscape';
                if index_r == length( evaloptions.Nrs ) + 1
                    titlename = [flag_country ', pool, ' results_eval.priors(1).R(1).horizon(h).name ' (in months)'] ;
                    for p = evaloptions.Npriorspecs
                        counter_q = 1 ; 
                        for q = indexstart : indexend
                            dens_medians(p,counter_q) = median( results_eval.priors(p).pool.horizon(h).dens{q} ) ;
                            dens_upper(p,counter_q) = prctile( results_eval.priors(p).pool.horizon(h).dens{q} , 75 ) ;
                            dens_lower(p,counter_q) = prctile( results_eval.priors(p).pool.horizon(h).dens{q} , 25 ) ;
                            counter_q = counter_q + 1 ; 
                        end
                    end
                else
                    titlename = [flag_country ', R=' num2str(r) ', ' results_eval.priors(1).R(1).horizon(h).name ' (in months)'] ;
                    counter_p = 1;
                    for p = evaloptions.Npriorspecs
                        counter_q = 1 ; 
                        for q = indexstart : indexend
                            dens_medians(counter_p,counter_q) = median( results_eval.priors(p).R(r).horizon(h).dens{q} ) ;
                            dens_upper(counter_p,counter_q) = prctile( results_eval.priors(p).R(r).horizon(h).dens{q} , 75 ) ;
                            dens_lower(counter_p,counter_q) = prctile( results_eval.priors(p).R(r).horizon(h).dens{q} , 25 ) ;
                            counter_q = counter_q + 1 ; 
                        end
                        counter_p = counter_p + 1;
                    end
                end

                % append B-AR benchmark
                dens_medians(counter_p + 1,:) = median( results_eval.benchmark_BAR.horizon(h).dens( indexstart : indexend , : ) , 2 ) ;
                %dens_upper(p+1,:) = prctile( results_eval.benchmark_BAR.horizon(h).dens( indexstart : indexend , : ) , 75 , 2 ) ;
                %dens_lower(p+1,:) = prctile( results_eval.benchmark_BAR.horizon(h).dens( indexstart : indexend , : ) , 25 , 2 ) ;

                % - create plot
                % ----------------- 
                f_plot_nowcast_realization( dens_medians , [] , [] , truegdp , flag_truegdp , titlename,ylabelname,xtickslabelname( indexstart : indexend ) )

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
%     ps2pdf('psfile', filename_ps , 'pdffile', filename_pdf , ...
%        'gspapersize', 'a4', ...
%        'gscommand', '/usr/bin/gs', ...
%        'gsfontpath', 'usr/share/ghostscript/9.07/lib', ...
%        'gslibpath', 'usr/share/ghostscript/9.07/lib')
    %delete(filename_ps) 
end

function f_plot_nowcast_realization( means , lower , upper , truegdp , flag_truegdp ,  titlename,ylabelname,xtickslabelname )

% benchmark comes last!
max_global = max(max(max(means)),max(truegdp)) ; 
min_global = min(min(min(means)),min(truegdp)) ; 
range_global = max_global - min_global ; 

%ylim([-2 2])
subplot(2,1,1)
% ----------------------------------------------------------------------- %
% - 2000Q1-2008Q4
% ----------------------------------------------------------------------- %

% p1 = plot(means(1,:),'Color',[0.9 0 0],'LineStyle','-','LineWidth',1.5);
p1 = plot(truegdp(1:length(truegdp)/2),'d','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0],'MarkerSize',6);
hold on
%plot(truegdp,'Color',[0 0 0],'LineStyle','-','LineWidth',0.1);

p2 = plot(means(1,1:length(truegdp)/2),'o','MarkerEdgeColor',[0.9 0 0],'MarkerSize',8);
plot(means(1,1:length(truegdp)/2),'Color',[0.9 0 0],'LineStyle',':','LineWidth',0.01);
%plot(upper(1,:),'Color',[0.9 0 0],'LineStyle',':','LineWidth',0.8);
%plot(lower(1,:),'Color',[0.9 0 0],'LineStyle',':','LineWidth',0.8);

p3 = plot(means(2,1:length(truegdp)/2),'o','MarkerEdgeColor',[1 .5 0],'MarkerSize',8);
plot(means(2,1:length(truegdp)/2),'Color',[1 .5 0],'LineStyle',':','LineWidth',0.01);
%plot(upper(2,:),'Color',[1 .5 0],'LineStyle',':','LineWidth',0.8);
%plot(lower(2,:),'Color',[1 .5 0],'LineStyle',':','LineWidth',0.8);

p4 = plot(means(3,1:length(truegdp)/2),'o','MarkerEdgeColor',[0.6 0 0.6],'MarkerSize',8);
plot(means(3,1:length(truegdp)/2),'Color',[0.6 0 0.6],'LineStyle',':','LineWidth',0.01);
%plot(upper(3,:),'Color',[0.6 0 0.6],'LineStyle',':','LineWidth',0.8);
%plot(lower(3,:),'Color',[0.6 0 0.6],'LineStyle',':','LineWidth',0.8);

p5 = plot(means(4,1:length(truegdp)/2),'o','MarkerEdgeColor',[0 0 1],'MarkerSize',8);
plot(means(4,1:length(truegdp)/2),'Color',[0 0 1],'LineStyle',':','LineWidth',0.01);
%plot(upper(4,:),'Color',[0 0 1],'LineStyle',':','LineWidth',0.8);
%plot(lower(4,:),'Color',[0 0 1],'LineStyle',':','LineWidth',0.8);

p6 = plot(means(5,1:length(truegdp)/2),'o','MarkerEdgeColor',[0 0.5 0.5],'MarkerSize',8);
plot(means(5,1:length(truegdp)/2),'Color',[0 0 0],'LineStyle',':','LineWidth',0.01);

p7 = plot(means(6,1:length(truegdp)/2),'o','MarkerEdgeColor',[0 0 0],'MarkerSize',8);
plot(means(6,1:length(truegdp)/2),'Color',[0 0 0],'LineStyle',':','LineWidth',0.01);


xt = get(gca, 'XTick');
set(gca, 'FontSize', 7)

box off
grid on
title(titlename,'Fontsize',10,'FontWeight','normal')
% lgd = legend([p1 p2 p3 p4 p5 p6], ...
%        [flag_truegdp ' release'],...
%        'NIG', ...
%        'MG',...
%        'PMNM',...
%        'HS+',...
%        'B-AR(1)',...
%        'Location','best') ;
% lgd.FontSize = 7;
ylim([min_global - range_global/10 max_global + range_global/10])
ylabel(ylabelname,'Fontsize',7)
xticks(1:1:length(truegdp)/2) 
xticklabels(xtickslabelname)

subplot(2,1,2)
% ----------------------------------------------------------------------- %
% - 2009Q1-2016Q4
% ----------------------------------------------------------------------- %

% p1 = plot(means(1,:),'Color',[0.9 0 0],'LineStyle','-','LineWidth',1.5);
p1 = plot(truegdp(length(truegdp)/2 + 1 : end),'d','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0],'MarkerSize',6);
hold on
%plot(truegdp,'Color',[0 0 0],'LineStyle','-','LineWidth',0.1);

p2 = plot(means(1,length(truegdp)/2 + 1 : end),'o','MarkerEdgeColor',[0.9 0 0],'MarkerSize',8);
plot(means(1,length(truegdp)/2 + 1 : end),'Color',[0.9 0 0],'LineStyle',':','LineWidth',0.01);
%plot(upper(1,:),'Color',[0.9 0 0],'LineStyle',':','LineWidth',0.8);
%plot(lower(1,:),'Color',[0.9 0 0],'LineStyle',':','LineWidth',0.8);

p3 = plot(means(2,length(truegdp)/2 + 1 : end),'o','MarkerEdgeColor',[1 .5 0],'MarkerSize',8);
plot(means(2,length(truegdp)/2 + 1 : end),'Color',[1 .5 0],'LineStyle',':','LineWidth',0.01);
%plot(upper(2,:),'Color',[1 .5 0],'LineStyle',':','LineWidth',0.8);
%plot(lower(2,:),'Color',[1 .5 0],'LineStyle',':','LineWidth',0.8);

p4 = plot(means( 3 , length(truegdp)/2 + 1 : end ),'o','MarkerEdgeColor',[0.6 0 0.6],'MarkerSize',8);
plot(means(3, length(truegdp)/2 + 1 : end ),'Color',[0.6 0 0.6],'LineStyle',':','LineWidth',0.01);
%plot(upper(3,:),'Color',[0.6 0 0.6],'LineStyle',':','LineWidth',0.8);
%plot(lower(3,:),'Color',[0.6 0 0.6],'LineStyle',':','LineWidth',0.8);

p5 = plot(means( 4 , length(truegdp)/2 + 1 : end ),'o','MarkerEdgeColor',[0 0 1],'MarkerSize',8);
plot(means(4, length(truegdp)/2 + 1 : end ),'Color',[0 0 1],'LineStyle',':','LineWidth',0.01);
%plot(upper(4,:),'Color',[0 0 1],'LineStyle',':','LineWidth',0.8);
%plot(lower(4,:),'Color',[0 0 1],'LineStyle',':','LineWidth',0.8);

p6 = plot(means(5 , length(truegdp)/2 + 1 : end ),'o','MarkerEdgeColor',[0 0.5 0.5 ],'MarkerSize',8);
plot(means(5, length(truegdp)/2 + 1 : end ),'Color',[0 0 0],'LineStyle',':','LineWidth',0.01);

p7 = plot(means(6 , length(truegdp)/2 + 1 : end ),'o','MarkerEdgeColor',[0 0 0],'MarkerSize',8);
plot(means(6, length(truegdp)/2 + 1 : end ),'Color',[0 0 0],'LineStyle',':','LineWidth',0.01);

xt = get(gca, 'XTick');
set(gca, 'FontSize', 7)

box off
grid on
%title(titlename,'Fontsize',10,'FontWeight','normal')
lgd = legend([p1 p2 p3 p4 p5 p6 p7], ...
       [flag_truegdp ' release'],...
       'Nd',...
       'NIG', ...
       'MG',...
       'PMNM',...
       'HS+',...
       'B-AR(1)',...
       'Location','best') ;
lgd.FontSize = 7;
ylim([min_global - range_global/10 max_global + range_global/10])
ylabel(ylabelname,'Fontsize',7)
xticks(1:1:length(truegdp)/2) 
xticklabels(xtickslabelname(length(xtickslabelname)/2+1:end))
end


