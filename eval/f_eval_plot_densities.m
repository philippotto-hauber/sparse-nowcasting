function f_eval_plot_densities(flag_survey, flag_sample, flag_truegdp, Np, flag_country) 

    % - options ---
    % ----------------------------
    evaloptions = load_evaloptions(flag_country) ; 

    % - directories ---
    % ----------------------------
    dir_load = ['C:\Users\Philipp\Documents\Dissertation\sparse nowcasting\eval\' flag_country '\' flag_survey ' ' flag_sample '\Np = ' num2str(Np) '\' flag_truegdp '\'] ;
    dir_save = [dir_load 'graphs\'] ; 
    if exist(dir_save, 'dir') ~= 7;mkdir(dir_save); end  
    dir_truegdp = 'C:\Users\Philipp\Documents\Dissertation\sparse nowcasting\data\out\' ;

    % - load forecast structure
    % -----------------
    load([dir_load 'results_eval.mat'])

    % - load true gdp mat-file ---
    % ----------------------------
    load([dir_truegdp 'truegdp' flag_country '.mat'])
    
    % - filenames and 
    % --------------------------
    filename_ps = [dir_save 'graphs_densities.ps'] ; 
    filename_pdf = [dir_save 'graphs_densities.pdf'] ; 
    
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
        
        for q = 1 : evaloptions.Nquarters
            % true gdps and benchmark and string to save   
            if strcmp(flag_truegdp,'first')
                truegdp = evaloptions.multfac*truegdp_strct.first(q) ; 
            elseif strcmp(flag_truegdp,'second')            
                truegdp = evaloptions.multfac*truegdp_strct.second(q) ; 
            elseif strcmp(flag_truegdp,'final')            
                truegdp = evaloptions.multfac*truegdp_strct.final(q) ; 
            end


            % open figure
            figure
            fig = gcf;
            fig.PaperOrientation = 'landscape';

            for h = 1 : evaloptions.Nhs

                % name for main title of plot
                maintitlename = [ results_eval.priors(1).R(r).horizon(h).name ] ; 

                % get benchmark densities
                densBAR = results_eval.benchmark_BAR.horizon( h ).dens( q , 1 : end  ) ;

                % get prior densities
                if index_r == length( evaloptions.Nrs ) + 1
                    densNIG = results_eval.priors( 1 ).pool.horizon( h ).dens{ q } ; 
                    densMG = results_eval.priors( 2 ).pool.horizon( h ).dens{ q } ;
                    densPMNM = results_eval.priors( 3 ).pool.horizon( h ).dens{ q } ;
                    densHS = results_eval.priors( 4 ).pool.horizon( h ).dens{ q } ;
                    densNd = results_eval.priors( 5 ).pool.horizon( h ).dens{ q } ;
                else
                    densNIG = results_eval.priors( 1 ).R( r ).horizon( h ).dens{ q } ; 
                    densMG = results_eval.priors( 2 ).R( r ).horizon( h ).dens{ q } ;
                    densPMNM = results_eval.priors( 3 ).R( r ).horizon( h ).dens{ q } ;
                    densHS = results_eval.priors( 4 ).R( r ).horizon( h ).dens{ q } ;
                    densNd = results_eval.priors( 5 ).R( r ).horizon( h ).dens{ q } ;
                end

                % plot
                subplot(2,2,h)
                f_plotdens(densNIG,densMG,densPMNM,densHS,densNd,densBAR,truegdp,[flag_truegdp ' release'],maintitlename)
            end
            % save
            if index_r == length( evaloptions.Nrs ) + 1
                dim = [0 0.9 1 0.1] ;
                str = [name_country_long ', ' truegdp_strct.quarters{q} ', pool' ];
                annotation('textbox',dim,...
                    'String',str,...
                    'EdgeColor', 'none', ...
                    'HorizontalAlignment', 'center',...
                    'FontSize',12,...
                    'FontWeight','bold');
            else
                dim = [0 0.9 1 0.1] ;
                str = [name_country_long ', ' truegdp_strct.quarters{q} ', R=' num2str(r) ];
                annotation('textbox',dim,...
                    'String',str,...
                    'EdgeColor', 'none', ...
                    'HorizontalAlignment', 'center',...
                    'FontSize',12,...
                    'FontWeight','bold');
            end
            
            print(filename_ps,'-dpsc','-fillpage','-append')
            close
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

function f_plotdens(densNIG,densMG,densPMNM,densHS,densNd,densBAR,truegdp,nametruegdp ,titlename)
% figure
% fig = gcf;
% fig.PaperOrientation = 'landscape';
n = 2^6 ; 
ylim( [ 0 1 ] )
max_global = max( [ max( densNIG ) , max( densMG ) , max( densPMNM ) , max( densHS ) , max( densNd ) , max( densBAR ) ] ) ; 
min_global = min( [ min( densNIG ) , min( densMG ) , min( densPMNM ) , min( densHS ) , min( densNd ) , min( densBAR ) ] ) ; 
range_global = max_global - min_global ; 

[ ~ , yNd , xs , ~ ] = kde( densNd , n , min_global - range_global/10 , max_global + range_global/10 ) ;
plot(xs,yNd,'Color',[0.9 0 0],'LineWidth',1.2)
hold on

[ ~ , yNIG , ~ , ~ ] = kde( densNIG , n , min_global - range_global/10 , max_global + range_global/10 ) ;
plot(xs,yNIG,'Color',[1 .5 0],'LineWidth',1.2)
hold on
[ ~ , yMG , ~ , ~ ] = kde( densMG , n , min_global - range_global/10 , max_global + range_global/10 ) ;
plot(xs,yMG,'Color',[0.6 0 0.6],'LineWidth',1.2)

[ ~ , yPMNM , ~ , ~ ] = kde( densPMNM , n , min_global - range_global/10 , max_global + range_global/10 ) ;
plot(xs,yPMNM,'Color',[0 0 1],'LineWidth',1.2)

[ ~ , yHS , ~ , ~ ] = kde( densHS , n , min_global - range_global/10 , max_global + range_global/10 ) ;
plot(xs,yHS,'Color',[0 0.5 0.5],'LineWidth',1.2)

[ ~ , yBAR , ~ , ~ ] = kde( densBAR , n , min_global - range_global/10 , max_global + range_global/10 ) ;
plot(xs,yBAR,'Color',[0 0 0],'LineWidth',1.2,'LineStyle',':')

line([truegdp truegdp],ylim,'Color',[0 0 0],'LineStyle','-','LineWidth',1.2);

title(titlename,'Fontsize',10,'FontWeight','normal')
lgd = legend(['Nd:\sigma=' num2str(round(std(densNd),2))],...
        ['NIG:\sigma=' num2str(round(std(densNIG),2))], ...
       ['MG:\sigma=' num2str(round(std(densMG),2))], ...
       ['PMNM:\sigma=' num2str(round(std(densPMNM),2))],...
       ['HS:\sigma=' num2str(round(std(densHS),2))],...       
       ['B-AR(1):\sigma=' num2str(round(std(densBAR),2))],...
        nametruegdp,...
        'Location','best') ;
lgd.FontSize = 5;    
end






