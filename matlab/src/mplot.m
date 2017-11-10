function mplot(info)
    %odir='../data'; 
    trim=0;
    EB=1; %errorbar
    c_index=info.cindex;
    params=initial_params;
    odir=info.odir;
    
    params.no_models=6;
    FONTSIZE = 6; %30;
    LINEWIDTH = 1.2;
    MARKERSIZE=4;
    loyolagreen = 1/255*[0,104,87];
    purple=[0.2857         0    1.0000];
    mcolors={'m','b','r','g','c',loyolagreen};
    
    
    params.insize=10;
   
    %set the xTickLabels
    extd=6;
    for i=1:extd
        if(i==1)
            xi=params.insize;
            xi=10;
            xe=0;
        else
            if(i<extd)
                xe=100*(i-1);
            else
                xe=params.prsize;
            end
        end
        strXl{1,i}=[int2str(xi),' + ',int2str(xe)];
       
    end
    

    h = figure;
    box on;
    cc=hsv(12);
    hold on;
    for j = 1:params.no_models
        [m, s] = get_result(j, odir);
       % m
        if(j==1)
            plot(m(1)*ones(1,6),'--s','color','k','MarkerSize',MARKERSIZE,'MarkerFaceColor','k'); %,'LineWidth',1);
            if(c_index==0)
                plot(0.0752*ones(1,6),'--s','color',cc(10,:),'MarkerSize',MARKERSIZE,'MarkerFaceColor',cc(10,:)); %,'LineWidth',1);
            else
                plot(0.51*ones(1,6),'--s','color',cc(10,:),'MarkerSize',MARKERSIZE,'MarkerFaceColor',cc(10,:)); %,'LineWidth',1);
            end
        end
        if(j==1)
            hold_m=m(1); hold_s=s(1);
        end
        if(j==4 || j==5 || j==6) % for params.modell==2
            m(1)=hold_m;
            s(1)=hold_s;
        end
        if(trim && j==5)
            s=zeros(size(s));
        end
        cindex=j;

        if(0 && (cindex==1 || cindex==4))
                errorbar(1:extd, m, s,...
                '--s','MarkerSize',MARKERSIZE,'color',mcolors{cindex},...
                            'MarkerFaceColor',mcolors{cindex}); %,'LineWidth',LINEWIDTH); 
         else
                    errorbar(1:extd, m, s,...
                '-.s','MarkerSize',MARKERSIZE,'color',mcolors{cindex},...
                            'MarkerFaceColor',mcolors{cindex},'LineWidth',LINEWIDTH); 
         end

      
        xlabel('Size of dataset (internal+external)', 'fontsize', FONTSIZE); %,'fontweight','b');
            if(c_index)
                ylabel('Wpc-index', 'fontsize', FONTSIZE); %,'fontweight','b');
            else    
                ylabel('Correlation coefficient', 'fontsize', FONTSIZE); %,'fontweight','b');
            end
        
    end
    
    legend('LR non-private data (dimension=10)',...
            'LR non-private data (dimension=65, best)',...
            'Linear Regression (LR)', 'Private LR', 'Robust private LR',...
            'Output perturbed LR',...
            'Functional Mechanism LR',...
            'Robust LR',...
            'location', 'Best'); 

    hold off;
    xlim([0.5 6.5]);
    set(gca,'XTick',1:6);
    set(gca,'XTickLabel',strXl,'fontsize',FONTSIZE);
    if(trim)
        if(c_index==0)
            ylim([-0.5 0.5]);
        end
    else
        if(c_index==0)
            ylim([-0.03 0.22]);
        else
            ylim([0.49 0.55]);
        end
    end
    
    set(gca,'FontSize',FONTSIZE);
    set(gca,'DefaultLineLineWidth',LINEWIDTH);
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 15.0, 12.0]);
  
    if(c_index)
        save_fname='../wpc.eps';
    else
        save_fname='../corr.eps';
    end
    print(h,'-depsc2',save_fname,'-r0');
   
end

function [m, s] = get_result(index,  odir)
    load([odir '/result.mat']);
    stand = NaN;
    if index == 1
        stand = result.lr;
    elseif index == 2
        stand = result.plr;
    elseif index == 3
        stand = result.rplr;
    elseif index == 6
        stand = result.rnplr;
    elseif index == 4
        stand = result.oplr;
    elseif index == 5
        stand = result.fmlr;
    end
    sLen = size(stand, 1);
    kFold = size(stand, 2);
   
    
    res = NaN(kFold, sLen);
    for i = 1:sLen
        res(:, i) = [stand{i,:}];
    end
    m = nanmean(res);
   s = std(res); 
 end

