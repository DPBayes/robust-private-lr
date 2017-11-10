function plot_precision_given_recall_ops(measure_by_column, params, info)
    
    FONTSIZE = 30;
    LINEWIDTH = 3;
    MARKERSIZE=25;
    mcolors=['m','b','r','g'];
    
    recall = [1 5 10];
    recall =[1 2];
    len_recall = length(recall);   
    if(params.rmse) len_recall=1;end
    
    if(info.sch==0) %in
        sdir=[params.tensor_path,'/sch0'];
        d=info.in_ind;
    end
    if(info.sch==1) %d
        sdir=[params.tensor_path,'/sch1'];
        d=info.ds_ind;
    end
    if(info.sch==2) %eps
        sdir=[params.tensor_path,'/sch2'];
        d=info.eps_ind;
    end
    if(info.sch==3) %c
        sdir=[params.tensor_path,'/sch3'];
        d=info.c_ind;
    end
    
    suffix = '_ranking_drugs';
    if measure_by_column == 1
        suffix = '_ranking_patients';
    end

     if params.rmse
        prefix='rmse_';
     end
    
    load([params.odir prefix 'SumResult.mat']);
    result=Sresult;
    
    if(info.sch>=0)
        if(params.corr)
            load([sdir,'/tensor_corr.mat']);
        else
            load([sdir,'/tensor_rmse.mat']);
        end
    end
    
    
    
    extd=6;
    for i=1:extd
        if(i==1)
            x=params.insize;
            
        else
            if(i<extd)
                x=100*(i-1)+params.insize;
            else
                x=params.prsize+params.insize;
            end
        end
       % strXl{1,i}=[int2str(xi),' + ',int2str(xe)];
        strXl{1,i}=int2str(x);
    end
    
  
    strX={'Size of data'};
    
    
    
    %sLen = size(result.rr, 1);
    [sLen, eLen, ~]= size(result.irr);
     co = [0 0 1;
      0 0.5 0;
      1 0 0;
      0 0.75 0.75;
      0.75 0 0.75;
      0.75 0.75 0;
      0.25 0.25 0.25];
    if(info.sch==-1)
        h = figure;
        box on;
        set(h, 'Position', [100, 100, 1500, 600]);
        %grid on;
       % set(gca,'ygrid','on')
        hold;
        set(gca,'FontSize',FONTSIZE);
        set(gca,'DefaultLineLineWidth',LINEWIDTH);
       
        
    end
    %set(groot,'defaultAxesColorOrder',co);
    for eIndex = 1:eLen
        for k = 1:len_recall
            subplot(eLen, len_recall, (eIndex-1)*len_recall+k);
             hold all;
            for j = 1:params.no_models %2 to ignore OPS, 3 with OPS
                [m, s] = concate_result(result, j, eIndex, recall(k), info);
                m
               % disp([int2str(j),' ',int2str(d)]);
                if(info.sch>=0)
                    x=squeeze(T(j,d,:));x=x';
                    x=x+[m zeros(1,info.width-length(m))];
                    T(j,d,:)=x;
                else

                    if(j==1)
                        plot(m(1)*ones(1,6),'--k'); %,'LineWidth',1);
                    end

                   if(params.cross_valid) 
                        if(j==1)
                            hold_m=m(1); hold_s=s(1);
                        end
                        if(j==4) % for params.modell==2
                            m(1)=hold_m;
                            s(1)=hold_s;
                        end
                        %subplot(1,2,1);
                        errorbar(1:sLen, m, s); 
                        %plot(m,'-.O','MarkerSize',MARKERSIZE,'color',mcolors(j),...
                           % 'MarkerFaceColor',mcolors(j)); 
                        
                   else
                       plot(m);
                   end
                   if j == params.no_models

                        if(params.synth_data)
           %                 set(gca,'XTickLabel', [10:100:(params.prsize+10)]);    
                        else
                         %   set(gca,'XTickLabel', [0,10,110,210,310,410,520]); %[10:100:(params.prsize+10)]);
                        set(gca,'XTick',[1:1:length(strXl)])
                        set(gca,'XTickLabel',strXl)
                        end
                        %%{
                        %xlabel(['Size of dataset (',int2str(params.insize),' non-private)'], 'fontsize', FONTSIZE); %,'fontweight','b');
                        %xlabel('Size of dataset (non-private + private)', 'fontsize', FONTSIZE); %,'fontweight','b');
                        xlabel('Size of dataset', 'fontsize', FONTSIZE); %,'fontweight','b');
                        if(params.corr)
                            %ylabel('Spearmann correlation coefficient, five-fold cross validation', 'fontsize', FONTSIZE); %,'fontweight','b');
                            ylabel('Correlation coefficient', 'fontsize', FONTSIZE); %,'fontweight','b');
                        else
                            ylabel('Mean squared error', 'fontsize',FONTSIZE);
                        end
                        %{
                        else
                            ylabel('precision', 'fontsize', 6);
                        end
                        %}
                        %hold all;
                  
                
                        if eIndex == 1 && k == 1
                     %   legend('OPT', 'DP-No-Clip', 'DP-Clip', 'DP-Denoise', 'location', 'northwest');
                      %  legend('OPT', 'DP-No-Clip', 'DP-Clip', 'DP-Clip-Denoise', 'location', 'northwest');
                        %legend('OPT', 'DP-No-Clip', 'DP-Clip', 'FM', 'location', 'northwest');
                       % legend('OPT', 'RR-Clip', 'DP-Clip', 'DP-Clip-HS', 'location', 'northwest');
                           if(info.sch==-1)
                             if(params.no_models==5)  
                                %legend('Baseline with non-private data','Ridge Regression', 'DP-No-Clip', 'DP-Clip', 'RR-Clip', 'Functional mechanism', 'location', 'northwest');
                                legend('Baseline with non-private data','Ridge Regression', 'DP-No-Clip', 'DP-Clip', 'Functional mechanism', 'location', 'northwest');
                             else
                                 
                                legend('Baseline with non-private data','Ridge Regression (RR)', 'Private RR', 'Robust private RR', 'Output perturbed RR', 'location', 'EastOutside');
                                %legend('Baseline with non-private data','Output perturbed RR', 'Private RR', 'Robust private RR', 'Functional mechanism RR', 'location', 'EastOutside');
                             end
                           end
            %            legend('NP', 'DP-fhs1', 'DP-vhs'); %, 'location', 'southwest');
                        %legend('RR', 'iRR', 'location', 'northwest');
                            if(params.corr)
                               % title('Ranking of patients averaged over 124 drugs, dimension 10, \epsilon=1','fontsize',FONTSIZE); %,'fontweight','b');
                            else
                               % title('MSE averaged over 124 drugs, dimension 10, \epsilon=1','fontsize',FONTSIZE); %,'fontweight','b');
                            end
                        end

                      hold off;
                   end
                    %{
                    if(params.rmse)
                        if measure_by_column == 1
                            title('RMSE over patients averaged across drugs');
                        else
                            title('RMSE over drugs averaged across patients');
                        end
                    end

                    if recall(k) == 5
                        if measure_by_column == 1
                            title('ranking patients');
                        else
                            title('ranking drugs');
                        end
                    end
                    %}
                end
            end
        end
    end
    xlim([1 6]);
   % ylim([7.7*1e-3 8.1*1e-3]);
    %if(params.rmse==0) ylim([0 1.1]); end
    
    %saveas(h, [result_dir 'plot_rr_ops1_ops2_mean' suffix '.png'], 'png');
    %c=info.c*1e2;
    
    %print(h,'-dpdf',[params.plotdir,'/',int2str(d),'_',int2str(c),'.pdf'])
   
    if(info.sch==-1)
        if(params.corr)
            fname=[params.plotdir,'/corr.eps'];
            print(h,'-depsc2',fname,'-r0');
            fname=[params.plotdir,'/corr.fig'];
            saveas(gcf,fname,'fig');
        else
            fname=[params.plotdir,'/rmse.eps'];
            print(h,'-depsc2',fname,'-r0');
            fname=[params.plotdir,'/rmse.fig'];
            saveas(gcf,fname,'fig');
        end
    end    
    
    %system(['epstopdf ',fname]); 
    %saveas(gcf,[params.plotdir,'/',int2str(d),'_',int2str(c),'.fig'],'fig')
    %close all;
      
   
    
    if(info.sch>=0)
    if(params.corr)
        save([sdir,'/tensor_corr.mat'], 'T');
       % save([params.odir,'/tensor_corr.mat'], 'T');
    else
        save([sdir,'/tensor_rmse.mat'], 'T');
    end
    end
    %{
    if eLen > 1
  %      saveas(h, [params.result_dir 'plot_rr_ops1_ops2' suffix '_multipleEpsilon.png'], 'png');
    else
   %     saveas(h, [params.result_dir 'plot_rr_ops1_ops2' suffix '.png'], 'png');
    end
    %}
end

function [m, s] = concate_result(result, index, eIndex, k, info)
    stand = NaN;

    if index == 1
        stand = result.orr;
    elseif index == 2
        stand = result.rr;
    elseif index == 3
        stand = result.dpvb;
    elseif index == 4
        stand = result.irr;
    elseif index == 5
        stand = result.fm;
    end
    sLen = size(stand, 1);
    kFold = size(stand{1, 1}, 1);
    if index == 1
        eIndex = 1;
    end
    res = NaN(kFold, sLen);
    for i = 1:sLen
        res(:, i) = stand{i, eIndex}(:, k);
    end
    
    m = nanmean(res);
    %n=(sqrt(info.prsize+info.insize));
    s = std(res); %./n;
end

