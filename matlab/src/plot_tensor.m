function plot_tensor(info)
 hold all;
 h = figure;
 %positions = {[1:2, 4:5], 3, 6};
%positions = {[0.1, 0.16, 0.3, 0.8], [0.65, 0.66, 0.25, 0.3], [0.65, 0.16, 0.25, 0.3]};
%left bottom  width height
positions = {[0.7, 0.55, 0.25, 0.25], [0.1, 0.2, 0.45, 0.6], [0.7, 0.2, 0.25, 0.25]};
FONTSIZE = 8; %20;

%set(0,'defaulttextinterpreter','latex');

%{
ylabels = {'y 1', 'y 2', 'y 3'};


for k=1:2,
  %subplot(2, 3, positions{k});
  subplot('Position', positions{k});
  set(gca, 'FontSize', FONTSIZE);
  plot(randn(10, 1), randn(10, 1), '.');
  xlabel('x axis')
  ylabel(ylabels{k})
end

k = 3;
subplot('Position', positions{k});
set(gca, 'FontSize', FONTSIZE);
plot([5.0], [5.0], '.');
axis([0 1 0 1]);
axis off;
%legend('Points', 'Location', 'West')

set(gcf, 'Units', 'centimeters')
set(gcf, 'OuterPosition', [0, 0, 10.0, 10.0])
set(gcf, 'PaperUnits', 'centimeters')
set(gcf, 'PaperPosition', [0, 0, 10.0, 10.0])
set(gcf, 'PaperPositionMode', 'manual')

%}

 for sch_i=1:3
       sch=sch_i-1;
       
      subplot('Position', positions{sch_i});
      set(gca, 'FontSize', FONTSIZE);
      
      
      
      

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % set values 
        corr=1;
        params=initial_params();
        odir=info.odir; %'../data7';
        params.corr=corr;
        
        
        np_model=1;
        chk_model=2;
        base_model=1;
        in_col=1;

        ds=5:5:40; %66;
        lenD=length(ds);
        start_d_ind=1;
        end_d_ind=8;
        
        epss=0.5:0.5:3;
        %epss=[0.5 epss];
        start_eps_ind=2;
        lenE=length(epss);

        cs= 0.01:0.01:0.1; %[0.5 0.1 0.05 0.01];
        lenC=length(cs);


         ins=[0 5 10 15 20 25 30]; %params.min_insize:5:params.max_insize;%[10 15 20 25]; % 25 30];
         %ins=[5 10 20 30];
         lenIn=length(ins);
  
        if(sch==3)
            tensor_cvsext=1;
            tensor_evsext=0;
             tensor_invsext=0;
             tensor_dvsext=0;
            fav_d=10;
             fav_in=10;
            extd=6;
            strXl=cell(1,length(cs));
            %strXl{1,1}='{';
            %j=-4;

            for i=1:lenC
                %j=j+5;
                strXl{1,i}=num2str(cs(i)); %[''',int2str(i),''']; %[strXl,' ',int2str(i)];
            end
            strT=['Non-private data ',int2str(fav_in),', dimension ',int2str(fav_d),', epsilon ',int2str(1)];
            strX={'Clipping threshold'};

            save_fname=['../tensor_cvsext.eps'];
     
        end
        if(sch==2)
            tensor_evsext=1;
             tensor_invsext=0;
             tensor_dvsext=0;
             tensor_cvsext=0;
             fav_d=10;
            fav_in=10;
            extd=6;
            strXl=cell(1,length(epss));
            %strXl{1,1}='{';
            %j=-4;

            ii=0;
            for i=start_eps_ind:1:lenE
                %j=j+5;
                ii=ii+1;
                strXl{1,ii}=num2str(epss(i)); %[''',int2str(i),''']; %[strXl,' ',int2str(i)];
            end
            %strT=['Non-private data ',int2str(fav_in),', dimension ',int2str(fav_d),', clipping ',num2str(0.05)];
            strT='c';
            %strx2= text('interpreter','latex','fontsize',22,'string','\epsilon');           
            strX={'Privacy parameter\epsilon'};
            save_fname=['../tensor_evsext.eps'];
        end

        if(sch==1)
            tensor_cvsext=0;
            tensor_evsext=0;
            tensor_invsext=0;
            tensor_dvsext=1;

            extd=6;
            fav_in=10;
            strXl=cell(1,lenD);
            %strXl{1,1}='{';
            %j=-4;

            for i=1:1:length(ds)
                %j=j+5;
                strXl{1,i}=int2str(ds(i)); %[''',int2str(i),''']; %[strXl,' ',int2str(i)];
            end
            %strXl{1,j+1}='}';
            %strT=['Non-private data ',int2str(fav_in),', epsilon ',int2str(1),', clipping ',num2str(0.05)];
            strT='a';
            %strX={'Dimension'};
            strX={'Reduced dimensionality'};
            save_fname=['../tensor_dvsext.eps'];
        end

        if(sch==0)
            tensor_cvsext=0;
            tensor_evsext=0;
            tensor_dvsext=0;
            tensor_invsext=1;

            extd=5;
           fav_d=10;


            strXl=cell(1,lenIn);
            %strXl{1,1}='{';
            %j=-4;

            for i=1:lenIn
                %j=j+5;
                strXl{1,i}=int2str(ins(i)); %[''',int2str(i),''']; %[strXl,' ',int2str(i)];
            end
            %strT=['Dimension ',int2str(fav_d),', epsilon ',int2str(1),', clipping ',num2str(0.05)];
            strT='b';
            strX={'Size of non-private data'};
             save_fname=['../tensor_invsext.eps'];
        end


        params.insize=params.max_insize;
        params.prsize = params.total_size - params.tesize - params.insize; 
         extd=6;
            for i=1:extd
                if(i==1)
                    y=0; %params.insize;
                else
                    if(i<extd)
                        y=100*(i-1); %+params.insize;
                    else
                        y=params.prsize; %+params.insize;
                    end
                end
                strYl{1,i}=int2str(y);
            end

        strY={'Size of private dataset'};

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % initialization 
        miner=1e3;
        maxer=-1e3;

        %Z=zeros(len_ins,extd);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %compute limits and extract matrix for correct dimension, model and IN size
        for j=1:1
           load([odir,'sch1/tensor.mat']); 
           Z=squeeze(T(chk_model,:,:));
            %Z=squeeze(T(chk_model,start_d_ind:end_d_ind,:)); % for DP-clip model in range of dimensions for IN (size i) and all sizes of EXT
            % for NP model in range of dimensions for fav_in size of IN and all sizes of EXT 
            W=squeeze(T(base_model,:,:)); % for DP-clip model in range of dimensions for IN (size i) and all sizes of EXT

            %u=W(2,1)
			u=max(W(:,1));
			v=max(max(Z/u));
			v=v*10; v=floor(v); v=v/10;
			cmax=v
            %cmax=1.5
            load([odir,'sch',int2str(sch),'/tensor.mat']); 
            
            if(tensor_cvsext)

                    Z=squeeze(T(chk_model,1:lenC,:)); % for DP-clip model in range of dimensions for IN (size i) and all sizes of EXT
                    % for NP model in range of dimensions for fav_in size of IN and all sizes of EXT 
                    W=squeeze(T(base_model,1:lenC,:)); % for DP-clip model in range of dimensions for IN (size i) and all sizes of EXT

            end

            if(tensor_evsext)

                    Z=squeeze(T(chk_model,start_eps_ind:lenE,:)); % for DP-clip model in range of dimensions for IN (size i) and all sizes of EXT
                    % for NP model in range of dimensions for fav_in size of IN and all sizes of EXT 
                    W=squeeze(T(base_model,start_eps_ind:lenE,:)); % for DP-clip model in range of dimensions for IN (size i) and all sizes of EXT

            end

            if(tensor_dvsext)

                    Z=squeeze(T(chk_model,start_d_ind:end_d_ind,:)); % for DP-clip model in range of dimensions for IN (size i) and all sizes of EXT
                    % for NP model in range of dimensions for fav_in size of IN and all sizes of EXT 
                    W=squeeze(T(base_model,start_d_ind:end_d_ind,:)); % for DP-clip model in range of dimensions for IN (size i) and all sizes of EXT
            end

            if(tensor_invsext)
                Z=squeeze(T(chk_model,1:lenIn,:)); % for DP-clip model in range of dimensions for IN (size i) and all sizes of EXT
                %Y=squeeze(X(fav_d,1:extd)); % for DP-clip model for fav_d dimension and all sizes of EXT
                %Y=squeeze(X(fav_d,:)); % for DP-clip model for fav_d dimension and all sizes of EXT
                %if(length(Y)<extd) Y=[Y Y(end)]; end % to match with various sizes of EXT batch
                %Z(j,:)=Y;

                W=squeeze(T(base_model,1:lenIn,:)); % for DP-clip model in range of dimensions for IN (size i) and all sizes of EXT
                %Y=squeeze(X(fav_d,1:extd)); % for DP-clip model for fav_d dimension and all sizes of EXT
                %W(j,:)=Y;
            end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % setting tensors
        cmin=0;
        l=0; % value for normalization 
        winter=0;    

        thr=1; % value for binarizing,
        binarize=0;
        if(binarize) colormap(gray(2)); end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % plot tensors
        limo=[chk_model];
        for i=1:1
            j=limo(i);
            if(j==1) K=1; else K=1; end
            for k=1:K
                clear X;
                 
                 X=Z;
                 X=(X-l)/(u-l);
               
                 Y=W; Y=(Y-l)/(u-l);
          
                if(binarize) X(X<thr)=-1; X(X>thr)=1; end

                set(gca,...
                    'fontsize',FONTSIZE,...
                    'CLim',[cmin cmax],...
                    'XMinorTick', 'off','YMinorTick', 'off',...
                    'XGrid','off','YGrid','off')
                %caxis([cmin cmax]) 
                if(params.corr)
                    colorDepth = 1000; colormap(jet(colorDepth));
                else
                    
                    colorDepth = 1000; colormap(jet(colorDepth));
                    %colormap(flipud(jet(colorDepth)));
                end

               if(winter)  colormap('Winter'); end
               axis tight; %no white borders
               axis image; %real x,y scaling
              % colormap(mmap)
              %X(X==0)=NaN;
                X=padup_nan(X,0);
                %X(:,end)=[];
                pcolor(X');
                
               %grey = [0.8,0.8,0.8];
               %set(gca, 'Color', grey)
                if(tensor_dvsext || tensor_evsext || tensor_cvsext || tensor_invsext)
                   %set(gca,'XTick',[1,6,11,16,21,26,31]) %start_d:1:end_d)
                   %set(gca,'XTickLabel',{'5','10','15','20','25','30','35'}) %start_d:end_d)
                   set(gca,'XTick',[1:1:length(strXl)])
                   set(gca,'XTickLabel',strXl,'fontsize',FONTSIZE)
                end


               % title(strT); %,'fontweight','b')
                set(gca,'YTick',[1,2,3,4,5,6])
                if(tensor_invsext==0)
                    %set(gca,'YTickLabel',{'13';'113';'213';'313';'413';'506'}) %the last is not always 500
                     %set(gca,'YTickLabel',{'0';'100';'200';'300';'400';'493'}) %the last is not always 500
                     set(gca,'YTickLabel',strYl,'fontsize',FONTSIZE);
                else
                    %set(gca,'YTickLabel',{'0';'100';'200';'300';'400';'493'}) %the last is not always 500
                    set(gca,'YTickLabel',strYl,'fontsize',FONTSIZE);
                end


               xlabel(strX,'fontsize',FONTSIZE)
               ylabel(strY,'fontsize',FONTSIZE)
               shading flat
               shading interp
               


               %view([0.5 90]);


            end
        end

        %axes('Position', [0 0.11 .82 0.82], 'Visible', 'off');
        %axes('Position', [0.05 -0.05 0.9 0.9], 'Visible', 'off');

        % colorbar('SouthOutside') % add colorbar
        %print(h,'-dpdf','tensor_rmse.pdf')
        %set(gca,'XTickLabel',5:15); 
        %set(gca,'YTickLabel',0:6);

        caxis([cmin cmax]); 
        if(tensor_dvsext)
            if(params.corr)
                c=colorbar('location','EastOutside',...
                     'Ytick',[0,1,cmax]); %,...
                    %'YTickLabel',{'0 correlation','1 maximum NP','maximum DP relative to NP'});
                    c.Label.String='Relative improvement on rank correlation'; 
                    %ylabel(c,'L');
            else
                 c=colorbar('location','EastOutside',...
                    'Ytick',[cmin,1,cmax]); %,...
                c.Label.String='R^2 relative improvement';
            end
            c.Label.FontSize = FONTSIZE;
        end
        
       

        %{
        if(corr)


            
            %{
            fname=['../plots/tensor/tensor_corr_',int2str(sch),'.fig'];
            saveas(gcf,fname,'fig');
            fname=['../plots/tensor/tensor_corr_',int2str(sch),'.png'];
            saveas(gcf, fname, 'png');
            %}
        else
            %print(h,'-dpdf','../plots/tensor_rmse.pdf','-r0');
            print(h,'-depsc2','../plots/tensor_rmse.eps','-r0');
        end
    %}
 end
 % save_fname=['../plots/tensor/tensor_corr.eps'];
 if(params.corr)
    save_fname=['../tensor.eps'];
 else
     save_fname=[odir,'/tensor_r2.eps'];
 end
  set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 15.0, 10.0]);

  print(h,'-depsc2',save_fname,'-r0');
end
