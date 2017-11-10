function run_rplr_exp()
    clear;
    Rpt=1; % repeats over Monte Carlo crossvalidation.
    %Set to 1 for quick check. Value used is 100 for plots on the paper.
    
    % geneate results on correlation co-efficients
    Sch=-1; % -1: plotting curves, 1: plotting tensor (use plot_tensor.m later)
    info.odir='../codir/';
    info.cindex=0;
    experiments(Sch,Rpt,info);
    
    %generate results on wpc-index
    Sch=-1;
    info.odir='../wodir/';
    info.cindex=1;
    experiments(Sch,Rpt,info);
    
    %generate tensors
    %Sch=1; % -1: plotting curves, 1: plotting tensor (use plot_tensor.m later)
    %info.odir='../todir/';
    %info.cindex=0;
    %experiments(Sch,Rpt,info);
end

function experiments(Sch,Rpt,info)
    %-----------------------------------------------------------------------
    % initialize 
    if(Sch>0)
        for sch=1:3
            initialize(sch-1,info);
        end
    end
    
    %-----------------------------------------------------------------------
    % execute 
    for rpt=1:Rpt
        % for tensor
        if(Sch>0)
            rng(floor(rand*5000));
            disp(['rpt ',int2str(rpt)]);
            disp('......................................');
            disp('......................................');
            for sch=1:3
               
                disp(['sch ',int2str(sch)]);
                disp('...................');
                execute(sch-1,info);

            end
            close all;
           
        else
            
            rng(500);
            execute(-1,info);
            %mplot(info)
        end
    end   
    
    if(Sch>0)
        for sch=1:3
            average_tensor(sch-1,Rpt,info);
        end
         plot_tensor(info)
    end
    
    
    %-----------------------------------------------------------------------
    function average_tensor(sch,Rpt,info) 
         params = initial_params();
         load([info.odir,'sch',int2str(sch),'/tensor.mat']);
         %{
        if(params.corr)
            load([params.tensor_path,'sch',int2str(sch),'/tensor_corr.mat']);
        else
            load([params.tensor_path,'sch',int2str(sch),'/tensor_rmse.mat']);
        end
        %}
         T=T/Rpt;
        save([info.odir,'sch',int2str(sch),'/tensor.mat'], 'T');
        %{
       if(params.corr)
            save([params.tensor_path,'sch',int2str(sch),'/tensor_corr.mat'], 'T');
      else
            save([params.tensor_path,'sch',int2str(sch),'/tensor_rmse.mat'], 'T');
        end
        %}
    end  
        
    
    %-----------------------------------------------------------------------
    %<<<<<<<<<<
    function initialize(sch,info) 
     
         params = initial_params();
         %The detail information of the field in params, please, refers to the initial_params.m.
         
        % vary dimension 
        if(sch==1)
            D=params.min_d:params.step_d:params.max_d;
        else
            D=10;
           % D=65;
        end
        lenD=length(D);
        lenSch(2)=lenD;
        
       %vary epsilon 
       if(sch==2) 
            epss=1:5;
            epss=[0.5 epss];
       else
            epss=2; %%% NB: only used in tensor
       end
        lenE=length(epss);
        lenSch(3)=lenE;
      
        %vary clipping threshold (not used anymore)
        if(sch==3)
            C=0.01:0.01:0.1; 
        else
            C=0.05;
        end
        lenC=length(C);
        lenSch(4)=lenC;
                    
       
       % vary internal data size 
       if(sch==0)
            ins=params.min_insize:5:params.max_insize; %[5 10 15 20 25 30];
       else
           ins=10;
       end
       lenIn=length(ins);
       lenSch(1)=lenIn;
        
       %-----------------------------------------------------------------------
       %<<<<<<<<<<
       % define the size of T to hold the tensor
        params.insize=params.min_insize;
        params.prsize = params.total_size - params.tesize - params.max_insize; 

        % to correctly set size of T
        if(params.no_internal)
            list_prsize = 100:100:params.prsize;
            sum_pr = list_prsize(length(list_prsize));
            if(sum_pr<params.prsize)
                list_prsize = [list_prsize params.prsize];
            end
        else if(params.only_internal)
                list_prsize = 0;
            else
                list_prsize = 0:100:params.prsize;
                sum_pr = list_prsize(length(list_prsize));
               if(sum_pr<params.prsize)
                list_prsize = [list_prsize params.prsize];
               end
            end
        end
        prsize_len = length(list_prsize);
        info.width=prsize_len;
        if(sch>0)
            T=zeros(2,lenSch(sch+1),prsize_len);
        else
            T=zeros(params.no_models,lenSch(1),prsize_len);
        end

        save([info.odir,'sch',int2str(sch),'/tensor.mat'], 'T');
        %{
        if(params.corr)
            save([params.tensor_path,'sch',int2str(sch),'/tensor_corr.mat'], 'T');
        else
            save([params.tensor_path,'sch',int2str(sch),'/tensor_rmse.mat'], 'T');
        end
        %}
        %>>>>>>>>>>>>>
        
    end
    %>>>>>>>>>>>>> initialize 
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % sch:= -1:baselines, 0:vary internal data, 1:vary dimension, 2:vary epsilon 
    function execute(sch,info)
        params = initial_params();
       
        comp=1; comp2=1; comp3=1;
        if(params.only_plot) comp=0; comp2=0; end
       
        info.sch=sch;
        
        %-----------------------------------------------------------------------
        %<<<<<<<<<< configure parameters (same as in initialize)
        %vary internal data
        if(sch==0)
            ins=params.min_insize:5:params.max_insize; 
        else
           ins=10;
        end
        lenIn=length(ins);
        lenSch(1)=lenIn;
        
       

        %vary dimension
        if(sch==1) 
            D=params.min_d:params.step_d:params.max_d;
        else
            % two different sets of experiments done
            D=10; 
            %D=65;
        end
        lenD=length(D);lenSch(2)=lenD;
        info.maxd=D(lenD);
        info.mind=D(1);
        
        if(sch==2) 
            epss=0.5:0.5:3;
        else
            % two different set of experiments are done
            %epss=1; 
            epss=2;	 %%% NB: used in comparison plots, run both cases epss=1 and epss=2
        end
        lenE=length(epss);
        lenSch(3)=lenE;
      
        %not used anymore
        if(sch==3)
           C=0.01:0.01:0.1; 
        else
           C=0.05;
        end
        lenC=length(C);
        lenSch(4)=lenC;
        
        for j=1:params.fold     
            info.seed(j)=floor(rand*5000);
        end
        %>>>>>>>>>>>>> configure parameters              
        
       
        
        
            %not used
            params.insize=params.insize;
            params.prsize = params.total_size - params.tesize - params.max_insize; 

            % to correctly set size of T
            if(params.no_internal)
                list_prsize = 100:100:params.prsize;
                sum_pr = list_prsize(length(list_prsize));
                if(sum_pr<params.prsize)
                    list_prsize = [list_prsize params.prsize];
                end
            else if(params.only_internal)
                    list_prsize = 0;
                else
                    list_prsize = 0:100:params.prsize;
                    sum_pr = list_prsize(length(list_prsize));
                   if(sum_pr<params.prsize)
                    list_prsize = [list_prsize params.prsize];
                   end
                end
            end
            prsize_len = length(list_prsize);
            info.width=prsize_len;
            
            %{
            if(sch>0)
                T=zeros(params.no_models,lenSch(sch+1),prsize_len);
            else
                T=zeros(params.no_models,lenSch(1),prsize_len);
            end
            if(params.corr)
                save([params.odir,'tensor_corr.mat'], 'T');
            else
                save([params.odir,'tensor_rmse.mat'], 'T');
            end
            %}
        
         
         %to set size of external data 
        params.insize=params.max_insize;
        params.prsize = params.total_size - params.tesize - params.insize; 

        
            
           
     %-----------------------------------------------------------------------
     %<<<<<<<<<< start execution
     for in=1:lenIn
            info.in_ind=in;
            info.in=ins(in);
            params.insize=ins(in);
             
            info.prsize=list_prsize;
            info.insize=params.insize;
           % params.prsize = params.total_size - params.tesize - params.insize; 

            for j=1:lenD

                d=D(j);
                data_gen=1;
                info.ds_ind=j;
                info.d=d;

                for e=1:lenE
                    eps=epss(e);
                    info.eps=eps;
                    info.eps_ind=e;
                    
                    for i=1:lenC
                        c=C(i);
                        info.c=c;
                        info.c_ind=i;
                        
                        
                        disp(['internal=',int2str(ins(in)),' dimension=',int2str(d),' eps=',num2str(eps)]);
                        
                        
                        if(data_gen) data_prepare(params,info); end
                        
                        for iter=1:params.ITER    
                           
                            % non-private LR
                            if(1)
                                params_lr = params; 
                                params_lr.privacy = 0;
                                params_lr.privacy_var = 0;
                                params_lr.noise_X = 0;
                                params_lr.hs=0;
                                params_lr.rplr=0;
                                params_lr.denoise=0;
                                params_lr.clip=0;
                                params_lr.modell=1;
                                params_lr.result_filename=['lr.mat'];
                                if(comp) lr_dp_main(params_lr,info); end
                                disp(['LR done'])
                            end
                            
                          

                            %rplr
                            if(1)
                                params_rplr = params;
                                params_rplr.privacy = 0;
                                params_rplr.privacy_var = 0;
                                params_rplr.noise_X = 0;
                                params_rplr.rplr=1;
                                %params_rplr.hs=1;
                                %params_rplr.denoise=0;
                                params_rplr.clip=1;
                                params_rplr.modell=1;
                                params_rplr.result_filename =['rplr.mat'];
                                if(comp && params.no_models>2) 
                                    lr_dp_main(params_rplr,info);
                                end
                                   disp('RPLR done');
                            end
                            
                            if(sch<0)
                                % private LR
                                params_plr = params;
                                params_plr.privacy = 0;
                                params_plr.privacy_var = 0;
                                params_plr.noise_X = 0;
                                params_plr.rplr=1;
                                params_plr.hs=1;
                                params_plr.denoise=0;
                                params_plr.clip=0;
                                params_plr.modell=1;
                                params_plr.result_filename=['plr.mat'];
                                if(comp && params.no_models>2) 
                                    lr_dp_main(params_plr,info);
                                end
                                disp('PLR done')
                            end
                            
                            %robust non-private  LR
                            %%{
                            if(sch<0)
                                params_rnplr = params;
                                params_rnplr.privacy = 0;
                                params_rnplr.privacy_var = 0;
                                params_rnplr.noise_X = 0;
                                params_rnplr.rplr=0;
                                params_rnplr.hs=0;
                                params_rnplr.clip=1;
                                params_rnplr.denoise=0;
                                params_rnplr.modell=1; %%
                                params_rnplr.result_filename=['rnplr.mat'];
                                if(comp && params.no_models>1)
                                    lr_dp_main(params_rnplr,info);
                                end
                                disp('RNPLR done')
                            end
                            %%}
                            
                            %OP perturbation
                            %%{
                            if(sch<0)
                                params_oplr = params;
                                params_oplr.privacy = 0;
                                params_oplr.privacy_var = 0;
                                params_oplr.noise_X = 0;
                                params_oplr.rplr=0;
                                params_oplr.hs=0;
                                params_oplr.clip=0;
                                params_oplr.denoise=0;
                                params_oplr.modell=3; %%
                                params_oplr.result_filename=['oplr.mat'];
                                if(comp && params.no_models>1)
                                    lr_dp_main(params_oplr,info);
                                end
                                
                            end
                            
                             %for FM
                            if(sch<0)
                                params_fmlr = params;
                                params_fmlr.privacy = 0;
                                params_fmlr.privacy_var = 0;
                                params_fmlr.noise_X = 0;
                                params_fmlr.rplr=0;
                                params_fmlr.hs=0;
                                params_fmlr.clip=0;
                                params_fmlr.denoise=0;
                                params_fmlr.modell=2; %%
                                params_fmlr.result_filename=['fmlr.mat'];
                                if(comp && params.no_models>1)
                                    lr_dp_main(params_fmlr,info);
                                end
                            end
                            
                            
                          
                            %%%%evaluate the result
                           if(comp2)
                               compute_error(params,info);
                               disp(['Error computed'])
                           end

                         
                           
                          % if(iter==1) action=0; else action=1; end
                           %average_rmse(params,action,0);


                        end
                        %action=2;
                        %average_rmse(params,action,0)

                        if(sch>=0)
                               save_tensor(params,info)
                        end
                    end
                end
            end
       end

    end
end

