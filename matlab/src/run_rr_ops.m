function run_rr_ops()
%%%This file is to run multiple model on the same setting up, evaluate the
%%%prediction and also plot the result. Currently, it runs Ridge
%%%Regression, Ridge Regression with input perturbation, OPS with schema 1,
%%%and OPS with schema 2 orderly.
%%%PS: make sure the params.partition_data_dir and params.data_filename
%%%have beedn specified correctly.
%%%input:
%%%     None
%%%output:
%%%     None
%%%

    clear;
  
    Rpt=1; % repeats over seeding
    Sch=-1; % -1: plotting curves, 1: plotting tensor (use plot_tensor.m later)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(Sch>0)
        for sch=1:3
            initialize(sch-1);
        end
    end
    
    for rpt=1:Rpt
        
        
        if(Sch>0)
            rng(floor(rand*5000));
            disp(['rpt ',int2str(rpt)]);
            disp('......................................');
            disp('......................................');
            for sch=1:3
                %sch=0;
                disp(['sch ',int2str(sch)]);
                disp('...................');
                execute(sch-1);

            end
            close all;
        else
            rng(500);
            execute(-1);
        end
    end   
    
    if(Sch>0)
        for sch=1:3
            average_tensor(sch-1,Rpt);
        end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function average_tensor(sch,Rpt) 
         params = initial_params();
        if(params.corr)
            load([params.tensor_path,'sch',int2str(sch),'/tensor_corr.mat']);
        else
            load([params.tensor_path,'sch',int2str(sch),'/tensor_rmse.mat']);
        end
        
         T=T/Rpt;
        
       if(params.corr)
            save([params.tensor_path,'sch',int2str(sch),'/tensor_corr.mat'], 'T');
            %save([params.odir,'/tensor_corr.mat'], 'T');
       else
            save([params.tensor_path,'sch',int2str(sch),'/tensor_rmse.mat'], 'T');
        end
    end  
        
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function initialize(sch) 
         params = initial_params();%The detail information of the field in params, please, refers to the initial_params.m.
         
        if(sch==1)
            D=params.min_d:params.step_d:params.max_d;
        else
            D=10;
            D=65;
        end
        lenD=length(D);
        lenSch(2)=lenD;
        
       if(sch==2) 
            epss=1:5;
            epss=[0.5 epss];
       else
            epss=1;
       end
        lenE=length(epss);
        lenSch(3)=lenE;
      
        if(sch==3)
            C=0.01:0.01:0.1; %[0.1 0.05 0.01]; % 1e-2 5*1e-2 1e-3 5*1e-3 1e-4 5*1e-4 1e-5]; 0.9 0.8 0.7 0.6
        else
            C=0.05;
        end
        lenC=length(C);
        lenSch(4)=lenC;
                    
       

       if(sch==0)
            ins=params.min_insize:5:params.max_insize; %[5 10 15 20 25 30];
       else
           ins=10;
       end
       lenIn=length(ins);
        lenSch(1)=lenIn;
        
       
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
                T=zeros(params.no_models,lenSch(sch+1),prsize_len);
            else
                T=zeros(params.no_models,lenSch(1),prsize_len);
            end

            
            if(params.corr)
                save([params.tensor_path,'sch',int2str(sch),'/tensor_corr.mat'], 'T');
            else
                save([params.tensor_path,'sch',int2str(sch),'/tensor_rmse.mat'], 'T');
            end
    end

    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function execute(sch)
        comp=1; comp2=1; comp3=1;
        
        
        info.sch=sch;
        params = initial_params();%The detail information of the field in params, please, refers to the initial_params.m.
       if(sch==1)
            D=params.min_d:params.step_d:params.max_d;
       else
            D=10;
            %D=65;
       end
        lenD=length(D);
        lenSch(2)=lenD;
        
       if(sch==2) 
            epss=0.5:0.5:3;
            %epss=[0.5 epss];
       else
            epss=1;
            epss=2;	
       end
      lenE=length(epss);
      lenSch(3)=lenE;
      
        if(sch==3)
            C=0.01:0.01:0.1; %[0.1 0.05 0.01]; % 1e-2 5*1e-2 1e-3 5*1e-3 1e-4 5*1e-4 1e-5]; 0.9 0.8 0.7 0.6
        else
            C=0.05;
        end
        lenC=length(C);
        lenSch(4)=lenC;
                    
        
        if(params.only_plot) comp=0; comp2=0; end
       if(sch==0)
            ins=params.min_insize:5:params.max_insize; %[5 10 15 20 25 30];
       else
           ins=10;
       end
       lenIn=length(ins);
        lenSch(1)=lenIn;
        
       %%{
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
        
%%}
        
           params.insize=params.max_insize;
            params.prsize = params.total_size - params.tesize - params.insize; 


            
            info.maxd=D(lenD);
            info.mind=D(1);

       for j=1:params.fold     
            info.seed(j)=floor(rand*5000);
       end
        
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

                  %  params = initial_params();%The detail information of the field in params, please, refers to the initial_params.m.

                    prefix = '';
                    %when real need transfer to rank, the classify_test option should be set to 0.
                    if params.real2rank prefix = 'rank_'; end
                    if params.rmse prefix='rmse_'; end

                   

                   
                    for i=1:lenC
                        c=C(i);
                        info.c=c;
                        info.c_ind=i;
                        disp(['internal=',int2str(ins(in)),' dimension=',int2str(d),' c=',num2str(c),' eps=',num2str(eps)]);
                        
                        
                        if(data_gen) data_prepare(params,info); end
                        for iter=1:params.ITER    
                            %RR

                            %params = initial_params();%The detail information of the field in params, please, refers to the initial_params.m.

                            %RR OPT
                            if(1)
                                params_orr = params; 
                                params_orr.privacy = 0;
                                params_orr.privacy_var = 0;
                                params_orr.noise_X = 0;
                                params_orr.hs=0;
                                params_orr.dpvb=0;
                                params_orr.denoise=0;
                                params_orr.clip=0;
                                params_orr.modell=1;
                                params_orr.result_filename=[prefix 'OPT_RR.mat'];
                               % if(comp) dr_DP_main(params_orr,info); end
                                if(comp) lr_dp_main(params_orr,info); end
                            end
                            
                            %FM
                            if(0)
                                params_orr = params;
                                params_orr.privacy = 0;
                                params_orr.privacy_var = 0;
                                params_orr.noise_X = 0;
                                params_orr.dpvb=0;
                                params_orr.hs=0;
                                params_orr.clip=0;
                                params_orr.denoise=0;
                                params_orr.modell=2; %%
                                params_orr.result_filename=[prefix 'OPT_RR.mat'];
                                if(comp && params.no_models>1)
                                    dr_DP_main(params_orr,info);
                                end
                            end
                            
                            %OP
                            if(0)
                                params_orr = params;
                                params_orr.privacy = 0;
                                params_orr.privacy_var = 0;
                                params_orr.noise_X = 0;
                                params_orr.dpvb=0;
                                params_orr.hs=0;
                                params_orr.clip=0;
                                params_orr.denoise=0;
                                params_orr.modell=3; %%
                                params_orr.result_filename=[prefix 'OPT_RR.mat'];
                                if(comp && params.no_models>1)
                                    dr_DP_main(params_orr,info);
                                end
                            end
                            %disp('..............................................................');

                            %{
                            %for np clip

                            params_rr = params; 
                            params_rr.privacy = 0;
                            params_rr.privacy_var = 0;
                            params_rr.noise_X = 0;
                            params_rr.hs=0;
                            params_rr.dpvb=0;
                            params_rr.clip=1;
                            params_rr.model=1;
                            params_rr.result_filename=[prefix 'noDP_RR.mat'];
                            if(comp && params.no_models>1) dr_DP_main(params_rr,info); end
                            %disp('..............................................................');
                            %}

                           

                            %dp clip no HS
                            params_ops = params;
                            params_ops.privacy = 0;
                            params_ops.privacy_var = 0;
                            params_ops.noise_X = 0;
                            params_ops.dpvb=1;
                            params_ops.hs=1;
                            params_ops.denoise=0;
                            params_ops.clip=1;
                            params_ops.modell=1;
                            params_ops.result_filename =[prefix 'dpvb.mat'];
                            if(comp && params.no_models>2) 
                               % dr_DP_main(params_ops,info);
                                 lr_dp_main(params_ops,info);
                            end
                            %   disp('..............................................................');

                             % DP no clip no HS
                            params_rr = params;
                            params_rr.privacy = 0;
                            params_rr.privacy_var = 0;
                            params_rr.noise_X = 0;
                            params_rr.dpvb=1;
                            params_rr.hs=1;
                            params_rr.denoise=0;
                            params_rr.clip=0;
                            params_rr.modell=1;
                            params_rr.result_filename=[prefix 'noDP_RR.mat'];
                            if(comp && params.no_models>2) 
                                %dr_DP_main(params_rr,info);
                                lr_dp_main(params_rr,info);
                            end
                            %   disp('..............................................................');

                            %DP clip ver 2/ dpdenoise
                            %%{
                            if(0)
                                params_irr = params;
                                params_irr.privacy = 0;
                                params_irr.privacy_var = 0;
                                params_irr.noise_X = 0;
                                params_irr.dpvb=1;
                                params_irr.hs=1;
                                params_irr.clip=2;
                                params_irr.denoise=0;
                                params_irr.modell=1; %%
                                params_irr.result_filename=[prefix 'iDP_RR.mat'];
                                if(comp && params.no_models>1)
                                    dr_DP_main(params_irr,info);
                                end
                            end
                            %%}

                            %RR NP clip
                            %%{
                            if(1)
                                params_irr = params;
                                params_irr.privacy = 0;
                                params_irr.privacy_var = 0;
                                params_irr.noise_X = 0;
                                params_irr.dpvb=0;
                                params_irr.hs=0;
                                params_irr.clip=1;
                                params_irr.denoise=0;
                                params_irr.modell=1; %%
                                params_irr.result_filename=[prefix 'iDP_RR.mat'];
                                if(comp && params.no_models>1)
                                    %dr_DP_main(params_irr,info);
                                    lr_dp_main(params_irr,info);
                                end
                            end
                            %%}
                            
                            %OP perturbation
                            %%{
                            if(0)
                                params_irr = params;
                                params_irr.privacy = 0;
                                params_irr.privacy_var = 0;
                                params_irr.noise_X = 0;
                                params_irr.dpvb=0;
                                params_irr.hs=0;
                                params_irr.clip=0;
                                params_irr.denoise=0;
                                params_irr.modell=3; %%
                                params_irr.result_filename=[prefix 'iDP_RR.mat'];
                                if(comp && params.no_models>1)
                                    dr_DP_main(params_irr,info);
                                end
                                
                            end
                            
                             %for FM
                            if(0)
                                params_irr = params;
                                params_irr.privacy = 0;
                                params_irr.privacy_var = 0;
                                params_irr.noise_X = 0;
                                params_irr.dpvb=0;
                                params_irr.hs=0;
                                params_irr.clip=0;
                                params_irr.denoise=0;
                                params_irr.modell=2; %%
                                params_irr.result_filename=[prefix 'iDP_RR.mat'];
                                if(comp && params.no_models>1)
                                    dr_DP_main(params_irr,info);
                                end
                            end
                               %for Lasso
                            if(0)
                                params_irr = params;
                                params_irr.privacy = 0;
                                params_irr.privacy_var = 0;
                                params_irr.noise_X = 0;
                                params_irr.dpvb=0;
                                params_irr.hs=0;
                                params_irr.clip=0;
                                params_irr.denoise=0;
                                params_irr.modell=4; %%
                                params_irr.result_filename=[prefix 'iDP_RR.mat'];
                                if(comp && params.no_models>1)
                                    dr_DP_main(params_irr,info);
                                end
                            end
                            %%}

                            %{
                            %RR with DP clip and HS
                            params_irr = params;
                            params_irr.privacy = 0;
                            params_irr.privacy_var = 0;
                            params_irr.noise_X = 0;
                            params_irr.dpvb=1;
                            params_irr.hs=2;
                            params_irr.clip=1;

                            params_irr.model=1; %%
                            %}

                            if(params.no_models>4)
                                %for FM
                                %%{
                                params_fm = params;
                                params_fm.privacy = 0;
                                params_fm.privacy_var = 0;
                                params_fm.noise_X = 0;
                                params_fm.dpvb=0;
                                params_fm.hs=0;
                                params_fm.clip=0;

                                params_fm.modell=2; %%
                                %%}

                                params_fm.result_filename=[prefix 'iDP_FM.mat'];
                                if(comp && params.no_models>1)
                                    dr_DP_main(params_fm,info);
                                end
                                %    disp('..............................................................');
                            end

                            %%%%evaluate the result
                           if(comp2)
                                calculate_precision_recall(0, params);%evaluate according to the ranking of drugs
                            %calculate_precision_recall(1, params);%evaluate according to the ranking of patients
                           end

                           if(iter==1) action=0; else action=1; end
                           average_rmse(params,action,0);


                        end
                        action=2;
                        average_rmse(params,action,0)

                        %%%%plot
                        if(comp3 && params.no_models>1)
                            plot_precision_given_recall_ops(0, params,info);%plot the result for ranking over drugs
                        %plot_precision_given_recall_ops(1, params);%plot the reuslt for ranking over patients
                        end 
                    end
                end
            end
       end

    end
end

