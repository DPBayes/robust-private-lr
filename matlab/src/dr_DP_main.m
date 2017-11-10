function lr_dp_main(params,info)
%%%
%%%This is the main entry for creating model on training data and give the prediction result for the tesing data.
%%%The data is loaded from preprocessed_all_data_random.mat which is generated by data_prepare*.mat
%%%This function mainly works as following:
%%%         1.load data from an existing file, get testing, internal, external and the noise data separately
%%%         2. l1normalize external data and then l2normalize the addition product of the external and noise data
%%%         3. provide the internal and normalized noised external data for training
%%%         4. use the trained model to predict the response of testing data
%%%         5. save the prediction result to result_filename
%%%Input:
%%%     params: struct, configuration parameters define by initial_params.m, please check the explanation in that file
%%%     rfile: string, result filename, which could also be defined in initial_params.model
%%%Output:
%%%     None
%%%
%%%
 
    %%%%load all_data
    load([params.odir params.processed_data_filename]); % for synthetic data     load('../data/all_data.mat');    
    %    load('preprocessed_all_data_random.mat');
 
    %{
    info.x_max=all_data.x_max;
    info.x_min=all_data.x_min;
    info.y_max=all_data.y_max;
    info.y_min=all_data.y_min;
    %}
    
    %%%%cross validation
   params.k_fold = 1;
   if(params.debug || params.sdebug) params.cross_valid=0; end
    if(params.cross_valid)
        params.k_fold = params.fold;
    end
    if(params.modell==2)
        params.l2normalize=0;
        params.l1normalize=0;
    end
    
     %%%%checking the effect of the size of private data
   %  if(params.sdebug) params.prsize=200; end
   %{
    list_prsize = [params.prsize];
    if(params.debug) 
       list_prsize = 100;
    else
        if(params.size_var)
           list_prsize = 0:100:params.prsize;
           %if(params.sdebug) list_prsize = 0:100:200; end
        end
    end
    %}
    
    if(params.no_internal)
        list_prsize = 100:100:params.prsize;
        sum_pr = sum(list_prsize);
        if(sum_pr<params.prsize)
            list_prsize = [list_prsize params.prsize];
        end
    else if(params.only_internal)
            list_prsize = 0;
        else
            list_prsize = 0:100:params.prsize;
            sum_pr = max(list_prsize);
           if(sum_pr<params.prsize)
            list_prsize = [list_prsize params.prsize];
           end
        end
    end
    prsize_len = length(list_prsize);
    
    %prsize_len = prsize_len + params.mult_size;
    
    
    %%%%privacy setting
    list_eps = 0;
    if(params.privacy)
    	%list_eps = [params.eps];  
        if(params.privacy_var)
        	list_eps = [0.2 1 10 100];
            list_eps = [0.0001 0.2 1 10 1000 100000];
            list_eps = [1];
        end
    else
        params.privacy_iter = 1;
    end

    Result = struct();

    %%%%cross validation
    Kindex=1:params.fold; %[4 2 1 5 3 ];
    for kindex = 1:params.k_fold
        kIndex=Kindex(kindex);
        params.extrapolate=0;
      %  disp(['fold ',num2str(kIndex)]);
       
        %train and test data preparation
	    tr = all_data.tr_data{kIndex};
        %internal data
        
        te = all_data.te_data{kIndex};
        
        if(params.l2normalize)
            tr.X = L2normalize(tr.X,params.l2norm_dim,params);
            te.X = L2normalize(te.X,params.l2norm_dim,params);
        end
        
         if(params.l1normalize)
            tr.X = L1normalize(tr.X,params.l1norm_dim,params);
            te.X = L1normalize(te.X,params.l1norm_dim,params);
         end
        
      
        
        %external(private data)
        pri_data = all_data.pri_data{kIndex};


        if params.corr_drug == 1
        %ranking by row during data_prepare
            tr.Y = tr.rY;
            te.Y = te.rY;
            pri_data.Y = pri_data.rY;
        else
        %ranking by column during data_prepare
            tr.Y = tr.cY;
            te.Y = te.cY;
            pri_data.Y = pri_data.cY;
        end

        %to check what can be learned when the external data only contain noise
        %{
        if params.noise_only == 1
            pri_data.X = pri_data.X*0;
        end
        %}
        %testing for the trade-off of utility and privacy
        eIndex = 1;%the index corresponds to the privacy budget
        for eps = list_eps
            params.eps = eps;
            % prData = struct();      prData.X = [];
                   
            if(params.privacy && params.noise_X && params.prsize > 0)
                    %extract noise data
                    for rIndex = 1:params.privacy_iter
                        %input perturbation
                        %disp('get input perturbation noise...');
                        %prData.X(rIndex,:,:) = all_data.noise_data{kIndex, eIndex, rIndex};
                    end
            end
             
            %the index corresponds to different private data size
           
            for sindex = 1:prsize_len;
                sIndex=sindex;
                 %disp(['private data extract ',int2str(sIndex)]);
                data = struct();
                %%%TODO:here may need change to sample
                %pri_index = all_data.pri_index{sIndex};
                
                
                if(sindex>numel(list_prsize)) sIndex=numel(list_prsize); params.extrapolate=params.extrapolate+1; end
                
                pri_index = 1:list_prsize(sIndex);
                 
                data.X = pri_data.X(pri_index, :);
                data.Y = pri_data.Y(pri_index, :);
                data.eps = params.eps;

                %currently, privacy_iter equals to 1. The origin aim of the for loop is to reduce the randomness of noise.
                %However, that randomness was tolerant during the pilot experiment.
                for rIndex = 1:params.privacy_iter
                    %{
                    %l1normalize is needed for input perturbation
                    if params.noise_only == 0 %(no need for OPS)
                                data.X = L1normalize(data.X);
                    end
                    if(params.noise_X)
                        %disp('private noised data generate...');
                        %when use the same epsilon value list
                        %data.X = data.X + squeeze(prData.X(rIndex, pri_index, :));
                        %when use different epsilon value from data_prepare
                        ipX = input_perturbation(data.X, params);
                        data.X=ipX;
                    end
                    %}
                    if(params.l2normalize)
                        data.X = L2normalize(data.X,params.l2norm_dim,params);
                    end
                      
                    if(params.l1normalize)
                        data.X = L1normalize(data.X,params.l1norm_dim,params);
                    end
                  %  keyboard;
                    
                    %[pred, ~] = train_prediction_test(tr, data, te, params);
                    Repeat=1;
                    Pred = zeros(size(te.Y));
                    for repeat=1:Repeat
                        pred = trpred(tr, data, te, params, info);
                        Pred = Pred + pred; 
                    end
                    Pred=Pred/Repeat;
                    
                    %save prediction result
                    Result.Y{kIndex, sindex, rIndex, eIndex} = Pred; %pred.Y; % rIndex and eIndex are fixed to one now
                end
            end
            eIndex = eIndex + 1;
        end
    end
    save(sprintf('%s/%s', params.odir, params.result_filename), 'Result','-v7.3');
   % disp('**finished**');   
end
