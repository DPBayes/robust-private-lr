function data_prepare(params, info)

    X = [];
    Y = [];
    %%%%read data
    if(params.synth_data)
            synth_data=gen_synth_data(params);
            save([params.odir params.data_filename], 'synth_data');
        X = synth_data.X; %.value;
        Y = synth_data.Y;
    else
       % GE = 4; CN = 3; MT = 2; TT = 1;
        
        % ver1: number of data points n = 650, number of drugs m = 116 
        % ver2: number of data points n = 613, number of drugs m = 124 
        
        %ver1 Y: 650 x 116 n x m
        %ver2 Y: 613 x 124 n x m
        %ver3 Y: 985 x 265 n x m
        load('../Data/ver3/DrugResponse.mat');
        Y=DrugResponse';
        Y=Y'; %for ver2, ver3
        
        % ver1 X: 650 x 13321 n x d
        % ver2 X: 613 x 13321 n x d
        % ver3 X: 985 x 17490 n x d
        load('../Data/ver3/GeneExpression.mat');
        X=GeneExpression; 
        
    end
    
    %%%% handling missing value
    if(params.mod_nan)
        % replace NaN with 0
        X=padup_nan(X,params.pad);
        Y=padup_nan(Y,params.pad);
    end
    
    
    %%%%log transform and rescaling
    %log transform and scale
    if(params.logscaling)
        disp('logscaling...');
        X = log(1+X);
        %X = rescaling(X);
        %Y = rescaling(Y);
    end

    %%%%dimension reduction
    if(params.dimension_reduce)
        if(0)
          feature_sets = feature_selection(X, Y, params);
        end
        
        %no sorting 
        %%{
        if(0)
            load('../Data/ver2/UsefulGenesIndexes.mat');
            feats = UsefulGenesIndexes;
        end        
        feats=findGeneIndex;
        
        %%}
        
        %lasso sorted 
        %{
        load('ver2/lasso_sorted_genes.mat');
        feats=sorted_genes;
        %}
        
        %feats = [feature_sets(1:10)' feature_sets(13:14)']';  
        %feats = [feature_sets(1:10)' feature_sets(15:16)']';  
         d=info.d;
         if(d>0) 
            X = X(:, feats);
            X = X(:,1:d);
          end
    end

    %{
    % to clamp response in -1,1 or 0,1
    if(0)
        Yr=clamp_data(Y,params,0);
        Yc=clamp_data(Y,params,1);
        Y=Yr;
        Xc=clamp_data(X,params,0);
        X=Xc;
    end
    %rng(params.seed);
    
     % permute data points
    if(0)
        [n,d]=size(X);
        z=rand(1,n);
        [s,i]=sort(z,2,'descend');
        X=X(i,:);
        Y=Y(i,:);
    end
    %}
    
    %drug_stats=get_drug_stats(Y,params); % NB: the drug statistics were computed using the Python implementation
    
    stats.x_max=nanmax(X,[],1);
    stats.y_max=nanmax(Y,[],1);
    stats.x_min=nanmin(X,[],1);
    stats.y_min=nanmin(Y,[],1);
    
  %  Y=ynormalize(Y);
   % Y=rescaling(Y,0); 
    data.X = X;%(index,:);
    data.Y = Y;%(index, :);
    
    
    
   %%%%%%%%%%%%%%%%
   if(params.anal_data)
    analyze_data(data,params);
   end
   %%%%%%%%%%%%%%%%
   
   
   
    %%%%privacy setting
    params.privacy_iter = 1;
   % 
    
    all_data = struct();
    list_eps = [0.0001 0.2 1 10 1000 100000];
    list_eps = [1];
    all_data.list_eps = list_eps;
    eFold = length(list_eps);
    
    kFold = params.fold;
    
    all_data.tr_data = cell(kFold);
    all_data.te_data = cell(kFold);
    all_data.pri_data = cell(kFold);
    all_data.noise_data = cell(kFold, eFold, params.privacy_iter);
     
    n = size(data.X, 1);
    selected = zeros(1, n);
    
    params.old_normalizey=0;
    
    for kIndex = 1:kFold
        rng(info.seed(kIndex));
        [test_index, selected] = split_train_and_test_random(n, selected, params);
        
        %split public and private data
        te.X = data.X(test_index, :);
        te.Y = data.Y(test_index, :);
        
        te.Y=ynormalize(te.Y,params);
        
        [rY, cY] = preprocessY(te.Y, params);
        
        stats.y_min=nanmin([stats.y_min;cY],[],1); %cY is considered for rmse
        stats.y_max=nanmax([stats.y_max;cY],[],1);
        
        te.cY = cY;
        te.rY = rY;
        
        tr.X = data.X(~test_index, :);
        tr.Y = data.Y(~test_index, :);
        tr.Y=ynormalize(tr.Y,params);
        
       [rY, cY] = preprocessY(tr.Y, params);
        
        tr.rY = rY;
        tr.cY = cY;
        stats.y_min=nanmin([stats.y_min;cY],[],1);
        stats.y_max=nanmax([stats.y_max;cY],[],1);
        
        total_tr_size = size(tr.X, 1);
       % tr_size = total_tr_size - params.prsize;
        tr_size = params.insize;
       
        tr_index = datasample(1:tr_size, tr_size, 'Replace', false);
        pri_index = true(1, total_tr_size);
        pri_index(tr_index) = false;
        pri_data.X = tr.X(pri_index, :);
        pri_data.cY = tr.cY(pri_index, :);
        pri_data.rY = tr.rY(pri_index, :);
        tr.X = tr.X(~pri_index,:);
        tr.cY = tr.cY(~pri_index,:);
        tr.rY = tr.rY(~pri_index,:);
        
        
        all_data.tr_data{kIndex} = tr;
        all_data.te_data{kIndex} = te;
        all_data.pri_data{kIndex} = pri_data;
        %pridataX = pri_data.X * 0;
        
      %  keyboard;
    end
    
    %%private index
    pr_size = 0:100:params.prsize;
    size_len = length(pr_size);
    pri_index = cell(size_len);
    pri_index{1} = ones(1, 0);
    selected = zeros(1, params.prsize);
    
    for sIndex = 2:size_len
        available_index = 1:params.prsize;
        available_index = available_index(~selected);
        sample_index = datasample(available_index, 100, 'Replace', false);
        pri_index{sIndex} = sample_index;
        selected(sample_index) = 1;
    end
    all_data.pri_index = pri_index;
    %{
    all_data.x_max=x_max;
    all_data.x_min=x_min;
    all_data.y_max=y_max;
    all_data.y_min=y_min;
    %}
    save([info.odir params.processed_data_filename], 'all_data');
    save([info.odir,'stats.mat'], 'stats');
    %save([params.odir,'drug_stats.mat'],drug_stats); % NB: the drug statistics were computed using the Python implementation
    %%% 
    get_testingY(params,info);
    %%% 
end

function Z=ynormalize(Y,params)
    [n, m]=size(Y);
    Z=Y;
    if(params.old_normalizey==0)
        Z=Y-ones(n,1)*nanmean(Y,1);
    end
    
    %Z=Z./sqrt(ones(n,1)*nansum(Z.^2,1));
   % Z=Z/n;
end
function [rY, cY] = preprocessY(Y, params)
    rY=Y; cY=Y;
    
    if(params.normalize && params.old_normalizey)
        if(params.l2normalize || 1)
            rY = ml2normalize(Y);
            cY = ml2normalize(Y,1);
        end

        if(params.l1normalize)
            rY = L1normalize(Y,params.l1norm_dim);
            cY = L1normalize(Y,params.l1norm_dim);
        end

        if params.real2rank
            rY = real2Rank(rY);
        end


        if params.real2rank
            cY = real2Rank_bycolumn(cY);
        end
    end
end
            

    
    


