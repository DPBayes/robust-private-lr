function [features] = feature_selection(X, Y, params)
    disp('feature selection starting...');
    [n, d] = size(X);
    
    features = [];
    if(params.dimension_reduce == 1)
        %after rescaling, variance based selection
        stdX=std(X, 1, 1);%for sanger data, max:0.3048, min:0.047
        R = [];
        C = [];
        k = 1;
        oldS = 0;
        std_arr = min(stdX):0.01:max(stdX);
        for st = std_arr
            features = [];
            s = 1;
            for i = 1:d
                if(stdX(i) >= st) 
                	features(s) = i;
                	s = s + 1;
                end
            end
            if oldS == 0 || oldS ~= s
                oldS = s;
                %train, predict and test, find the feature set reach the maximum
                %accuracy
		data.X = X(:, features);
		data.Y = Y;
		if(params.l1norm)
			data.X = l1normalize(data.X);
			data.Y = l1normalize(data.Y);

		end

		if(params.l2norm)
			data.X = l2normalize(data.X);
			data.Y = ml2normalize(data.Y);
		end

		cr = [];
		cc = [];	
		for iter = 1:params.cross_valid
			[tr, te] = split_train_and_test(data, iter, params);
			
			[~, result] = train_prediction_test(tr, te, params);
			cr(iter) = result.C.R;
			cc(iter) = result.C.avg_rho_p;
			if params.corr_drug
				cc(iter) = result.C.avg_rho; 
			end
		end
		r = mean(cr);
		c = mean(cc);
                R(k) = r;
                C(k) = c;
                save(sprintf('%s/features_%f.mat', params.result_dir, st), 'features', 'r', 'c');
                k = k + 1;
            end
        end
        
        [c, index] = max(C);
        load(sprintf('%s/features_%f.mat', params.result_dir, std_arr(index)));
        display([c, std_arr(index), length(features)]);
        save(sprintf('%s/features_max.mat', params.result_dir), 'features', 'c');
    elseif(params.dimension_reduce == 2)
        %%%load the gene set that on the top of variance rank list and the
        %%%set lead to optimal accuracy based on the wrapper method for feature
        %%%selection.
       load('features_sanger.mat');
    elseif params.dimension_reduce == 4
       %%%load the pre-selected gene set based on lasso
	   load(params.GE_feature_file);
    elseif(params.dimension_reduce == 3)
        %%%using lasso to find the optimal feature set for each task, and
        %%%choose the largest feature set as the one for all tasks, since
        %%%there is a paper showed that randonly choose a equal size of
        %%%gene set could get a result with no significant different result
        %%%from that based on the comprehensively selected gene set.
        m = size(Y, 2);		
        X = l2normalize(X);
            Y = ml2normalize(Y, 1);
        mDF = [];
        parfor i = 1:m
            [B stats] = lasso(X, Y(:,i), 'CV', 5);
            mDF(i) = stats.DF(stats.IndexMinMSE);
            save_file(params.result_dir, i, B, stats);
        end
        [v, i]=max(mDF);
        res = load(sprintf('%s/lassofs_result_%d.mat', params.result_dir, i));
        features = find( res.b(:,res.stats.IndexMinMSE) ~= 0);
        save(sprintf('%s/lasso_features_max.mat', params.result_dir), 'features');
    end
    disp(['feature selection done. Feature size:', num2str(length(features))]);
end
