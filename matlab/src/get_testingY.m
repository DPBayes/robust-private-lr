function get_testingY(params,info)
    %load all_data
    %load('preprocessed_all_data_random_bycolumn.mat');
    params = initial_params();
    load([info.odir params.processed_data_filename]);
    [kFold, ~] = size(all_data.te_data);
    
    if(params.chk_model)
       load([info.odir params.data_filename]);
        W=synth_data.W;
    end
        
    
    test=struct();
    for kIndex = 1:kFold
        
        if(params.chk_model)
            lY=W; clY=W;
        else
            te = all_data.te_data{kIndex};
            [n, m] = size(te.rY);
            lY = NaN(n, m);
            clY = NaN(n, m);

            if(params.rmse)
                lY=te.rY;
                clY=te.rY;
            else
                for i = 1:n
                    if sum(~isnan(te.rY(i,:))) >= params.topK
                        lY(i,:) = label_topk(te.rY(i,:), params.topK);
                    end
                end
                for i = 1:m
                    if sum(~isnan(te.cY(:,i))) >= params.topK
                        clY(:, i) = label_topk(te.cY(:,i), params.topK);
                    end
                end
            end  
        end
        %test.Y{kIndex}=te.Y;
        test.lY{kIndex} = lY;
        test.clY{kIndex} = clY;

    end

    %save('../result_consis_median/testY.mat', 'test');
    save([info.odir 'testY.mat'], 'test');
end

function get_testingY_fixed(params,info)
    if ~exist('params', 'var')
        disp('load params');
    	params = initial_params();
    end
    all_data = data_prepare(params);

    %%%%cross validation
    params.k_fold = 1;
    if(params.cross_valid)
        params.k_fold = 5;
    end
    test=struct();
    for kIndex = 1:params.k_fold
	
        %train and test data preparation
        disp('train and test data splitting');
        [~, te] = split_train_and_test(all_data, kIndex, params);
        [n, m] = size(te.Y);
        lY = NaN(n, m);
        clY = NaN(n, m);
        for i = 1:n
            if sum(~isnan(te.Y(i,:))) > 10
                lY(i,:) = label_topk(te.Y(i,:), 10);
            end
        end
        for i = 1:m
            if sum(~isnan(te.Y(:,i))) > 10
                clY(:, i) = label_topk(te.Y(:,i), 10);
            end
        end
        test.Y{kIndex}=te.Y;
        test.lY{kIndex} = lY;
        test.clY{kIndex} = clY;
       
    end

    [~, data_filename, ~] = fileparts(params.data_filename);
    save([info.odir data_filename '_testY.mat'], 'test');
end

function [z] = label_topk(x, k)
    [~, rank_index] = sort(x);
    nanIndex = isnan(x);
    z = repmat(-1, [1, length(x)]);
    z(rank_index(1:k)) = 1;
    z(nanIndex) = NaN;
end