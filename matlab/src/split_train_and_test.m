function [train, test]=split_train_and_test(data, k, params)
    num = size(data.X, 1);
    tsize = round(num * (1 - params.frac));

    endIndex = tsize * k;
    if(endIndex > num)
        endIndex = num;
    end
    test_index = ((k - 1)*tsize+1):endIndex;
    train_index = [1:((k - 1)*tsize), (endIndex + 1):num];  
    test.X = data.X(test_index, :);
    test.Y = data.Y(test_index, :);
    
    train.X = data.X(train_index, :);
    train.Y = data.Y(train_index, :);
end

    
    
    