function [test_index, selected]=split_train_and_test_random(num, selected, params)
    tsize = params.tesize; 
    rng(floor(rand*5000));
    available_index = 1:num;
    availabel_index = available_index;
    %availabel_index = available_index(~selected);
    
    sample_index = datasample(availabel_index, tsize, 'Replace', false);
    test_index = false(num, 1);
    test_index(sample_index) = true;
    selected(sample_index) = 1;
end

    
    
    