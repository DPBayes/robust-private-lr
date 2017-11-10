function save_tensor(params,info)

    if(info.sch==0) %in
        sdir=[info.odir,'/sch0'];
        d=info.in_ind;
    end
    if(info.sch==1) %d
        sdir=[info.odir,'/sch1'];
        d=info.ds_ind;
    end
    if(info.sch==2) %eps
        sdir=[info.odir,'/sch2'];
        d=info.eps_ind;
    end
    if(info.sch==3) %c
        sdir=[info.odir,'/sch3'];
        d=info.c_ind;
    end
    
    load([info.odir 'result.mat']);
    load([sdir,'/tensor.mat']);

    for j = 1:2
        [m, s] = concate_result(result, j);
            x=squeeze(T(j,d,:));x=x';
            x=x+[m zeros(1,info.width-length(m))];
            T(j,d,:)=x;
    end
    save([sdir,'/tensor.mat'], 'T');
     
  
end

function [m, s] = concate_result(result, index)
    stand = NaN;

    if index == 1
        stand = result.lr;
    elseif index == 2
        stand = result.rplr;
    end
    
     sLen = size(stand, 1);
    kFold = size(stand, 2);
   
    
    res = NaN(kFold, sLen);
    for i = 1:sLen
        res(:, i) = [stand{i,:}];
    end
    m = nanmean(res);
   s = std(res);
end
