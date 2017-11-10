function [Z]=L1normalize(X,dim,params)
    [n, d]=size(X);
    Z = X;
    if(params.normalize)
        if(dim==1)
            X = X - ones(n,1)*nanmean(X,1);
           Z=X./(ones(n,1)*nansum(X,1));
        else
            X = X - ones(n,1)*nanmean(X,1);
            Z=X./(nansum(X,2)*ones(1,d));
        end
    end
end