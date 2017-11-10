function [Z]=L2normalize(X,dim,params)
    [n, d]=size(X);
    Z = X;
    
    if(params.normalize)
        if(dim==1)
             Z = X;
            if n > 1
                X = X - ones(n,1)*nanmean(X,1);
                Z=X./(ones(n,1)*sqrt(nansum(X.^2,1)));
            end

        else
            X = X - ones(n,1)*nanmean(X,1);

            Z=X./(sqrt(nansum(X.^2,2))*ones(1,d));
        end
    end
end