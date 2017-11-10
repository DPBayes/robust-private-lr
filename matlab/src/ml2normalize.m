function [Z]=ml2normalize(X, byCol)
    [n, d]=size(X);
   % Z=X-nanmean(nanmean(X));
    Z = X - ones(n,1)*nanmedian(X);
    if(1)
    if exist('byCol', 'var')
      %  disp('normalize by col');
        Z=Z./sqrt(ones(n,1)*nansum(Z.^2));
    else
        Z=Z./sqrt(nansum(Z.^2,2)*ones(1,d));
    end
    end
end
