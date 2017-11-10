function [X1]=clamp_data(X,c)

% X(X>c)=c;
 X1=X; X1(X>c)=c; 
 X1(X<-c)=-c;
 
 
  %{
    if(byColumn)
            X=X';
    end    
    [n,d]=size(X);
    params.lb=0;
    
    %X=X-nanmean(X,2)*ones(1,d);
    
    if(params.lb==-1) % range [-1,1]
        l=-1*nanmin(X,[],2)*ones(1,d);
        u=nanmax(X,[],2)*ones(1,d);
        ind1=X<0; ind2=X>=0;
        l(ind2)=1;
        u(ind1)=1;
        X=X./l;
        X=X./u;
    else % range [0,1]
        if(0 && params.clamp_exp)
            X=ones(n,d)./(1+exp(-X));
        else
            l=-1*nanmin(X,[],2)*ones(1,d);
            X=X+l;
            u=nanmax(X,[],2)*ones(1,d);
            X=X./u;
        end
    end
    if(byColumn)
            X=X';
    end 
    %}