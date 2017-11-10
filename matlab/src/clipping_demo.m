function [err_rmse, err_rho]=clipping_demo(d,n)

    synth=1;
    
    if nargin < 2,
      n=500;
    end
    eps = 1;
    theta = 1e0; %1/0.1^2;
    lambda = 1e0; %1/0.3^2;
    
    L1=10;
    L2=10;
    
    st=[0.05, 0.1:0.1:2];
    lenC1=length(st);
    lenC2=lenC1;
    
    err_rho=zeros(lenC1, lenC2);
    err_rmse=zeros(lenC1, lenC2);
    params=initial_params;
    
    
    fc1=0; fc2=0; 
    fci1=0; fci2=0;
    for l1=1:L1 %repeating over data
        fprintf('doing repeat %d/%d\n', l1, L1);
        if(synth)
            X = randn(n,d);
        else
            X=oX;
        end
        
        X = L2normalize(X,2,params);
        sx=std(X(:));
        
        if(synth)
            w = randn(d,1)./sqrt(lambda); 
            y = X*w + randn(n,1)./sqrt(theta);
        else
            y=oY(:,1);
              [X,y]=ignore_nan(X,y);
        end
        
        y=ynormalize(y);
        sy=std(y);

        cs1=st*sx;
        cs2=st*sy;

        for l2=1:L2 %repeating over DP noise
            
            for ci1=1:lenC1
                c1=cs1(ci1);  
        
                for ci2=1:lenC2
                  c2=cs2(ci2);  
                  
                  U = wishrnd((2*d*c1*c1*eye(d))/eps,d+1);
                  v=sqrt(2)*d/eps;
                  v=v*(4*c1*c2);
                  V=laprnd(d,1,0,v);

                  Xc=clamp_data(X,c1); 
                  yc=clamp_data(y,c2);

                  xx=(Xc'*Xc+U); xy=Xc'*yc+V;
                  w2=(xx+lambda*eye(d))\xy;
                  py=X*w2;
                  rho = corr(py,y,'type','Spearman','rows','complete');
                  rmse = mean((py-y).^2);
                  err_rho(ci1, ci2)=err_rho(ci1, ci2) + rho;
                  err_rmse(ci1, ci2)=err_rmse(ci1, ci2) + rmse;
                end
            end
        end
    end
    err_rho=err_rho/(L1*L2);
    err_rmse=err_rmse/(L1*L2);
end


function Z=ynormalize(Y)
    [n, m]=size(Y);
     Z=Y-ones(n,1)*nanmean(Y,1);
   % Z=Z./sqrt(ones(n,1)*nansum(Z.^2,1));
   % Z=Z/n;
end
