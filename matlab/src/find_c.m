function [c]=find_c(d,n,oX,oY)

    synth=1;
    
    n=500;
    eps = 1;
    theta = 1e-1; %1/0.1^2;
    lambda = 1e0; %1/0.3^2;
    
    L1=10;
    L2=10;
   %{ 
    cs1=0.01:0.01:0.1;
    cs2=1:1:1e1;
    cs2=[0.1 0.5 cs2];
    %}
    st=0.1:0.1:1;
    lenC1=length(st);
    lenC2=lenC1;
    
    
    params=initial_params;
    
    fc1=0; fc2=0; 
    fci1=0; fci2=0;
    sch=0;	

    if(sch==0)		
    	err=zeros(lenC1,lenC2);
    end
    for l1=1:L1 %repeating over data
    if(sch==1)		
    	err=zeros(lenC1,lenC2);
    end
       
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
                  
                  U = wishrnd((2*c1*c1*eye(d))/eps,d+1);
                  v=sqrt(2)*d/eps;
                  v=v*(4*c1*c2);
                  V=laprnd(d,1,0,v);

                  Xc=clamp_data(X,c1); 
                  yc=clamp_data(y,c2);

                  xx=(Xc'*Xc+U); xy=Xc'*yc+V;
                  w2=(xx+lambda*eye(d))\xy;
                  ind=(ci1-1)*lenC2+ci2;
                  py=X*w2;
                  rho = corr(py,y,'type','Spearman','rows','complete');
                  %rmse = mean((py-y).^2);  
                  err(ci1,ci2)=err(ci1,ci2) + rho;
                end
            end
        end
        if(sch==1) 
    		err=err/(L2);
		%[ci1,ci2]=find(err==min(err(:))); 
		[ci1,ci2]=find(err==max(err(:))); 
        	fc1=fc1+cs1(ci1);
        	fc2=fc2+cs2(ci2);
        	fci1=fci1+ci1;
        	fci2=fci2+ci2;
	end
    end
    if(sch==0)
    	err=err/(L1*L2);
	%[ci1,ci2]=find(err==min(err(:))); 
	[ci1,ci2]=find(err==max(err(:))) 
    end 	
    if(sch==1)	
    	fc1=fc1/L1
    	fc2=fc2/L1
    	fci1=fci1/L1
    	fci2=fci2/L1
    end	
    keyboard;
end


function Z=ynormalize(Y)
    [n, m]=size(Y);
     Z=Y-ones(n,1)*nanmean(Y,1);
   % Z=Z./sqrt(ones(n,1)*nansum(Z.^2,1));
   % Z=Z/n;
end
    
