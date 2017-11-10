function [synth_data]=gen_synth_data(params)
    n=params.total_size;%2*1e2;
    d=1e1;
    m=116;
    K=1;
    %n=n+params.prsize;
         
    
    config=1;
    params.synth_dist=1;
    
        if(config==1)
            X=rand(n,d); %*1e2+1e1;
            W=randn(d,m);%*1e1;
            E=randn(n,m);%*1e3;
            Y=X*W + E;
            
            if(0)
                Z=(rand(n,d)<0.2);
                X(Z)=0; 
                X=X+rand(n,d);

                Z=(rand(n,m)<0.2);
                Y(Z)=0;
                Y=Y+rand(n,m);

                Z=(rand(n,m)<0.2);
                Y(Z)=NaN;

                Y=clamp_data(Y,1);

                Z=(rand(n,d)<0.2);
                X(Z)=X(Z)*1e1;
            end            
            %X=clamp_data(X,params,0);
            %Y=randn(n,m); %*1e2;
        end
    
        if(config==2)
            mu=rand(1,K)*1e1;
            b=(rand(1,K) < 0.5);
            mu=mu+b*1e1;
            W=rand(d,m*K);
            A=ones(1,d)*5;
            B=ones(1,d)*2;
            for i=1:n
                z=round(rand*(K-1))+1;
                if(params.synth_dist==1)
                    x=randn(1,d)*1e1+mu(1,z);
                end
                if(params.synth_dist==2)
                    x=betarnd(A,B)*1e1;
                end
               if(params.synth_dist==3)
                    x=rand(1,d)*1e2;
                end 
                w=W(:,(z-1)*m+1:z*m);
                Y(i,:)=x*w;
                %y=(rand(1,d)<0.0);
                %x(y)=NaN;
                X(i,:)=x;
            end
            
        
            %E=randn(n,m)*1e1;
            %Y=X*W + E;
        end
   
    
        
        
        
        synth_data.X=X;
        synth_data.Y=Y;
        synth_data.W=W;
