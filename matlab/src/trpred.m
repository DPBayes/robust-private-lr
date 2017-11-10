function pred = trpred(intr, extr, te, params, info)
    if(params.modell==2)
        if(size(extr.X,1)>0)
             load([info.odir,'stats.mat']);
             datamin=[stats.x_min stats.y_min];
             datamax=[stats.x_max stats.y_max];
   
             S=extr.X;
             T=extr.Y;
   
             if(params.clip==1)
                 max_val=nanstd(S(:))*0.3;
                 max_val2=nanstd(T(:))*0.4;
   
                 S=clamp_data(S,max_val);
                 T=clamp_data(T,max_val2);
                 datamin=[-max_val*ones(size(stats.x_min)) -max_val2*ones(size(stats.y_min))];
                 datamax=[max_val*ones(size(stats.x_max)) max_val2*ones(size(stats.y_max))];
             end
             
             fm_data = [S T]; %extr.X extr.Y];
             fm_data=preprocess_fm(fm_data,datamin,datamax);
             bnd=size(S,2);
             extr.X=fm_data(:,1:bnd);
             extr.Y=fm_data(:,bnd+1:end);
             clear fm_data;

             datamin=stats.x_min;
             datamax=stats.x_max;
             fm_data = te.X;
             fm_data=preprocess_fm(fm_data,datamin,datamax);
             te.X=fm_data;
             model= train(intr, extr, params, info);
             pred = predict(model, te, params);
        else
           model= train(intr, extr, params, info);
           pred = predict(model, te, params);
        end
    else
        model= train(intr, extr, params, info);
        pred = predict(model, te, params);
    end
    
  
end

function model = train(intr, extr, params, info)

    X=intr.X; % n x d
    Y=intr.Y; % n x m
    
    S=extr.X;
    T=extr.Y;
    
    laplace=1; % use Laplace mechanism
    wishart=0;
   
    eps=info.eps; %*2; % multiplication by 2 to make epsilon comparable with unbounded DP
   
    %{
    % c=0.1;
    max_val=info.c;
    %max_val=0.07; %found by find_c
    max_val=0.05; %found by find_c
    
    max_val2=info.c;
    %max_val2=1; %found by find_c
    max_val2=5; %found by find_c
    max_val2=1.5; %based on s.d.
    %}
       
    max_val=nanstd(S(:))*0.3;
    max_val2=nanstd(T(:))*0.4; %0.4;
    max_val3=max_val;
    
    if(max_val==0 || params.clip==0)
        %need to set to compute sensitivity
    %    bound01=0;
        max_val=1;
        max_val2=1;
        max_val2=20; % new normalization of y
        max_val3=max_val;
    end
    
   if(params.clip==1)
       %apply clipping
        S=clamp_data(S,max_val);
        T=clamp_data(T,max_val2);
   end
   
   if(params.intercept) %better without intercept
        [in, d] = size(X);
        [en, d] = size(S);

        X=[ones(in, 1) X];
        S=[ones(en, 1) S];
    end
    
    [in, d] = size(X);
    [en, d] = size(S);
    
    m=size(Y,2);
  
    lambda = params.lambda;
    sigma = params.sigma;
    theta = params.theta;
    
    
    W=zeros(d,m);
    
    for i=1:m
        % consider i-th drug
        y=Y(:,i);
        
        % ignore missing values
        [iX,iy]=ignore_nan(X,y);
        
        % process ext data
        if(en>0)
           
                
            t=T(:,i); 
            [iS,it]=ignore_nan(S,t);
           %{
           if(params.clip==2)
                sm=nanmax(iS,[],2); cm=nanmax(sm,it);
                cm(cm<max_val)=max_val;
                cma=(cm*ones(1,d));
                iS=(iS./cma)*max_val;
                it=(it./cm)*max_val;
           end
           %}     
                if(params.rplr)
                    
   
                    %<<< Laplace
                    if(laplace)
                        nq=d*(d+1); % !
                        var=sqrt(2)*nq/eps;
                        var=var*(max_val*max_val)*2; % !
                        u=laprnd(d,d,0,var); 
                        U = diag(diag(u)) + tril(u,-1) + tril(u,-1)';
                    end
                    %>>> Laplace

                    %<<< Wishart 
                    if(wishart)
                        U = wishrnd((2*d*max_val*max_val*eye(d))/eps,d+1); %epsilon/2 bounded DP, epsilon/4 unbounded DP
                    end
                    
                    SS = iS'*iS;
                    Sigma=sigma*(SS+U);


                    %<<< perturb 2
                    nq=d;
                    var_xy=sqrt(2)*nq/(eps);
                    var_xy=var_xy*(4*max_val2*max_val3); %epsilon/2 bounded DP, epsilon/4 unbounded DP
              
                    if(laplace || wishart)
                        V=laprnd(d,1,0,var_xy);
                    end

                     ST = iS'*it; %clamp_data(iS,max_val3)'*it;
                     Mu=sigma*(ST+V);
                else

                    %<<< without privacy
                   U=zeros(d,d);
                   SS = (iS'*iS);
                   
                   Sigma=sigma*(SS+U);

                   V=zeros(d,1);
                   ST =  iS'*it;
                   Mu=sigma*(ST+V);
                end

         else
            %<<< with only internal data
            Sigma = zeros(d,d); %lambda * eye(d);
            Mu = zeros(d,1);
         end
        

        if(params.no_internal)
            theta=0;
        end
            
            
        if(params.modell==1)
            Psi = theta*(iX'*iX) + (Sigma + lambda*eye(d));
            Nu = Psi\(Mu + theta * iX'*iy);
            w=Nu;
        end
        
        if(params.modell==2)
           
            if(en>0)
                 fm_data=[iS it];
                [w, b] = Functional_Linear(fm_data, eps);
            else
                Psi = theta*(iX'*iX) + (Sigma + lambda*eye(d));
                Nu = Psi\(Mu + theta * iX'*iy);
                w=Nu;
            end
        end
        
        if(params.modell==3)
            if(en>0)
                addpath(genpath(pwd));
                w=LR_output_perturbation(eps,iS,it);
            else
                Psi = theta*(iX'*iX) + (Sigma + lambda*eye(d));
                Nu = Psi\(Mu + theta * iX'*iy);
                w=Nu;
            end
        end

        if(params.modell==4)
            if(en>0)
                b=lasso([iX;iS],[iy;it],'Lambda',0.1);
            else
                b=lasso(iX,iy,'Lambda',0.1);
            end
            w=b(:,1);
        end
        W(:,i)=w;
    end
    
    model.W=W;
end


function pred = predict(model, te, params)
    W = model.W;
    if(params.intercept)
        [tn, d] = size(te.X);
        X=[ones(tn,1) te.X];
    else
        X=te.X;
    end
    pred = X*W;
end
