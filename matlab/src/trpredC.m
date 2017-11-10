
function pred = trpredC(intr, extr, te, params, info)
    if(params.modell==2)
        if(size(extr.X,1)>0)
             load([params.odir,'stats.mat']);
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
             pred = predict(model, te, params, info);
        else
           model= train(intr, extr, params, info);
           pred = predict(model, te, params, info);
        end
    else
        model= train(intr, extr, params, info);
        pred = predict(model, te, params, info);
    end
    
    %size(extr.X)
  

end

function model = train(intr, extr, params, info)

    X=intr.X; % n x d
    Y=intr.Y; % n x m
    
    S=extr.X;
    T=extr.Y;
    %{
    varT=nanvar(T,[],1);
    [val,ind]=max(varT);
    ind=1;
    best_t=T(:,ind); %not used now
    %}
    
    laplace=0;
    wishart=1;
   
    eps=info.eps; %*2; % multiplication by 2 to make epsilon comparable with unbounded DP
   % c=0.1;
    max_val=info.c;
    %max_val=0.07; %found by find_c
    max_val=0.05; %found by find_c
    
    max_val2=info.c;
    %max_val2=1; %found by find_c
    max_val2=5; %found by find_c
    max_val2=1.5; %based on s.d.
    
       
    max_val=nanstd(S(:))*0.3;
    max_val2=nanstd(T(:))*0.3; %0.4;
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
       % X=clamp_data(X,max_val);
        S=clamp_data(S,max_val);
       % Y=clamp_data(Y,max_val);
        T=clamp_data(T,max_val2);
   else
       %{
       X = L1normalize(X);
       S = L1normalize(S);
       Y = L1normalize(Y);
       T = L1normalize(T);
       %}
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
  
    if(params.extrapolate)
        scale=params.scale*(params.extrapolate+1);
        %scale=50;
     else
        scale=1;
     end

    lambda = params.lambda;
    sigma = params.sigma;
    theta = params.theta;
    
    hs = params.hs; hs_frac=1; %params.hs;
    
    hs_frac_in=1;
    
    %if(hs==1) hs_frac=0.5; end %[0.01 for synth. ]
    if(hs==1 || hs==0) hs_frac=1; end % 0.05, 0.1 work for real, [0.01 for synth. ]
   
   % if(params.debug) m=1; end
    
    if(hs==2)
      hs_frac = 0.5;  
      %hs_frac = dpseqinf_opt_noise(1, in, en, d, 10);
    end
    
   % disp(['en ',int2str(en),', hs_frac ',num2str(hs_frac),', scale ',int2str(scale),', extrp ',int2str(params.extrapolate)]);       
    W=zeros(d,m);
    
    for i=1:m
       % disp(['drug ',int2str(i)]);

        % consider i-th drug
        y=Y(:,i);
        
        % ignore missing values
        [iX,iy]=ignore_nan(X,y);
        
        % process ext data
        if(en>0)
           
                
            t=T(:,i); 
            [iS,it]=ignore_nan(S,t);
           
           if(params.clip==2)
                sm=nanmax(iS,[],2); cm=nanmax(sm,it);
                cm(cm<max_val)=max_val;
                cma=(cm*ones(1,d));
                iS=(iS./cma)*max_val;
                it=(it./cm)*max_val;
           end
                
                %%
                if(params.dpvb)
                    
                   
                    %{
                    %% Gaussian
                    if(0)
                        delta=1e-1;
                        Delta=1; 
                        nq=d*(d+1)/2;
                        var=sqrt(nq*log(1/delta)); 
                        u=randn(d,d)*var; 
                    end
                    %%
                    %}

                    %% Laplace
                    if(laplace)
                        nq=d*(d+1)/2;

                        var=sqrt(2)*nq/eps;



                        var=var*(max_val*max_val);
                        u=laprnd(d,d,0,var); 
                        %U=(u+u')/2;
                        U = diag(diag(u)) + tril(u,-1) + tril(u,-1)';

                    end
                    %%

                    %% Wishart 
                    if(wishart)
                        U = wishrnd((2*d*max_val*max_val*eye(d))/eps,d+1); %epsilon/2 bounded DP, epsilon/4 unbounded DP
                    end

                    %{
                    %% diagonal 
                    if(0) 
                        nq=d; var=sqrt(2)*nq;
                        u=laprnd(d,1,0,var);
                        U=diag(u);
                    end
                    %%
                    %}


                    %U=zeros(d,d);

                    SS = scale*(iS'*iS);

                    %{
                    %only diagonal
                    if(0) %% not-working
                        diag_SS=diag(SS);
                        SS=diag(diag_SS);
                    end
                    %}
                    
                    Sigma=hs_frac*sigma*(SS+U);


                    %%
                    nq=d;
                    %{
                    if(0)
                        var=sqrt(d*log(1/delta));
                    end
                    %}
                    
                    var_xy=sqrt(2)*nq/(eps);

                    var_xy=var_xy*(4*max_val2*max_val3); %epsilon/2 bounded DP, epsilon/4 unbounded DP
                    %v=randn(d,1)*var;

                    if(laplace || wishart)
                        V=laprnd(d,1,0,var_xy);
                    end

                    %V=zeros(d,1);
                    %%
                    %%
                     ST = scale * iS'*it; %clamp_data(iS,max_val3)'*it;
                     Mu=hs_frac*sigma*(ST+V);
                     
                    % Mu=zeros(d,1);
                    %%
                else

                    %% without privacy
                   %hs_frac=0.1;%scale=1;
                   %hs_frac=1; %[1 better] for similar model it should not be same as DP case not 1

                   U=zeros(d,d);
                   SS = scale*(iS'*iS);
                   %only diagonal
                    if(0) %% not-working
                        diag_SS=diag(SS);
                        SS=diag(diag_SS);
                    end

                   Sigma=hs_frac_in*sigma*(SS+U);

                   V=zeros(d,1);
                   ST = scale * iS'*it;
                   Mu=hs_frac_in*sigma*(ST+V);

                   %%
                end

         else
            %% with only internal data
            Sigma = zeros(d,d); %lambda * eye(d);
            Mu = zeros(d,1);
            % W=Mu;
            %%
        end
        
        
            scale_lambda=1; %[10 for synth. ]
            if(en>0) scale_lambda=1; end

            %scale_lambda=10; %[10 is working]
            %scale_lambda=1e2; %stronger prior
            %scale_lambda=0.01; %weaker prior
            
            if(en>0)
                post_scale=1; %1/(scale*en); %[1 better] 1/(scale*en);
            else
                post_scale=1;
            end
            
            if(params.no_internal)
                theta=0;
            end
            
            %%{
            if(en>0 && params.dpvb && params.denoise)
                val=dpdenoise(d, 1, in, en, (iX'*iy)', Mu', var_xy*var_xy*eps);
                Mu=val; %Mu*scale_xy;
            end
            %%}
            
        if(params.modell==1)
            Psi = theta*(iX'*iX) + (post_scale*Sigma + scale_lambda*lambda*eye(d));
            Nu = Psi\(Mu + theta * iX'*iy);
            w=Nu;
        end
        
        if(params.modell==2)
           
            if(en>0)
                 fm_data=[iS it];
                [w, b] = Functional_Linear(fm_data, eps);
            else
                Psi = theta*(iX'*iX) + (post_scale*Sigma + scale_lambda*lambda*eye(d));
                Nu = Psi\(Mu + theta * iX'*iy);
                w=Nu;
            end
        end
        
        if(params.modell==3)
            if(en>0)
                addpath(genpath(pwd));
                w=LR_output_perturbation(eps,iS,it);
            else
                Psi = theta*(iX'*iX) + (post_scale*Sigma + scale_lambda*lambda*eye(d));
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
           % w
        end
        W(:,i)=w;
    end
    
    model.W=W;
end


function pred = predict(model, te, params, info)
    W = model.W;
    if(params.intercept)
        [tn, d] = size(te.X);
        X=[ones(tn,1) te.X];
    else
        X=te.X;
    end
    pred = X*W;
end
