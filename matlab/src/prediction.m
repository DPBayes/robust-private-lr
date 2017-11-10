function [pred]=prediction(model,tr,te,params)
    trY = tr.Y;
    trX = tr.X;
    teX = te.X;
    
    nte = size(teX, 1);
    P = [];
    if params.model_mean
        disp('using model mean for prediction ...');
        p = nanmean(trY);
        P = ones(nte, 1)*p;
    elseif params.model_nn
        disp('using model NN for prediction ...');
        for i=1:nte
            l=findnn(trX,teX(i,:));
            %originally from mrinal, here is P(i,:)=Ytr(i,:);
            P(i,:)=trY(l,:);
        end
    else
        W = model.W;
        P = teX*W;
        if params.model_logistic
            sigmoid=@(x)(1./(1+exp(-x)));
            P = sigmoid(P);
        end
    end
    
    if(params.chk_model)
        pred.Y = model.W;
    else
        pred.Y = P;
    end
end
