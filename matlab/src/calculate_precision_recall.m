function calculate_precision_recall(byColumn, params, info)
    result = struct();
    predict = struct();
    
    prefix = '';
    
    if params.real2rank == 1
        prefix = 'rank_';
    end
    
    if params.rmse
        prefix='rmse_';
    end
    
    orr = [prefix 'OPT_RR.mat'];
    [result.orr, predict.orr]= single_file(orr, 0, byColumn, params, 1);
    
    %disp('rr');
    rr = [prefix 'noDP_RR.mat'];
    [result.rr, predict.rr]= single_file(rr, 0, byColumn, params, 2);
    %disp('==========================================================');
    %disp('irr');
    irr = [prefix 'iDP_RR.mat'];
    [result.irr, predict.irr] = single_file(irr, 0, byColumn, params, 4);
    %disp('==========================================================');
    %disp('dpvb');
   % ops1 = [prefix 'ops_sch1.mat'];
   dpvb = [prefix 'dpvb.mat'];
    if(params.no_models>2)
        [result.dpvb,predict.dpvb]=single_file(dpvb, 0, byColumn, params, 3);
    end
    %ops2 = [prefix 'ops_sch2.mat'];
    %[result.ops2,predict.ops2]=single_file(ops2, 0, byColumn, params);
    if(0)
        fm = [prefix 'iDP_FM.mat'];
        [result.fm, predict.fm]= single_file(fm, 0, byColumn, params, 1);
    end
    if byColumn == 1
            %save([params.odir 'result_rr_ops_ranking_patients.mat'], 'result', 'predict');
            save([params.odir prefix 'result.mat'], 'result', 'predict');
    else
             %save([params.odir 'result_rr_ops_ranking_drugs.mat'], 'result', 'predict');
             save([params.odir prefix 'result.mat'], 'result', 'predict');
    end
end


function [Prec, Predict] = single_file(filename, classify, byColumn, params, modelid)
    show_res=0;
    
    res = load([params.odir filename]);
    predY = res.Result.Y;

    [kFold, sLen, rLen, eLen] = size(predY);
    
    test = load([params.odir 'testY.mat']);
    testlY = test.test.lY;
    testclY = test.test.clY;
    
    %get prediction for each data size and epsilon value with kFold
    Predict = cell(kFold, sLen, rLen, eLen);
    Prec = cell(sLen, eLen, rLen);
    for sIndex = 1:sLen
        if(show_res)
            disp('----------------------------');
            disp(filename);
            disp(['ext size=',int2str(sIndex)]);
        end
        for eIndex = 1:eLen
            for rIndex = 1:rLen
                kprec = NaN(kFold, params.topK);
                for kIndex = 1:kFold
                    if(show_res)
                        disp(['fold=',int2str(kIndex)]);
                    end
                    pY = predY{kIndex, sIndex, rIndex, eIndex};
                    Y = testlY{kIndex};
                    if byColumn == 1
                        Y = testclY{kIndex};
                    end
                    
                    if(params.rmse)
                        [prec, Predict{kIndex, sIndex, rIndex, eIndex}] = calculate_rmse(pY, Y, classify, byColumn,params, modelid);
                        if(show_res)
                            disp(['prec ',num2str(prec)]); 
                        end
                    else
                        [prec, Predict{kIndex, sIndex, rIndex, eIndex}] = calculate_precision(pY, Y, classify, byColumn,params);
                    end
                    kprec(kIndex, :) = nanmean(prec);
                end 
                
                Prec{sIndex, eIndex, rIndex} = kprec;
            end
            
        end
    end
end

function [prec, pred]=calculate_rmse(pY, Y, classify, byColumn,params, modelid)
        [tSize, nCol] = size(pY);
        %prec = NaN(tSize, 10);
        prec=0;r2=params.r2;
        pred = NaN(tSize, nCol);
        if(modelid==4 && 0)
            load([params.odir,'stats.mat']);
            datamin=stats.y_min;
            datamax=stats.y_max;
            Y=preprocess_fm(Y,datamin,datamax);
        end
            
        if(params.corr && params.cindex==0)
            prec=0;
            
            if(params.rank_pat)
                for j=1:nCol
                    x=pY(:,j);
                    y=Y(:,j);
                    rho = corr(x,y,'type','Spearman','rows','complete');
                    prec=prec+rho;
                end
                prec=prec/nCol;
            else
                for j=1:tSize
                    x=pY(j,:)';
                    y=Y(j,:)';
                    rho = corr(x,y,'type','Spearman','rows','complete');
                    prec=prec+rho;
                end
                prec=prec/tSize;
            end
        elseif(params.cindex)
            % c-index
            Y=-Y;
            pY=-pY;
            load([params.odir,'/drug_stats.mat']);
            wpc=params.wpc;
            prec=0;
            prec_num=0; prec_den=0;
                for j=1:nCol 
                    rho=0;
                    x=pY(:,j);
                    
                    y=Y(:,j);
                    if(wpc)
                        sd=drug_stats.sd(j);
                        if(sd<0) sd=-sd; end
                        wd=drug_stats.wd(j);
                        if(wd<0) wd=-wd; end
                       % sd=1;
                        %wd=1;
                    else
                        sd=1;
                        wd=1;
                    end
                    
                    for i1=1:tSize
                        for i2=i1+1:tSize
                       
                            ci=0.5;
                            if(wpc)
                                if(~isnan(y(i1)+y(i2)))
                                    if(x(i1)>x(i2)) ci=0.5*(1+erf((y(i1)-y(i2))/(2*sd))); end
                                    if(x(i1)<x(i2)) ci=0.5*(1+erf((y(i2)-y(i1))/(2*sd))); end
                                end
                            else
                                if((y(i1)>y(i2)) && (x(i1)>x(i2))) ci=1; end
                                if((y(i1)<y(i2)) && (x(i1)<x(i2))) ci=1; end
                                if((y(i1)>y(i2)) && (x(i1)<x(i2))) ci=0; end
                                if((y(i1)<y(i2)) && (x(i1)>x(i2))) ci=0; end
                            end
                            
                            rho=rho+ci;
                        end
                    end
                    rho=(2/(tSize*(tSize-1)))*rho; 
                    
                    
                    prec_num=prec_num+wd*rho;
                    prec_den=prec_den+wd;
                end
                prec=prec_num/prec_den;
                
                %prec=prec/nCol;
        else
            if(r2==0)
                sz=pY./pY;
            else
                sz=Y.^2;
            end
            prec=nansum(nansum((pY-Y).^2))/(nansum(nansum(sz))); %tSize*nCol);
            if(r2)
                prec=1-prec;
            end
            %prec=sqrt(prec);
        end
        
        %%%%
        %if(prec>1) prec=1; end
       %%%%
       
       
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
        %{
        if byColumn == 1 
            % RMSE over patients averaged over drugs 
             for i = 1:nCol
               p = pY(:,i)';
               y = Y(:,i)';
               nan_index = isnan(y);
               p(nan_index) = NaN;
               
               if classify == 1
                   p = 1 - p;
               end
               pred(:,i)=p;
               
               dif=nansum((p-y).^2)/(tSize);
               prec=prec+sqrt(dif);
             end
             prec=prec/nCol;   
        else
            % RMSE over drugs averaged over patients
            for i = 1:tSize
               p = pY(i,:);
               y = Y(i, :);
               nan_index = isnan(y);
               p(nan_index) = NaN;
               if classify == 1
                   p = 1 - p;
               end
               pred(i,:)=p;
               
               
               dif=nansum((p-y).^2)/(nCol);
               prec=prec+sqrt(dif);
             end
             prec=prec/tSize;   
        end
        %}
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
end

function [prec, pred]= calculate_precision(pY, Y, classify, byColumn,params)
        [tSize, nCol] = size(pY);
        prec = NaN(tSize, params.topK);
        pred = NaN(tSize, nCol);
        if byColumn == 1
             for i = 1:nCol
               p = pY(:,i)';
               y = Y(:,i)';
               nan_index = isnan(y);
               p(nan_index) = NaN;
               
               if classify == 1
                   p = 1 - p;
               end
               ly = sum(nan_index == 0);
               
               if ly > params.topK   
                   p = real2Rank(p);
                   for k = 1:params.topK
                        n = find_k(p, y, k);
                        lp = label_topk(p, n);
                        if k == 5
                            pred(:,i)=lp;
                        end
                        
                        [pr, ~, ~] = other_metrics(y,  lp);
                        prec(i, k) = pr;
                   end
               else
                  n = find_k(p, y, params.topK/2);
                  pred(:,i) = label_topk(p, n);
               end 
            end

        else
            for i = 1:tSize
               p = pY(i,:);
               y = Y(i, :);
               nan_index = isnan(y);
               p(nan_index) = NaN;
               if classify == 1
                   p = 1 - p;
               end
               ly = sum(nan_index == 0);
               if ly > params.topK  
                   p = real2Rank(p);
                   for k = 1:params.topK
                        n = find_k(p, y, k);
                        lp = label_topk(p, n);
                        if k == 5 % why 5 ??
                            pred(i, :) = lp;
                        end
                        [pr, ~, ~] = other_metrics(y,  lp);
                        prec(i, k) = pr;
                   end   
               else
                  n = find_k(p, y, params.topK/2);
                  pred(i, :) = label_topk(p, n);
               end 
             end 
        end
end
function n = find_k(p, y, k)
    n = 1;
    while n < length(y)
         lp = label_topk(p, n);
         if sum(lp==1 & y==1) == k
             return;
         end
         n = n + 1;
    end
end
function [z] = label_topk(x, k)
    [~, rank_index] = sort(x);
    nanIndex = isnan(x);
    z = repmat(-1, [1, length(x)]);
    z(rank_index(1:k)) = 1;
    z(nanIndex) = NaN;
end
 
function [precision, recall, accuracy] = other_metrics(y, p)
	TP = sum(y== 1 & p == 1);
    recall = TP/sum(y==1);
	precision =  TP/sum(p == 1); 
	accuracy = sum(y==p)/length(y);	
end