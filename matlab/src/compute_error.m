function compute_error(params,info)
    result = struct();
    result.lr=get_model_error('lr.mat',params,info);
    
    result.rplr=get_model_error('rplr.mat',params,info);
    if(info.sch<0)
        result.rnplr=get_model_error('rnplr.mat',params,info);
        result.plr=get_model_error('plr.mat',params,info);
        result.oplr=get_model_error('oplr.mat',params,info);
        result.fmlr=get_model_error('fmlr.mat',params,info);
    end
    save([info.odir 'result.mat'], 'result');
end

function [Result]=get_model_error(filename,params,info)
    res = load([info.odir filename]);
    predY = res.Result.Y;

    [kFold, sLen] = size(predY);
    Result=cell(sLen,kFold);
    test = load([info.odir 'testY.mat']);
    testlY = test.test.lY;
    
    for sIndex = 1:sLen
        for kIndex = 1:kFold
                 pY = predY{kIndex, sIndex};
                 Y = testlY{kIndex};
                 result=get_error(pY,Y,params,info);
                 Result{sIndex,kIndex}=result;
                
        end
         
    end
    
end
function [result]=get_error(pY, Y, params,info)
        [tSize, nCol] = size(pY);
        %r2=params.r2;
        %{
        if(modelid==4 && 0)
            load([params.odir,'stats.mat']);
            datamin=stats.y_min;
            datamax=stats.y_max;
            Y=preprocess_fm(Y,datamin,datamax);
        end
        %}    
        if(params.corr && info.cindex==0)
            result=0;
            
            if(params.rank_pat)
                for j=1:nCol
                    x=pY(:,j);
                    y=Y(:,j);
                    rho = corr(x,y,'type','Spearman','rows','complete');
                    result=result+rho;
                end
                result=result/nCol;
            else
                for j=1:tSize
                    x=pY(j,:)';
                    y=Y(j,:)';
                    rho = corr(x,y,'type','Spearman','rows','complete');
                    result=result+rho;
                end
                result=result/tSize;
            end
        elseif(info.cindex)
            % c-index
            Y=-Y;
            pY=-pY;
            load([info.odir,'/drug_stats.mat']);
            wpc=params.wpc;
            result=0;
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
                result=prec_num/prec_den;
                
                %prec=prec/nCol;
       %{
        else
            if(r2==0)
                sz=pY./pY;
            else
                sz=Y.^2;
            end
            result=nansum(nansum((pY-Y).^2))/(nansum(nansum(sz))); %tSize*nCol);
            if(r2)
                result=1-result;
            end
            %prec=sqrt(prec);
         %}
        end
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
end


