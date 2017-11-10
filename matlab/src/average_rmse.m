function average_rmse(params,action,byColumn)
    
    prefix = '';
    
    if params.real2rank == 1
        prefix = 'rank_';
    end
    
    if params.rmse
        prefix='rmse_';
    end
    load([params.odir prefix 'result.mat']);
    
  % action
   
    if(action>0)
        load([params.odir prefix 'SumResult.mat']);
    else
        Sresult = struct();
        Sresult.orr = result.orr;
        if(params.no_models>1)
            Sresult.rr = result.rr;
            Sresult.irr = result.irr;
            Sresult.dpvb = result.dpvb;
           % Sresult.fm = result.fm;
        end
        %predict = struct();
    end
    
    %{
    in1=Sresult.rr;
    in2=result.rr;
    out=sumup(in1,in2,params,action);
    Sresult.rr=out;
   %}
    
    [Sresult.orr]=sumup(Sresult.orr,result.orr,params,action);
    if(params.no_models>1)
        [Sresult.rr]=sumup(Sresult.rr,result.rr,params,action);
        [Sresult.irr]=sumup(Sresult.irr,result.irr,params,action);
        [Sresult.dpvb]=sumup(Sresult.dpvb,result.dpvb,params,action);
        %[Sresult.fm]=sumup(Sresult.fm,result.fm,params,action);
    end

    
    
    if byColumn == 1
            %save([params.odir 'result_rr_ops_ranking_patients.mat'], 'result', 'predict');
            save([params.odir prefix 'SumResult.mat'], 'result', 'predict');
    else
             %save([params.odir 'result_rr_ops_ranking_drugs.mat'], 'result', 'predict');
                save([params.odir prefix 'SumResult.mat'], 'Sresult', 'predict');
    end
end


function [SPrec] = sumup(SPrec, Prec, params, action)
    show_res=0;
    
    
    %get prediction for each data size and epsilon value with kFold
    [kFold, sLen, rLen, eLen] = size(Prec);
    for sIndex = 1:sLen
        if(show_res)
            disp('----------------------------');
            disp(['ext size=',int2str(sIndex)]);
        end
        for eIndex = 1:eLen
            for rIndex = 1:rLen
                kprec = NaN(kFold, params.topK);
                
                for kIndex = 1:kFold
                    if(show_res)
                        disp(['fold=',int2str(kIndex)]);
                    end
                    
                    %{
                    pY = predY{kIndex, sIndex, rIndex, eIndex};
                    Y = testlY{kIndex};
                    if byColumn == 1
                        Y = testclY{kIndex};
                    end
                    
                    if(params.rmse)
                        [prec, Predict{kIndex, sIndex, rIndex, eIndex}] = calculate_rmse(pY, Y, classify, byColumn,params);
                        if(show_res)
                            disp(['prec ',num2str(prec)]); 
                        end
                    else
                        [prec, Predict{kIndex, sIndex, rIndex, eIndex}] = calculate_precision(pY, Y, classify, byColumn,params);
                    end
                    kprec(kIndex, :) = nanmean(prec);
                    %}
                    
                end 
                 kprec = Prec{sIndex, eIndex, rIndex};
                 skprec = SPrec{sIndex, eIndex, rIndex};
                 
                 if(action==0)
                    SPrec{sIndex, eIndex, rIndex} = kprec;
                 end
                 if(action==1)
                    SPrec{sIndex, eIndex, rIndex} = skprec + kprec;
                 end
                 if(action==2)
                    SPrec{sIndex, eIndex, rIndex} = skprec/params.ITER;
                 end
            end
           
        end
    end
end
