function [drug_stats]=get_drug_stats(Y,params)
    Y=-Y; % take negative
    [n,m]=size(Y);
    for i=1:m
        y=Y(:,i);
        sd(i)=nanstd(y);
    end
    drug_stats.sd=sd;
    
    for i=1:m
        y=Y(:,i);
        
        for iter=1:1000
            g=rand(n,1);
            [sg,ig]=sort(g);
            pc(iter)=get_pcindex(y,ig,sd(i));
        end
        m=mean(pc);s=std(pc);
        [sy,iy]=sort(y,'descend');
        pcd=get_pcindex(y,iy,sd(i));
        wd(i)=(pcd-m)/s;
        disp(['done ',int2str(i)]);
    end
    drug_stats.wd=wd;
    save([params.odir,'drug_stats.mat'],'drug_stats');
end

function [fpc]=get_pcindex(y,o,sd)

    [n,n1]=size(y);
   pc=0;
    for i=1:n
       
        if(~isnan(y(i)))
            for j=i+1:n
                ci=0.5;
                    
                if(~isnan(y(j)))
                                
                    if(o(i)>o(j)) ci=0.5*(1+erf((y(j)-y(i))/(2*sd))); end
                    if(o(i)<o(j)) ci=0.5*(1+erf((y(i)-y(j))/(2*sd))); end
                end
                pc=pc+ci;
            end
            
        end
        
    end
    fpc=(2/(n*(n-1)))*pc; 
end

    
    