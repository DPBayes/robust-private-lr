function analyze_data(data,params)

    X=data.X;
    Y=data.Y;
    [n,d]=size(X);
    [n,m]=size(Y);

    if(1)
        n=500;
        find_c(d,n,X,Y);
    end

    if(0)

        X = L2normalize(X,params.l2norm_dim,params);
        Y=ynormalize(Y);

        s=1;
        dr1=d*(d+1)/2;
        nr1=sum(sum(X'*X));
        rmax=0;

        cs1=0.01:0.01:0.1;
        cs2=1:1:1e1;
        cs2=[0.1 0.5 cs2];
        lenC1=length(cs1);
        lenC2=length(cs2);


        for i=1:lenC1
            c1=cs1(i);
            X=clamp_data(X,c);

            for j=1:lenC2
                c2=cs2(j);
                nr2=sum(sum(X'*X));
                dr2=c*d*(d+1)/2;

                r=(nr2/dr2)/(nr1/dr1);
                if(rmax<r) rmax=r; cmax=c; end
                R(i)=r;
            end
            rmax
            cmax
            plot(R,'LineWidth',2);  
            xlabel('Clipping threshold', 'fontsize', 16);
            ylabel('Effective gain', 'fontsize', 16);
        end

    end
end

function Z=ynormalize(Y)
    [n, m]=size(Y);
    Z=Y-nanmean(nanmean(Y));
    %Z=Z./sqrt(ones(n,1)*nansum(Z.^2,1));
   % Z=Z/n;
end