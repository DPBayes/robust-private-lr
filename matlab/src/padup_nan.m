function [X]=padup_nan(X,pad)
[r,c]=size(X);

for i=1:r
    for j=1:c
        x=X(i,j);
        if(isnan(x)==1)
	    if(pad==-1) r=rand; r=r*2; else r=pad; end	
            X(i,j)=r;
        end
    end
end
