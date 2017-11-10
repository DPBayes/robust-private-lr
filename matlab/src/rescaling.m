function [Z] = rescaling(X, byrow)
        %default, X will be escaled by column
	if ~exist('byrow', 'var')
		byrow = 0;
	end
	[r, c] = size(X);
	if(byrow)
		Z = X - nanmin(X, [], 2)*ones(1, c);
		Z = Z./(nanmax(Z, [], 2)*ones(1, c));
	else
		Z = X - ones(r, 1)*nanmin(X, [], 1);
		Z = Z./(ones(r, 1)*nanmax(Z, [], 1));
	end
	
