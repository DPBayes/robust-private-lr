function [Zx, Zy] = ignore_nan(X, y)
	valid_index = ~isnan(y);
	Zx = X(valid_index, :);
	Zy = y(valid_index);
end