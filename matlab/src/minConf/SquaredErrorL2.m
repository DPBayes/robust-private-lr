function [f,g,H,T] = SquaredErrorL2(w,X,y,lambda)
% w(feature,1)
% X(instance,feature)
% y(instance,1)

[f,g,H,T] = SquaredError(w,X,y);
[f2,g2,H2,T2] = SquaredError(w,ones(1,length(w)),0);
f = f + (lambda/2)*f2;
g = g + g2;
H = H + H2;
T = T + T2;

end