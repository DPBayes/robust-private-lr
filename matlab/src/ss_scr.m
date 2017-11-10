function [c]=ss_scr(w,D,N)
epsilon = 1;
%D = 10;
%N = 500;
theta = 1; %1/0.1^2;
lambda = 1; %1/0.3^2;
X = randn(N,D);
%w = randn(D,1)./sqrt(lambda);
y = X*w + randn(N,1)./sqrt(theta);
Xt = randn(N,D);
yt = Xt*w + randn(N,1)./sqrt(theta);
% [mean(X,1) mean(y)]
% [std(X,[],1) std(y,[],1)]

cs=0.01:0.01:0.1;
choose_thr(epsilon, lambda, theta, cs, X, y, Xt, yt);

norr = sqrt(max(sum([X;Xt].^2,2)));
X = X./norr; Xt = Xt./norr;
norry = max(abs([y;yt]));
y = y./norry; yt = yt./norry;

c=choose_thr(epsilon, lambda, theta, cs, X, y, Xt, yt);
end