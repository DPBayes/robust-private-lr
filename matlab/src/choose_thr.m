function c = choose_thr(epsilon, lambda, theta, cs, X, y, Xt, yt)
% chooses the best clipping threshold among the list cs
%
% epsilon: differential privacy budget
% lambda: assumed prior precison on the weights w ~ N(0,I/lambda)
% theta: assumed prior precison on the output y ~ N(Xw,I/theta)
% cs: list of candidate thresholds, data clipped to [-c,c]
% X: training input set N (samples) x D (dimensions)
% y: training output set N (samples) x 1
% X: training input set N (samples) x D (dimensions)
% y: training output set N (samples) x 1
% Xt and yt: (optional) test set
if nargin<7, Xt = X; yt = y; end
num_trials = 10;
[N,D] = size(X);
errs = zeros(length(cs),1);
errst = zeros(length(cs),1);
errs_dp = zeros(length(cs),1);
errst_dp = zeros(length(cs),1);
for i = 1:length(cs)
  c = cs(i);
  X1 = min(X,repmat(c,N,D)); X1 = max(X1,repmat(-c,N,D));
  y1 = min(y,repmat(c,N,1)); y1 = max(y1,repmat(-c,N,1));
  w1 = (theta*(X1'*X1)+lambda*eye(D))\(theta*X1'*y1);
  errs(i) = nanmean((X*w1-y).^2);  
  errst(i) = nanmean((Xt*w1-yt).^2);  
  
  for j = 1:num_trials
    XX1 = gamrnd(1,(2*c^2*(D*(D+1)/2+D))./epsilon, [D D]);
    XX1 = XX1.*(2*(rand(D)>0.5)-1);
    XX1 = triu(XX1) + triu(XX1)' + diag(diag(XX1));
    Xy1 = gamrnd(1,(2*c^2*(D*(D+1)/2+D))./epsilon, [D 1]);
    Xy1 = Xy1.*(2*(rand(D,1)>0.5)-1);
    w2 = (theta*(X1'*X1)+XX1+lambda*eye(D))\(theta*X1'*y1 + Xy1);
    errs_dp(i) = errs_dp(i) + nanmean((X*w2-y).^2)/num_trials;
    errst_dp(i) = errst_dp(i) + nanmean((Xt*w2-yt).^2)/num_trials;
  end
end
[~,ind] = min(errst_dp);
c = cs(ind);

%{
figure
plot(cs,errs)
hold on
plot(cs,errs_dp,'r')
title('Training set')
ylabel('MSE')
xlabel('clipping threshold')
ylim([0 2*max(errst)])
legend('NP','DP','Location','NorthWest')

figure
plot(cs,errst)
hold on
plot(cs,errst_dp,'r')
title('Test set')
ylabel('MSE')
xlabel('clipping threshold')
ylim([0 2*max(errs)])
legend('NP','DP','Location','NorthWest')
%}