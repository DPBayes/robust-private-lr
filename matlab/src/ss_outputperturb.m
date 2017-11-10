epsilon = 1;
addpath(genpath(pwd))
D = 10;
N = 1000;
X = randn(N,D);
w = randn(D,1).*(rand(D,1) > .5);
y = X*w + randn(N,1)*0.05;
X = X./(max(sqrt(sum(X.^2,2))));
y = y./(max(y));

w_est = LR_output_perturbation(10000,X,y);
w_dp = LR_output_perturbation(epsilon,X,y);

figure;
subplot(1,2,1);
imagesc(w_est);colormap gray;
title('Original Weights');
subplot(1,2,2);
imagesc(w_dp);colormap gray;
title('DP Weights');

figure
y_est = X*w_est;
plot(y_est(1:50))
hold on
y_dp = X*w_dp;
plot(y_dp(1:50),'r')
