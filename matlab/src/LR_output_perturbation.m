function w_dp = LR_output_perturbation(epsilon,X,y)
% Differentially private linear regression with L2 regularization
% based on output perturbation
%
% epsilon: privacy budget
% X: training input set, N (samples) x D (dimensions), x_{nd} \in [-1,1]
% y: training output set, N x 1, y_n \in [-1,1]
% w_dp: weights to be returned differentialy privately

[N,D] = size(X);
R = 1;
lambda = sqrt(D/(N*epsilon)); % R and lambda are chose in a data-independent manner

funObj = @(w)SquaredErrorL2(w,X,y,lambda);
funProj = @(w)sign(w).*projectRandom2C(w.^2,R);
w_est = minConf_PQN(funObj,zeros(D,1),funProj);

temp_norm = gamrnd(D,(12*R+8)/(lambda*N*epsilon));
temp_angle = rand(D,1)-0.5;
kappa = temp_norm*temp_angle./sqrt(sum(temp_angle.^2));
w_dp = w_est + kappa;

end

