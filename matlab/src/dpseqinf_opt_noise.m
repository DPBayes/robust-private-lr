function bestnoise = dpseqinf_opt_noise(dp_eps, int_N, ext_N, ndims, repeats, noisevals)
% DPSEQINF_OPT_NOISE  Find optimal noise precision multiplier
% Returns multiplier for optimal noise precision for DP perturbed data
% for inferring the mean of ndims dimensional Gaussian.
%
% inputs:
% dp_eps: DP epsilon value
% int_N: internal data size
% ext_N: protected data size
% ndims: dimensionality
% repeats: number of repeats to average over
% noisevals: noise precision values to try (default: logspace(-6, 0, 61),
%   no need to change unless more precision is needed
%   or return value hits lower bound 1e-6)
%
% (c) 2016 Antti Honkela
%

  if nargin < 6,
    noisevals = logspace(-6, 0, 61);

    res = repeat_simulate_laplace(dp_eps, int_N, ext_N, ndims, repeats, noisevals);

    [~, bestindex] = min(res);
    bestnoise = noisevals(bestindex);
  end
end


function res = errorfun(x, y),
  res = sum(abs(x-y));
end


function r = mygauss(mu, lambda)
  r = struct('mu', mu, 'lambda', lambda);
end


function x = randlap(b, sz)
  u = rand(sz) - 0.5;
  x = b*sign(u) .* log(1-2*abs(u));
end


function posterior = gauss_posterior(prior, lambda_noise, obsmean, obs_N),
  post_lambda = prior.lambda + obs_N * lambda_noise;
  post_mu = (obs_N * lambda_noise * obsmean + prior.lambda * prior.mu) / post_lambda;
  posterior = struct('mu', post_mu, 'lambda', post_lambda);
end


function errors = simulate_laplace_mechanism(dp_eps, ext_mean, ext_N, intposterior, fullposterior, noisevals),
  dp_sens = 2*prod(size(ext_mean))/ext_N;
  errors = zeros(size(noisevals));

  dp_mean = ext_mean + randlap(dp_sens/dp_eps, size(ext_mean));
  for k=1:length(noisevals),
    posterior = gauss_posterior(intposterior, noisevals(k), dp_mean, ext_N);
    errors(k) = errorfun(fullposterior.mu, posterior.mu);
  end
end

function ave_errors = repeat_simulate_laplace(dp_eps, int_N, ext_N, ndims, repeats, noisevals)
  errors = zeros(repeats, length(noisevals));
  LAMBDA_NOISE = 1;
  truemean = 0;
  for k=1:repeats,
    int_mean = truemean + randn(ndims, 1) / sqrt(int_N);
    ext_mean = truemean + randn(ndims, 1) / sqrt(ext_N);
    prior = mygauss(zeros(ndims, 1), 0.1);
    intposterior = gauss_posterior(prior, LAMBDA_NOISE, int_mean, int_N);
    fullposterior = gauss_posterior(intposterior, LAMBDA_NOISE, ext_mean, ext_N);
    errors(k,:) = simulate_laplace_mechanism(dp_eps, ext_mean, ext_N, intposterior, fullposterior, noisevals);
  end
  ave_errors = mean(errors, 1)';
end
