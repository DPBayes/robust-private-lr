import autograd.numpy as np
from autograd import grad, value_and_grad
from autograd.scipy.special import erfc
import numpy as np_orig
import numpy.random as npr
import scipy.stats as sps
import matplotlib.pyplot as plt
from scipy.optimize import minimize


def generate_data(D_x, D_y, N_int, N_ext, beta = [], sigma=0.1, DP_eps = 1.0):
    if len(beta) == 0:
        beta = npr.randn(D_x, D_y)
    x_int = npr.randn(N_int, D_x)
    x_ext = npr.randn(N_ext, D_x)
    y_int = x_int.dot(beta) + sigma*npr.randn(N_int, D_y)
    y_ext = x_ext.dot(beta) + sigma*npr.randn(N_ext, D_y)
    xx_int = x_int.T.dot(x_int)/N_int
    xx_ext = x_ext.T.dot(x_ext)/N_ext
    xy_int = x_int.T.dot(y_int)/N_int
    xy_ext = x_ext.T.dot(y_ext)/N_ext
    DP_xx_scale = 1. / (DP_eps / 2.)
    DP_xy_scale = 2. * D_x**2 / (DP_eps / 2.)
    xx_ext_DP = xx_ext + sps.wishart.rvs(df=D_x+1, scale=np.eye(D_x)*DP_xx_scale/N_ext)
    xy_ext_DP = xy_ext + npr.laplace(scale=DP_xy_scale/N_ext, size=(D_x, D_y))
    return ({'N_ext': N_ext, 'N_int': N_int, 'D_x': D_x,'D_y': D_y,
             'DP_xx_scale': DP_xx_scale, 'DP_xy_scale': DP_xy_scale,
             'xx_int': xx_int, 'xx_ext_DP': xx_ext_DP,
             'xy_int': xy_int, 'xy_ext_DP': xy_ext_DP,},
            {'x_int': x_int, 'x_ext': x_ext,
             'y_int': y_int, 'y_ext': y_ext,
             'xx_ext': xx_ext, 'xy_ext': xy_ext,
             'beta':beta})

def lgausspdf(x, mu, sigma):
    return -0.5 * np.sum((x-mu)**2 / sigma**2)

def llaplacepdf(x, sigma):
    return -np.sum(np.abs(x) / sigma)

def lsumpdf(x, a, b):
    return np.sum(np.log(4*b) + ((a**2 - 2*b*x) / (2 * b**2)) \
                  + np.log(np.exp(2*x/b) * 
                           erfc((a**2 + b*x)/(np.sqrt(2) * a * b))
                           + erfc((a**2 - b*x)/(np.sqrt(2) * a * b))))

def xy_true_logl(data, xy):
    return lgausspdf(np.reshape(data['xy_int'], -1),
                     xy, np.sqrt(1.0/data['N_int'])) + \
        lsumpdf(np.reshape(data['xy_ext_DP'], -1) - xy,
                1/data['N_ext'],
                data['DP_xy_scale']/data['N_ext'])

def xy_logl(data, xy):
    return lgausspdf(np.reshape(data['xy_int'], -1),
                     xy, np.sqrt(1.0/data['N_int'] + 1.0/data['N_ext'])) + \
        llaplacepdf(np.reshape(data['xy_ext_DP'], -1) - xy,
                    data['DP_xy_scale']/data['N_ext'])

def xy_logl_1d(data, xy, dim):
    return lgausspdf(np.reshape(data['xy_int'], -1)[dim],
                     xy,
                     np.sqrt(1.0/data['N_int'] + 1.0/data['N_ext'])) + \
        llaplacepdf(np.reshape(data['xy_ext_DP'], -1)[dim] - xy,
                    data['DP_xy_scale']/data['N_ext'])

def optimize_xy(data):
    D_x = data['D_x']
    D_y = data['D_y']
    val = np.zeros(D_x*D_y)
    allfits = []
    for dim in range(D_x*D_y):
        objective = lambda x: -xy_logl_1d(data, x, dim)
        init_params = (np.reshape(data['xy_ext_DP'], -1) \
                       + 0.1 * (np.reshape(data['xy_int'] - data['xy_ext_DP'], -1)))[dim]
        #init_params = np.reshape(data['xy_int'], -1)
        opt_params = minimize(value_and_grad(objective), init_params, jac=True, method='BFGS')
        if not opt_params['success']:
            opt_params = minimize(objective, init_params, method='Nelder-Mead')
        allfits.append(opt_params)
        val[dim] = opt_params['x']
    return val, allfits

def evaluate_solution(data, rest, opt_params):
    D_x = data['D_x']
    D_y = data['D_y']
    N_int = data['N_int']
    N_ext = data['N_ext']
    xy_fit = N_int * data['xy_int'] + N_ext * np.reshape(opt_params, (D_x, D_y))
    xx = N_int * data['xx_int'] + N_ext * data['xx_ext_DP'] + np.eye(D_x)
    xy = N_int * data['xy_int'] + N_ext * data['xy_ext_DP']
    nxx_int = N_int * data['xx_int'] + np.eye(D_x)
    nxy_int = N_int * data['xy_int']
    beta1 = np.linalg.solve(xx, xy)
    beta2 = np.linalg.solve(xx, xy_fit)
    beta3 = np.linalg.solve(nxx_int, nxy_int)
    print('raw:', np.sum((rest['beta'] - beta1)**2))
    print('new:', np.sum((rest['beta'] - beta2)**2))
    print('int:', np.sum((rest['beta'] - beta3)**2))

def test_me(D_x, D_y, N_int, N_ext):
    data, rest = generate_data(D_x, D_y, N_int, N_ext)
    #data, rest = generate_data(2, 1, 20, 200)
    opt_params, all_fits = optimize_xy(data)
    evaluate_solution(data, rest, opt_params)
    return data, rest, opt_params, all_fits


# model {
#   xx_int ~ wishart(N_int, xx/N_int);
#   xx_ext ~ wishart(N_ext, xx/N_ext);
#   for (i in 1:D_x) {
#     xy_int[i,] ~ normal(xy[i,], sqrt(1.0/N_int));
#     xy_ext[i,] ~ normal(xy[i,], sqrt(1.0/N_ext));
#   }
#   xx_perturb ~ wishart(D_x+1, diag_matrix(rep_vector(1.0, D_x))*DP_xx_scale/N_ext);
#   #xx_ext_DP = xx_ext + xx_perturb;
#   #increment_log_prob(wishart_log(xx_ext_DP - xx_ext, 
#   #				 D_x+1, diag_matrix(rep_vector(1.0, D_x))*DP_xx_scale/N_ext));
#   for (i in 1:D_x) {
#     xx_ext_DP[i,] ~ normal(xx_perturb[i,] + xx_ext[i,], 1e-6);
#     xy_ext_DP[i,] ~ double_exponential(xy_ext[i,], DP_xy_scale/N_ext);
#   }
# }


#def get_mean(fit, param):
#    return np.mean(fit.extract(param)[param], 0)

#dpseqinf_dat, rest = generate_data(10, 10, 20, 200)

#fit = pystan.stan(file='dpseqinf_lr.stan', data=dpseqinf_dat, iter=1000, chains=4)

#mu = get_mean(fit, 'mu')
#stanextmean = get_mean(fit, 'extmean')
