import autograd.numpy as np
from autograd import grad, value_and_grad
import numpy as np_orig
import numpy.random as npr
import scipy.stats as sps
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.special import erfc


# > y := integrate(1/sqrt(2*Pi*a^2) * exp(\                                      
# > -(t-x)^2/(2*a^2)) * exp(-t/b) / (2*b), t=-infinity..0);
#                                                /
#                                                |
# y :=       lim         -1/4 csgn(a) exp(- x/b) |
#      t -> (-infinity)                          |
#                                                \

#                            2               1/2   2
#            2              a  - 2 b x      2    (a  + b t - b x)
#     csgn(a)  exp(x/b) exp(----------) erf(---------------------)
#                                 2                 2 a b
#                              2 b

#             1/2   2               2  \
#            2    (a  - b x)       a   |
#      - erf(---------------) exp(----)|/b
#                 2 a b              2 |
#                                 2 b  /

# > z := integrate(1/sqrt(2*Pi*a^2) * exp(\                                      
# > -(t-x)^2/(2*a^2)) * exp(t/b) / (2*b), t=0..infinity);    
#                                       2
#                                      a  + 2 b x
# z :=      lim       -1/4 csgn(a) exp(----------)
#      t -> infinity                         2
#                                         2 b

#     /     1/2   2                               1/2   2        \
#     |    2    (a  - b t + b x)         2       2    (a  + b x) |
#     |erf(---------------------) csgn(a)  - erf(---------------)|/b
#     \            2 a b                              2 a b      /

# > y + z assuming a > 0, b > 0;                                                 
#          1/2   2               2           2
#         2    (a  - b x)       a           a
#     erf(---------------) exp(----) + exp(----)
#              2 a b              2           2
#                              2 b         2 b
# 1/4 ------------------------------------------
#                     b exp(x/b)

#                           2           2                 1/2   2
#                          a           a                 2    (a  + b x)
#            exp(x/b) exp(----) + exp(----) exp(x/b) erf(---------------)
#                            2           2                    2 a b
#                         2 b         2 b
#      + 1/4 ------------------------------------------------------------
#                                         b

# > simplify(y + z assuming a > 0, b > 0);
#          2
#         a  - 2 b x
# 1/4 exp(----------)
#               2
#            2 b

#     /              1/2   2                           1/2   2            \
#     |    2 x      2    (a  + b x)        2 x        2    (a  - b x)     |
#     |exp(---) erf(---------------) + exp(---) + erf(---------------) + 1|/b
#     \     b            2 a b              b              2 a b          /


# > y := integrate(1/sqrt(2*Pi*a^2) * exp(\                                      
# > -(t-x)^2/(2*a^2)) * exp(t/b) / (2*b), t=-infinity..0);
#                                          2
#                                         a  + 2 b x
# y :=       lim         -1/4 csgn(a) exp(----------)
#      t -> (-infinity)                         2
#                                            2 b

#     /      1/2   2                               1/2   2        \
#     |     2    (a  - b t + b x)         2       2    (a  + b x) |
#     |-erf(---------------------) csgn(a)  + erf(---------------)|/b
#     \             2 a b                              2 a b      /

# > z := integrate(1/sqrt(2*Pi*a^2) * exp(\                                      
# > -(t-x)^2/(2*a^2)) * exp(-t/b) / (2*b), t=0..infinity);
# memory used=308.9MB, alloc=149.4MB, time=4.95
#                                            /
#                                            |
# z :=      lim       1/4 csgn(a) exp(- x/b) |
#      t -> infinity                         |
#                                            \

#                            2               1/2   2
#            2              a  - 2 b x      2    (a  + b t - b x)
#     csgn(a)  exp(x/b) exp(----------) erf(---------------------)
#                                 2                 2 a b
#                              2 b

#             1/2   2               2  \
#            2    (a  - b x)       a   |
#      - erf(---------------) exp(----)|/b
#                 2 a b              2 |
#                                 2 b  /

# > q := simplify(y + z assuming a > 0, b > 0);                                  
#                2
#               a  - 2 b x
# q := -1/4 exp(----------)
#                     2
#                  2 b

#     /              1/2   2                           1/2   2            \
#     |    2 x      2    (a  + b x)        2 x        2    (a  - b x)     |
#     |exp(---) erf(---------------) - exp(---) + erf(---------------) - 1|/b
#     \     b            2 a b              b              2 a b          /




def lsumpdf(x, a, b):
    return -np.log(4*b) + ((a**2 - 2*b*x) / (2 * b**2)) \
        + np.log(np.exp(2*x/b) * 
                 erfc((a**2 + b*x)/(np.sqrt(2) * a * b))
                 + erfc((a**2 - b*x)/(np.sqrt(2) * a * b)))
