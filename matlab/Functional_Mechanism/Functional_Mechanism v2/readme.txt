This version is slightly updated from the description in the paper (Functional Mechanism in VLDB 2012), as we have found several heuristics that considerably enhance the performance of our algorithm:

1.  We normalize each attribute domain to [-1, 1] instead of [0, sqrt(1/d)]. For example, assume that 0 and 200 are the minimum and maximum values in the domain of AGE, then an age value 20 would be normalized to (20-100)/(200-100) = -0.8. Here 100 is the mean of 0 and 200.

2. The sensitivities of the cost functions are modified to accommodate the new normalization method.

3. The code is designed for the more general types of linear and logistic regressions that include an intercept b in the regression functions. For example, for linear regression, our code supports regression function of the form y = w*x + b (our paper only discusses the case when y = w*x).


NOTICE:  Training and testing data should be normalized (see the comments in the code)!