
class Optimizer:
    def __init__(self):
        pass


class Minimizer:
    def __init__(self, method,
                       max_iters=200, # Minimizer iterations.
                       max_evals=200,  # Funtion evaluations.
                       tol=1e-08):
        super().__init__()
        self.max_iters = max_iters
        self.max_evals = max_evals
        self.tol = tol
        self.method = method

        # Set scipy.minimize options
        self.bounds = None
        max_eval_str = 'maxfev' # Different string notation for different methods
        if self.method == 'L-BFGS-B':
            self.bounds = [(0,2*np.pi) for i in theta]
            max_eval_str = 'maxfun'
        self.options = {'disp':True,
                        'maxiter': max_iters,
                        max_eval_str: self.max_evals}
        #if method == 'Nelder-Mead':
       #     if len(theta) > 5:
       #         self.options['adaptive'] = True



    def __call__(self,theta):
        from scipy.optimize import minimize
        result = minimize(self.L,
                          theta,
                          bounds=self.bounds,
                          method=self.method,
                          options=self.options,
                          tol=self.tol)
        params = result.x
        return params


    def set_loss_function(self,loss_function):
        self.L = loss_function

    #def set_theta(self,theta):
    #    self.theta = theta

