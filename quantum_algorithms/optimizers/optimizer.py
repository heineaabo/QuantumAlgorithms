
class Optimizer:
    def __init__(self):
        pass


class Minimizer:
    def __init__(self, method,
                       max_iter=200, # Minimizer iterations.
                       max_eval=200,  # Funtion evaluations.
                       tol=1e-08,
                       disp=True,
                       adapt=False):
        #super().__init__()
        self.max_iter = max_iter
        self.max_eval = max_eval
        self.tol = tol
        self.method = method
        self.disp = disp
        self.gradient = None

        # Set scipy.minimize options
        self.bounds = None
        max_eval_str = 'maxfev' # Different string notation for different methods
        if self.method == 'L-BFGS-B':
            self.bounds = [(0,2*np.pi) for i in theta]
            max_eval_str = 'maxfun'
        self.options = {'disp':self.disp,
                        'maxiter': self.max_iter,
                        max_eval_str: self.max_eval}
        if method.lower() == 'cobyla':
            del self.options[max_eval_str]
        if method.lower() == 'nelder-mead' and adapt == True:
            self.options['adaptive'] = True



    def __call__(self,theta):
        from scipy.optimize import minimize
        result = minimize(self.L,
                          theta,
                          bounds=self.bounds,
                          method=self.method,
                          options=self.options,
                          jac=self.gradient,
                          tol=self.tol)
        params = result.x
        return params


    def set_loss_function(self,loss_function):
        self.L = loss_function

    def set_gradient_function(self,gradient_function):
        self.gradient = gradient_function

    def set_option(self,option,value):
        if option == 'tol':
            self.tol = value
        else:
            if option == 'xtol' and self.method.lower() == 'nelder-mead':
                option = 'xatol'
            if option == 'ftol' and self.method.lower() == 'nelder-mead':
                option = 'fatol'
            if option == 'ftol' and self.method.lower() == 'cobyla':
                option = 'tol'
            self.options[option] = value

class QKSPSA:
    def __init__(self,
                   max_iter=200, # Minimizer iterations.
                   n_g=1, # Averaging number
                   ):
        from qiskit .aqua.components.optimizers import SPSA
        self.optimizer = SPSA(max_trials=max_iter,last_avg=n_g)

    def set_loss_function(self,loss_function):
        self.L = loss_function

    def __call__(self,theta):
        #import numpy as np
        params = self.optimizer.optimize(len(theta),self.L,initial_point=theta) 
        return params[0]


