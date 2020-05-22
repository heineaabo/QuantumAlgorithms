import numpy as np
from tools import print_state

class QuantumGradientDescent:
    """
    Gradient Descent



    TODO: Add gradient averaging
    """
    def __init__(self,loss_function,theta,step_length=1e-1,max_iter=200,tol=1e-08,feedback=True):
        self.L = loss_function
        self.theta = theta
        self.step_length = step_length
        self.max_iter = max_iter
        self.tol = tol
        self.feedback = feedback

    def __call__simple(self):
        new_theta = self.theta
        num_params = len(new_theta)
        last_eval = 0
        for i in range(self.max_iter):
            main_eval = self.L(new_theta)
            # Check convergence
            if i > 0:
                if np.abs(main_eval,last_eval) < tol:
                    break
            param_evals = np.zeros(num_params)
            for k in range(num_params):
                copy = new_theta
                copy[k] += np.pi/2
                param_evals[k] = self.L(copy)
            grads = main_eval - param_evals
            new_theta -= grads*step_length
            last_eval = main_eval
        return new_theta

    def __call__(self):
        new_theta = self.theta
        num_params = len(new_theta)
        for i in range(self.max_iter):
            grads = np.zeros(num_params)
            for k in range(num_params):
                e_k = np.zeros(num_params)
                e_k[k] = 1
                e_k *= 0.5*np.pi
                left = self.L(new_theta+e_k)
                right = self.L(new_theta-e_k)
                print(left,right)
                grads[k] = 0.5*(left - right)
                #grads[k] = 0.5*(self.L(new_theta+e_k) - self.L(new_theta-e_k))
            new_theta -= grads*self.step_length
            if self.feedback:
                E = self.L(new_theta)
                print('Iteration {}\n⟨E⟩ = {}\n{}'.format(i+1,E,new_theta))
        return new_theta
            
