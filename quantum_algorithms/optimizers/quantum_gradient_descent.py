import numpy as np
from tools import print_state

class QuantumGradientDescent:
    """
    Gradient Descent



    TODO: Add gradient averaging
    """
    def __init__(self,
                 method='',
                 step_length=1e-1, 
                 max_iter=200, # Number of iterations (2 func_evals per iter)
                 tol=1e-08, # Convergence tolerance
                 avg_num=0, # Number of gradients averaged, 0 if not.
                 feedback=True):
        self.method = method
        self.step_length = step_length
        self.max_iter = max_iter
        self.tol = tol
        self.feedback = feedback
        self.energies = []
        self.thetas = []
        self.grads = []
        self.avg_num = avg_num

    def __gradient_simple(self,theta):
        main_eval = self.L(new_theta) 
        param_evals = np.zeros(num_params)
        for k in range(num_params):
            copy = theta
            copy[k] += np.pi*0.5
            param_evals[k] = self.L(copy)
        grad = main_eval - param_evals
        return grad

    def __call__(self,theta):
        if self.method.lower() == 'adam':
            new_theta = self._adam(theta)
        elif self.method.lower() == 'adagrad':
            new_theta = self._adagrad(theta)
        elif self.method.lower() == 'rmsprop':
            new_theta = self._RMSprop(theta)
        else:
            new_theta = self._vanilla(theta)
        return new_theta

    def _gradient(self,theta):
        grad = np.zeros_like(theta)
        for k in range(len(theta)):
            e_k = np.zeros_like(theta)
            e_k[k] = 0.5*np.pi
            left = self.L(theta+e_k)
            right = self.L(theta-e_k)
            grad[k] = 0.5*(right - left)
        # Gradient averaging
        self.grads.append(grad)
        if len(self.grads) > self.avg_num and self.avg_num != 0:
            grad = np.mean(self.grads[-self.avg_num:],axis=0)
        return grad

    def _adam(self,theta):
        alpha = 1 #self.step_length
        eps = 1e-08
        beta_1 = 0.9
        beta_2 = 0.999
        M = np.zeros_like(theta)
        R = np.zeros_like(theta)
        m_hat = 1
        r_hat = 1

        new_theta = theta
        num_params = len(new_theta)
        for i in range(self.max_iter):
            t = i+1
            grad = self._gradient(new_theta)
            M = beta_1*M + (1 - beta_1)*grad
            R = beta_2*R + (1 - beta_2)*np.power(grad,2)
            m_hat = M / (1 - np.power(beta_1,t))
            r_hat = R / (1 - np.power(beta_2,t))
            new_theta += alpha*m_hat/(np.sqrt(r_hat) + eps)
            if self.feedback:
                E = self.L(new_theta)
                self.energies.append(E)
                print('Iteration {}: ⟨E⟩ = {}, t = {}'.format(i+1,E,new_theta))
        return new_theta

    def _vanilla(self,theta):
        alpha = self.step_length
        new_theta = theta
        num_params = len(new_theta)
        for i in range(self.max_iter):
            grad = self._gradient(new_theta)
            new_theta += grad*alpha
            if self.feedback:
                E = self.L(new_theta)
                self.energies.append(E)
                print('Iteration {}: ⟨E⟩ = {}, t = {}'.format(i+1,E,new_theta))
        return new_theta
    
    def _adagrad(self,theta):
        alpha = self.step_length
        eps = 1e-08
        new_theta = theta
        num_params = len(new_theta)
        cache = np.zeros_like(theta)
        for i in range(self.max_iter):
            grad = self._gradient(new_theta)
            cache += np.power(grad,2)
            new_theta += alpha*grad / (np.sqrt(cache)+eps)
            if self.feedback:
                E = self.L(new_theta)
                self.energies.append(E)
                print('Iteration {}: ⟨E⟩ = {}, t = {}'.format(i+1,E,new_theta))
        return new_theta

    def _RMSprop(self,theta):
        alpha = self.step_length
        eps = 1e-08
        gamma = 0.9
        new_theta = theta
        num_params = len(theta)
        cache = np.zeros_like(theta)
        last_theta = theta
        last_update = 100
        self.thetas.append(theta[0])
        for i in range(self.max_iter):
            grad = self._gradient(theta)
            cache = gamma*cache + (1. - gamma)*np.power(grad,2)
            update = alpha*grad / (np.sqrt(cache)+eps)
            theta += update
            self.thetas.append(theta[0])
            if self.feedback:
                E = self.L(theta)
                self.energies.append(E)
                print('Iteration {}: ⟨E⟩ = {}, t = {}'.format(i+1,E,theta))
            if all(np.abs(update)) < self.tol:
                print('Converged!')
                break
            #if last_update + update < np.abs(last_update) and\
            #    last_update + update < np.abs(update):
            #        print('Updated!',alpha)
            #        alpha *= 0.9
            last_update = update
        return new_theta
            
    def set_loss_function(self,loss_function):
        self.L = loss_function


