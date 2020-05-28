import numpy as np

class SPSA:
    """

    Ref: 
    - https://www.jhuapl.edu/SPSA/PDF-SPSA/Spall_Implementation_Of_The_Simultaneous.pdf
    - https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=8005699&tag=1
    """
    def __init__(self,
                 max_iter=200,
                 min_change=0.1,
                 noise_var = 0.01,  # c
                 alpha=0.602,
                 gamma=0.101,
                 feedback=1,
                 grad_avg=0):
        """
        Simultaneous Perturbation Stochastic Approximation.

        Input:
            - loss_funcion: 
                The objective function to be minimized.
            - theta_init: 
                Vector with parameters to be optmimized.
            - max_iter: 
                Maximum number of iterations.
            - min_change: 
                Smallest initial change desired in 
                gain parameter 'a', used for theta parameters update.
            - noise_var: 
                Gain parameter 'c' for gradient approximation. 
                Approximate to the noise variance?
            - alpha: 
                Coefficient for determination of 'a'.
            - gamma: 
                Coefficient for determination of 'c'.
            - options:
                Dictionary for other options.
            - feedback (int):
                Number of iterations before printing loss evaluation.
            - grad_avg (int): 
                To implement gradient average. Number of previous gradients to be considerd.
        """
        self.gamma = gamma
        self.alpha = alpha
        self.A = int(max_iter * 0.1)
        self.c = noise_var
        self.min_change = min_change
        self.max_iter = max_iter

        self.feedback = feedback
        self.n_g = grad_avg

    def set_loss_function(self,loss_function):
        self.L = loss_function

    def __call__(self,theta):
        self.p = len(theta)
        self.a = self._step_length(theta)
        last_grad = [] 
        for k in range(self.max_iter):
            # Set gain variables
            a_k = self.a/np.power(self.A+k+1,self.alpha)
            c_k = self.c/np.power(k+1,self.gamma)

            # Evaluate gradient
            g = self._gradient(theta,c_k)

            # gradient averaging
            if k > self.n_g:
                for i in range(self.n_g):
                    g += last_grad[i]
                g /= self.n_g + 1
                last_grad.pop(0)
                last_grad.append(g)
            else:
                last_grad.append(g)
            
            # Update parameters
            theta = theta - a_k*g
            # Print loss evaluation.
            if k%self.feedback == 0:
                print('Iteration {}: \u2329E\u232A = {}'.format(k,self.L(theta)))
                print(theta,g)

    def _gradient(self,theta,c_k):
        # Sample from Bernoulli distribution
        delta = 2*np.random.binomial(n=1,p=0.5,size=self.p) - 1

        # Approximate gradient
        y_p = self.L(theta + c_k*delta)
        y_m = self.L(theta - c_k*delta)
        gradient = (y_p - y_m)/(2*c_k*delta)

        return gradient

    def _step_length(self,theta):
        g0 = self._gradient(theta,self.c)
        a0 = self.min_change*(np.power(self.A+1,self.alpha)\
                /np.sum(np.abs(g0)))
        return a0

