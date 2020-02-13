import numpy as np

class SPSA:
    def __init__(self,
                 loss_function,
                 theta_init,
                 max_iter=100,
                 step_length = 0.2,
                 noise_var = 0.01,
                 alpha=0.602,
                 gamma=0.101
                 options = {}):
        """
        Simultaneous Perturbation Stochastic Algorithm.

        Input:
            - loss_funcion : The objective function to be minimized.
            - theta_init   : Vector with parameters to be optmimized.
            - max_iter     : Maximum number of iterations.
            - step_length  : Gain parameter 'a' for theta parameters update.
            - noise_var    : Gain parameter 'c' for gradient approximation.
            - alpha        : Coefficient for determination of 'a'.
            - gamma        : Coefficient for determination of 'c'.
            - options      : Dictionary for other options.
                    -> 'feedback' (int) : Number of iterations before printing
                                          loss evaluation.
                    -> 'grad_avg' (int) : To implement gradient average.
                                          Number of previous gradients 
                                          to be considerd.
                    -> 'grad_bnds'(list): [lower,upper]
        
        Ref. 
            - https://www.jhuapl.edu/SPSA/PDF-SPSA/Spall_Implementation_Of_The_Simultaneous.pdf
            - https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=8005699&tag=1
        """
        self.L = loss_function
        self.theta = theta_init
        self.gamma = gamma
        self.alpha = alpha
        self.p = len(self.theta)
        self.A = int(max_iter * 0.1)
        self.a = step_length
        self.c = noise_var
        self.max_iter = max_iter
        self.options = options

        ### Options
        # Check feedback
        self.feedback = 1
        if 'feedback' in self.options:
            self.feedback = self.options['feedback']
        # Check if gradient averaging
        self.n_g = 0
        if 'grad_avg' in self.options:
            self.n_g = self.options['grad_avg']
        # Check if gradient bounds
        self.lb = 0
        self.ub = 0
        if 'grad_bnds' in self.options:
            self.lb = self.options['grad_bnds'][0]
            self.ub = self.options['grad_bnds'][1]

    def set_loss_function(self,loss_function):
        self.L = loss_function

    def run(self,prnt=True):
        theta = self.theta
        last_grad = [] 
        for k in range(self.max_iter):
            # Set gain variables
            a_k = self.a/np.power(self.A+k+1,self.alpha)
            c_k = self.c/np.power(k+1,self.gamma)

            # Sample from Bernoulli distribution
            delta = 2*np.random.binomial(n=1,p=0.5,size=self.p) - 1

            # Approximate gradient
            y_p = self.L(theta + c_k*delta)
            y_m = self.L(theta - c_k*delta)
            g = (y_p - y_m)/(2*c_k*delta)
            # gradient bounds
            if self.lb,self.ub != 0,0:
                for i,g_i in enumerate(g):
                    if g_i < 0:
                        sign = -1
                    else:
                        sign = 1
                    mag = np.abs(g_i)
                    if mag > self.ub:
                        g_i = sign*self.ub
                    elif mag < self.lb:
                        g_i = sign*self.lb
                    g[i] = g_i 
            # gradient averaging
            if k > self.n_g:
                for i for range(self.n_g):
                    g += last_grad[i]
                g /= self.n_g + 1
                last_grad.pop[0]
                last_grad.append(g)
            else:
                last_grad.append(g)

            # Update parameters
            theta = theta - a_k*g
            # Print loss evaluation.
            if k%self.feedback == 0:
                print('\u2329E\u232A =',self.L(theta))





