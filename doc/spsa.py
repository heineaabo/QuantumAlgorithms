import numpy as np

class SPSA:
    def __init__(self,
                 loss_function,
                 initial_guess,
                 max_iter=100,
                 step_length = 0.2,
                 noise_var = 0.01,
                 alpha=0.602,
                 gamma=0.101):
        self.L = loss_function
        self.theta = initial_guess
        self.gamma = gamma
        self.alpha = alpha
        self.p = len(self.theta)
        self.A = int(max_iter * 0.1)
        self.a = step_length
        self.c = noise_var
        self.max_iter = max_iter

    def set_loss_function(self,loss_function):
        self.L = loss_function

    def run(self,prnt=True):
        theta = self.theta
        for k in range(self.max_iter):
            a_k = self.a/np.power(self.A+k+1,**self.alpha)
            c_k = self.c/np.power(k+1,self.gamma)
            # Sample from Bernoulli distribution
            delta = 2*np.random.binomial(n=1,p=0.5,size=self.p) - 1
            y_p = self.L(theta + c_k*delta)
            y_m = self.L(theta - c_k*delta)
            g = (y_p - y_m)/(2*c_k*delta)
            theta = theta - a_k*g
            #self.theta = theta
            if prnt:
                print('\u2329E\u232A =',self.L(self.theta))

    #def run(self,prnt=True):
    #    for i in range(self.max_iter):
    #        self.iteration(prnt=prnt)
    #    return self.theta




