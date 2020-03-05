import numpy as np
import qiskit as qk

def UCCS(system,qc,qb,theta,depth=1):
    n_qubits = system[0]
    n_fermi = system[1]
    theta = theta.reshape(n_fermi,n_qubits-n_fermi)

    # Prepare reference state
    for k in range(n_fermi,n_qubits):
        qc.x(qb[k])

    # UCCS
    for rho in range(depth):
        for i in range(n_fermi):
            for a in range(n_fermi,n_qubits):
                # Pauli-z
                for k in range(i+1,a):
                    qc.z(qb[k])

    return qc

def UCCD(system,qc,qb,theta,depth=1):
    n_qubits = system[0]
    n_fermi = system[1]
    #theta = theta.reshape(n_fermi,n_fermi,n_qubits-n_fermi,n_qubits-n_fermi)

    # Prepare reference state
    for k in range(n_fermi,n_qubits):
        qc.x(qb[k])

    t = 0
    # UCCD
    for rho in range(depth):
        for i in range(n_fermi):
            for j in range(i+1,n_fermi):
                for a in range(n_fermi,n_qubits):
                    for b in range(a+1,n_qubits):
                        # Pauli-z
                        for k in range(i+1,j):
                            qc.z(qb[k])
                        for l in range(a+1,b):
                            qc.z(qb[l])
                        
                        t += 1

    return qc

def sort_theta(system,theta):
    l = system[0] # n_qubits
    n = system[1] #n_fermi
    singles = n*(l-n)
    #assert len(theta) == singles+doubles
    new_theta = [theta[:singles],theta[singles:]]
    return new_theta

def UCC(system,qc,qb,theta,terms='SD',depth=1):
    if terms == 'S':
        qc = UCCS(system,qc,qb,theta,depth=depth)
    elif terms == 'D':
        qc = UCCD(system,qc,qb,theta,depth=depth)
    elif terms == 'SD':
        theta = sort_theta(system,theta)
        qc = UCCS(system,qc,qb,theta[0],depth=depth)
        qc = UCCD(system,qc,qb,theta[1],depth=depth)
    else:
        raise ValueError('Need to spesify correct UCC type')
    return qc

class UnitaryCoupledCluster:

    def __init__(self,n,l,truncation,depth):
        self.n = n
        self.l = l
        self.trunc = truncation
        self.depth = depth

    def __call__(self,qc,qb,theta):
        """
        system: [n_qubits,n_fermi]
        """
        system = [self.l,self.n]
        if self.trunc == 'S':
            qc = UCCS(system,qc,qb,theta,depth=self.depth)
        elif self.trunc == 'D':
            qc = UCCD(system,qc,qb,theta,depth=self.depth)
        elif self.trunc == 'SD':
            theta = sort_theta(system,theta)
            qc = UCCS(system,qc,qb,theta[0],depth=self.depth)
            qc = UCCD(system,qc,qb,theta[1],depth=self.depth)
        else:
            raise ValueError('Need to spesify correct UCC type')
        return qc

    def get_num_parameters(self):
        S = self.n*(self.l-self.n)
        D = 0
        for i in range(self.n):
            for j in range(i+1,self.n):
                for a in range(self.n,self.l):
                    for b in range(a+1,self.l):
                        D += 1
        return S, D

