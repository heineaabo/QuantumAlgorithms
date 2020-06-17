###############################################
#            FOR IBM 5-qubit computer         #
###############################################


import numpy as np
import qiskit as qk

class UCC5:
    """
    Unitary Coupled Cluster ansatz, with optimized singles and doubles circuits.
    Jordan-Wigner transformed operators.
    Assumes the convention:
        ❘0⟩ - occupied orbital.
        ❘1⟩ - unoccupied orbital.

    Input:
            n (int): Number of particles.    
            l (int): Number of spin orbitals.
        trunc (str): CC truncation.
                     - 'S' for singles.
                     - 'D' for doubles.
                     - 'SD' for singles and doubles.
        depth (int): Depth of Trotter approximation. 
        occupied   : Convention for occupied orbital.
    """

    def __init__(self,n,l,trunc='SD',depth=1,occupied=1):
        self.n = n
        self.l = l
        self.trunc = trunc
        self.depth = depth
        self.occupied = occupied
        self.mp2 = False # MP2 intial parameters

    def __call__(self,theta,qc,qb,qa=None,i=None):
        if 'S' not in self.trunc and 'D' not in self.trunc:
            raise ValueError('Need to spesify correct UCC type')
        if not isinstance(theta,(list,tuple,np.ndarray)):
            theta = [theta]
        qc = self.prepare_Hartree_state(qc,qb)
        theta = self.sort_theta(theta)
        #qc.u2(0,np.pi,qb[1])
        ##qc.cx(qb[1],qb[0])
        #qc.cx(qb[1],qb[2])
        ##qc.cx(qb[1],qb[3])
        #qc.u2(0,np.pi,qb[1])
        for rho in range(self.depth):
            if self.mp2:
                rho = 0
            if 'S' in self.trunc:
                qc = self.UCCS(qc,qb,theta[0],rho,qa=qa,d_j=i)
            if 'D' in self.trunc:
                qc = self.UCCD(qc,qb,theta[1],rho,qa=qa,d_j=i)
        return qc

    def new_parameters(self,h=[],v=[]):
        """
        Define new parameters. If one-body and two-body integrals 
        will return the MP2 cluster amplitudes. Else return random parameteres.
        """
        if len(h)!=0 and len(v)!=0:
            self.mp2 = True
        S = 0
        D = 0
        params = []
        if 'S' in self.trunc:
            for i in range(self.n):
                for a in range(self.n,self.l):
                    S += 1
                    if self.mp2:
                        params.append(0)

        if 'D' in self.trunc:
            for i in range(self.n):
                for j in range(i+1,self.n):
                    for a in range(self.n,self.l):
                        for b in range(a+1,self.l):
                            D += 1
                            if self.mp2:
                                t = (v[i,j,b,a]-v[i,j,a,b])/(h[i,i]+h[j,j]-h[a,a]-h[b,b])
                                params.append(t)
        self.num_S = S
        self.num_D = D
        if self.mp2:
            params = np.asarray(params)
            if all([x.imag == 0 for x in params]):
                params = params.real
            else:
                print('Imaginary parameters?')
            return np.asarray(params)
        else:
            return 2*np.pi*np.random.randn((self.num_S+self.num_D)*self.depth)

    def prepare_Hartree_state(self,qc,qb):
        """
        Prepares Hartree-Fock state. Assume initialized to ❘0...0⟩ before 
        Example for n=2 and l=4 -> ❘1100⟩
        """
        # Reference state
        if self.occupied: #❘1...0...⟩
            for k in range(self.n):
                #qc.x(qb[k])
                qc.u3(np.pi,0,np.pi,qb[k])
        elif not self.occupied: #❘0...1...⟩
            for k in range(self.n,self.l):
                #qc.x(qb[k])
                qc.u3(np.pi,0,np.pi,qb[k])
        return qc

    def sort_theta(self,theta):
        """
        Sort parameters to two channels for singles and doubles.
        TODO: 
            - Sort for depth > 1.
        """
        singles = self.num_S
        doubles = self.num_D
        thetaS = []
        thetaD = []
        if self.mp2:
            assert len(theta) == (singles+doubles),'WRONG LENGTH OF {} != {}'.format(theta,singles+doubles)

            thetaS.append(theta[:singles])
            thetaD.append(theta[singles:])
        else:
            assert len(theta) == (singles+doubles)*self.depth
            for i in range(self.depth):
                thetaS.append(theta[(i*singles):(i+1)*singles])
                thetaD.append(theta[(self.depth*singles + (i)*doubles):\
                                    (self.depth*singles + (i+1)*doubles)])
        new_theta = [thetaS,thetaD]
        return new_theta

    def UCCS(self,qc,qb,theta,rho,qa=None,d_j=-1):
        """
        Optimized UCC singles implementation.
        """
        # UCCS
        return qc

    def UCCD(self,qc,qb,theta,rho,qa=None,d_j=-1):
        """
        Reduced UCC doubles implementation.
        qa -> Ancilla qubit register
        d_j -> parameter index of parameter to diffrentiate
        """
        # UCCD
        # YXXX
        qc.u3(-np.pi/2,-np.pi/2,np.pi/2,qb[0]) # Rx(-pi/2)
        #qc.rx(np.pi/2,qb[0]) # Rx(-pi/2)
        qc.u2(0,np.pi,qb[1]) # H gate 
        #qc.h(qb[1])
        qc.u2(0,np.pi,qb[2]) # H gate
        #qc.h(qb[2])
        qc.u2(0,np.pi,qb[3]) # H gate
        #qc.h(qb[3])

        qc.cx(qb[0],qb[1])
        qc.cx(qb[2],qb[1])
        qc.cx(qb[3],qb[1])

        qc.u1(-theta[rho][0]/8,qb[1]) # Rz
        #qc.rz(-theta[rho][0]/8,qb[1]) # Rz

        qc.cx(qb[3],qb[1])
        qc.cx(qb[2],qb[1])
        qc.cx(qb[0],qb[1])

        qc.u3(np.pi/2,-np.pi/2,np.pi/2,qb[0]) # Rx(pi/2)
        #qc.rx(-np.pi/2,qb[0]) # Rx(-pi/2)
        qc.u2(0,np.pi,qb[1]) # H gate 
        #qc.h(qb[1])
        qc.u2(0,np.pi,qb[2]) # H gate
        #qc.h(qb[2])
        qc.u2(0,np.pi,qb[3]) # H gate
        #qc.h(qb[3])

        return qc
