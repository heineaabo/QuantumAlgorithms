import numpy as np
import qiskit as qk

class RY:
    """
    RY ansatz.

    Assumes the convention:
        ❘0⟩ - occupied orbital.
        ❘1⟩ - unoccupied orbital.

    Input:
            n (int): Number of particles.    
            l (int): Number of spin orbitals.
        occupied   : Convention for occupied orbital.
    """

    def __init__(self,n,l,occupied=1,pairing=False):
        self.n = n
        self.l = l
        self.occupied = occupied
        self.pairing = pairing

    def __call__(self,theta,qc,qb,qa=None,i=None):
        qc = self.prepare_Hartree_state(qc,qb)
        #theta = self.sort_theta(theta)
        if self.n == 1 and self.l == 2:
            qc = self._n1_l2(theta,qc,qb)
        elif self.n == 2 and self.l == 4:
            if self.pairing:
                qc = self._n1_l2(theta,qc,qb)
            else:
                qc = self._n2_l4(theta,qc,qb)
        else:
            print('WRONG INPUT')
        return qc

    def _n1_l2(self,theta,qc,qb):
        s = 1 # Second orbital
        if self.pairing:
            s = 2
        # Occupied
        qc.ry(theta[0],qb[0])

        # Singles
        qc.x(qb[0])
        if self.pairing:
            qc.cx(qb[0],qb[1])
            qc.cx(qb[0],qb[s+1])
        qc.cx(qb[0],qb[s])
        qc.x(qb[0])

        return qc

    def _n2_l4(self,theta,qc,qb):
        s1 = 1
        s2 = 2
        s3 = 3
        if self.pairing:
            s1 = 2
            s2 = 4
            s3 = 6
        # Occupied
        qc.ry(theta[0],qb[0])
        qc.rx(theta[1],qb[0])
        qc.ry(theta[2],qb[s1])
        qc.rx(theta[3],qb[s1])
        if self.pairing:
            qc.cx(qb[0],qb[1])
            qc.cx(qb[s1],qb[s1+1])

        # Singles
        qc.u3(theta[4]/2,0,0,qb[s2])
        #qc.u3(np.pi/4,0,0,qb[s2]) # For Hadamard
        qc.cx(qb[0],qb[s2])
        qc.cx(qb[s1],qb[s2])
        qc.u3(-theta[4]/2,0,0,qb[s2])
        #qc.u3(-np.pi/4,0,0,qb[s2]) # For Hadamard
        # Remove the below to have C-Hadamard instead of C-Ry
        qc.cx(qb[0],qb[s2])
        qc.cx(qb[s1],qb[s2])
        ######

        qc.cx(qb[0],qb[s3])
        qc.cx(qb[s1],qb[s3])

        if self.pairing:
            qc.cx(qb[s2],qb[s2+1])
            qc.cx(qb[s3],qb[s3+1])

        ## Doubles
        qc.ch(qb[0],qb[s2])
        qc.x(qb[s1])
        qc.cx(qb[s1],qb[s2])
        qc.x(qb[s1])
        qc.ch(qb[0],qb[s2])

        qc.cx(qb[s2],qb[s3])
        if self.pairing:
            qc.cx(qb[s2],qb[s2+1])
            qc.cx(qb[s3],qb[s3+1])
        return qc

    def new_parameters(self):
        """
        Define new parameters. 
        """
        num_params = 1
        if self.n == 2 and self.l == 4:
            #num_params = 3 
            num_params = 5
        #return 2*np.pi*np.random.randn(num_params)
        return np.zeros((num_params,))

    def prepare_Hartree_state(self,qc,qb):
        """
        Prepares Hartree-Fock state. Assume initialized to ❘0...0⟩ before 
        Example for n=2 and l=4 -> ❘1100⟩
        """
        # Reference state
        if self.occupied: #❘1...0...⟩
            for k in range(self.n):
                qc.x(qb[k])
        elif not self.occupied: #❘0...1...⟩
            for k in range(self.n,self.l):
                qc.x(qb[k])
        return qc

    def sort_theta(self,theta):
        """
        Sort parameters.
        """
        return theta


