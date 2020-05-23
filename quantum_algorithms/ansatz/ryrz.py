import numpy as np
import qiskit as qk

class RYRZ:
    """
    """

    def __init__(self,n,l,depth=1,occupied=1,pair=True):
        self.n = n
        self.l = l
        self.depth = depth
        self.occupied = occupied
        self.pair = pair # If particles should be in pairs

    def __call__(self,theta,qc,qb,qa=None,i=None):
        """
        Pairing ansatz. For 2 particles and 4 orbitals.
        """
        theta = self.sort_theta(theta)
        qc = self.prepare_Hartree_state(qc,qb)
        #for d in range(self.depth):
        for d in range(1):
            qc.ry(theta[0],qb[0])
            qc.rz(theta[1],qb[0])
            # Entanglers
            qc.cx(qb[0],qb[1])
            qc.x(qb[0])
            qc.cx(qb[0],qb[2])
            qc.cx(qb[0],qb[3])
            qc.x(qb[0])
        return qc
    #def __call__(self,theta,qc,qb,qa=None,i=None):
    #    theta = self.sort_theta(theta)
    #    qc = self.prepare_Hartree_state(qc,qb)
    #    for d in range(self.depth):
    #        for i in range(self.l):
    #            qc.ry(theta[d,i,0],qb[i])
    #            qc.rz(theta[d,i,1],qb[i])
    #        for i in range(0,self.l,2):
    #            qc.cx(qb[i],qb[i+1])
    #        if not self.pair:
    #            for i in range(0,self.l,2):
    #                qc.cx(qb[i],qb[i+1])
    #    return qc

    def new_parameters(self):
        """
        New initial parameters.
        """
        #return 2*np.pi*np.random.randn(2*self.l*self.depth)
        return 2*np.pi*np.random.randn(2)

    def prepare_Hartree_state(self,qc,qb):
        """
        Prepares Hartree-Fock state. Assume initialized to ❘0...0⟩ before 
        Example for n=2 and l=4 -> ❘1100⟩
        """
        # Reference state
        if self.occupied: #❘1...0...⟩
            for k in range(self.n):
                qc.x(qb[k])
        else: #❘0...1...⟩
            for k in range(self.n,self.l):
                qc.x(qb[k])
        return qc

    def sort_theta(self,theta):
        return theta 
        #return theta.reshape(self.depth,self.l,2)


