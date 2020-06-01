import numpy as np
import qiskit as qk

class RYpairing:
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
        #if theta == None:
        #    theta = self.new_parameters()
        theta = self.sort_theta(theta)
        qc = self.prepare_Hartree_state(qc,qb)
        #for d in range(self.depth):
        for d in range(1):
            qc.u3(theta[0],0,0,qb[0]) # Ry gate
            # Entanglers
            qc.u3(np.pi,0,np.pi,qb[0])# X gate
            qc.cx(qb[0],qb[1])        # CNOT gate
            qc.cx(qb[0],qb[2])        # CNOT gate
            qc.cx(qb[0],qb[3])        # CNOT gate
            qc.u3(np.pi,0,np.pi,qb[0])# X gate
        return qc

    def new_parameters(self):
        """
        New initial parameters.
        """
        return 2*np.pi*np.random.randn(1)

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

def AND(qc,qb,ctrl,targ):
    c1 = qb[ctrl[0]]
    c2 = qb[ctrl[1]]
    t = qb[targ]
    # CH gate
    qc.u3(-np.pi/4,0,0,t) # Ry gate
    qc.cx(c1,t)
    qc.u3(np.pi/4,0,0,t) # Ry gate
    # CoX
    qc.u3(np.pi,0,np.pi,c2)# X gate
    qc.cx(c2,t)
    qc.u3(np.pi,0,np.pi,c2)# X gate
    # CH gate
    qc.u3(-np.pi/4,0,0,t) # Ry gate
    qc.cx(c1,t)
    qc.u3(np.pi/4,0,0,t) # Ry gate
    return qc    

class RYRZ:
    """
    """

    def __init__(self,n,l,occupied=1,pair=True,depth=1):
        self.n = n
        self.l = l
        self.occupied = occupied
        self.depth = depth

    def __call__(self,theta,qc,qb,qa=None,i=None):
        """
        Pairing ansatz. For 2 particles and 4 orbitals.
        """
        #if theta == None:
        #    theta = self.new_parameters()
        theta = self.sort_theta(theta)
        qc = self.prepare_Hartree_state(qc,qb)
        ### Rotation layer
        for j in range(self.depth):
            for i in range(self.n):
                qc.u3(theta[i,2*j],0,0,qb[i]) # Ry gate
                qc.u1(theta[i,2*j+1],qb[i])     # Rz gate
                if i == 0:
                    # For coupling map
                    #qc.swap(qb[1],qb[2])
                    #qc.cx(qb[0],qb[2])
                    #qc.swap(qb[1],qb[2])
                    qc.cx(qb[i],qb[(i+1)%2])
        ### Singles excitations
        for i in range(self.n,self.l,2):
            for j in range(self.depth):
                qc.u1(-0.5*theta[i,-2*j-1],qb[i])     # Rz gate
                qc.u3(-0.5*theta[i,-2*j-2],0,0,qb[i]) # Ry gate
        for i in range(self.n):
            qc.cx(qb[i],qb[self.n])
        for i in range(self.n,self.l,2):
            for j in range(self.depth):
                qc.u3(0.5*theta[i,2*j],0,0,qb[i]) # Ry gate
                qc.u1(0.5*theta[i,2*j+1],qb[i])     # Rz gate
        for i in range(self.n):
            qc.cx(qb[i],qb[self.n+1])
        ### Doubles excitation
        qc = AND(qc,qb,[0,1],2)
        qc.cx(qb[self.n],qb[self.n+1])# CNOT gate
        return qc

    def new_parameters(self):
        """
        New initial parameters.
        """
        return np.zeros((self.n+self.l)*self.depth)
        #return 2*np.pi*np.random.randn((self.n+self.l)*self.depth)

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
        assert len(theta)%2 == 0
        return theta.reshape((int(len(theta)/(2*self.depth)),2*self.depth)) 


