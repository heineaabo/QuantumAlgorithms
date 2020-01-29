import numpy as np
import qiskit as qk

class UnitaryCoupledCluster:
    """
    Unitary Coupled Cluster ansatz, with optimized singles and doubles circuits.
    Assumes the convention:
        ❘0⟩ - occupied orbital.
        ❘1⟩ - unoccupied orbital.

    Initialize class input parameters below.
    Initialize state to Hartree-Fock state with prepare_Hartree_state(circuit,qregister)
    Call class to perform UCC.

    Input:
            n (int): Number of particles.    
            l (int): Number of spin orbitals.
        trunc (str): CC truncation.
                     - 'S' for singles.
                     - 'D' for doubles.
                     - 'SD' for singles and doubles.
        depth (int): Depth of Trotter approximation. 
    """

    def __init__(self,n,l,trunc,depth):
        self.n = n
        self.l = l
        self.trunc = trunc
        self.depth = depth
        self.num_S, self.num_D = self.get_num_parameters()
        if 'S' not in self.trunc:
            self.num_S = 0
        if 'D' not in self.trunc:
            self.num_D = 0

    def __call__(self,qc,qb,theta):
        theta = self.sort_theta(theta)
        for rho in range(self.depth):
            if 'S' in self.trunc:
                qc = self.UCCS(qc,qb,theta[0],rho)
            if 'D' in self.trunc:
                qc = self.UCCD(qc,qb,theta[1],rho)
        if 'S' not in self.trunc and 'D' not in self.trunc:
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
        return S,D

    def new_parameters(self):
        return 2*np.pi*np.random.randn((self.num_S+self.num_D)*self.depth)


    def prepare_Hartree_state(self,qc,qb):
        """
        Prepares Hartree-Fock state. Assume initialized to ❘0...0⟩ before 
        Example for n=2 and l=4 -> ❘0011⟩
        """
        # Reference state
        for k in range(self.n):
            qc.x(qb[k])
        return qc

    def sort_theta(self,theta):
        """
        Sort parameters to two channels for singles and doubles.
        TODO: 
            - Sort for depth > 1.
        """
        singles = self.num_S
        doubles = self.num_D
        assert len(theta) == (singles+doubles)*self.depth
        #new_theta = [theta[:singles],theta[singles:]]
        thetaS = []
        thetaD = []
        for i in range(self.depth):
            thetaS.append(theta[(i*singles):(i+1)*singles])
            thetaD.append(theta[(self.depth*singles + (i)*doubles):(self.depth*singles + (i+1)*doubles)])
        #new_theta = [theta[:singles],theta[singles:]]
        new_theta = [thetaS,thetaD]
        return new_theta

    def UCCS(self,qc,qb,theta,rho):
        """
        Optimized UCC singles implementation.
        """
        n_qubits = self.l
        n_fermi = self.n
        #theta = theta.reshape(n_fermi,n_qubits-n_fermi)

        t = 0
        # UCCS
        for i in range(n_fermi):
            for a in range(n_fermi,n_qubits):
                # Pauli-z
                for k in range(i+1,a):
                    qc.z(qb[k])
                # Pauli strings
                qc.h(qb[i])
                qc.rx(np.pi/2,qb[i])
                qc.h(qb[a])

                qc.cx(qb[i],qb[a])

                #qc.rz(theta[i,a-n_fermi]/2,qb[a])
                qc.rz(theta[rho][t]/2,qb[a])

                qc.cx(qb[i],qb[a])

                qc.rx(-np.pi/2,qb[i])
                qc.rx(np.pi/2,qb[a])

                qc.cx(qb[i],qb[a])

                #qc.rz(-theta[i,a-n_fermi]/2,qb[a])
                qc.rz(-theta[rho][t]/2,qb[a])

                qc.cx(qb[i],qb[a])

                qc.h(qb[i])
                qc.rx(-np.pi/2,qb[a])
                qc.h(qb[a])

                t += 1

        return qc

    def UCCD(self,qc,qb,theta,rho=0):
        """
        Optimized UCC doubles implementation.
        """
        n_qubits = self.l
        n_fermi = self.n
        #theta = theta.reshape(n_fermi,n_fermi,n_qubits-n_fermi,n_qubits-n_fermi)

        t = 0
        # UCCD
        for i in range(n_fermi):
            for j in range(i+1,n_fermi):
                for a in range(n_fermi,n_qubits):
                    for b in range(a+1,n_qubits):
                        # Pauli-z
                        for k in range(i+1,j):
                            qc.z(qb[k])
                        for l in range(a+1,b):
                            qc.z(qb[l])

                        # Pauli strings
                        qc.h(qb[i])
                        qc.h(qb[j])
                        qc.h(qb[a])
                        qc.rx(np.pi/2,qb[a])
                        qc.h(qb[b])

                        qc.cx(qb[i],qb[b])
                        qc.cx(qb[j],qb[b])
                        qc.cx(qb[a],qb[b])

                        #qc.rz(theta[i,j,a-n_fermi,b-n_fermi]/8,qb[b]) #(xxyx)
                        qc.rz(theta[rho][t]/8,qb[b]) #(xxyx)
                        
                        qc.cx(qb[a],qb[b])
                        qc.cx(qb[j],qb[b])

                        qc.rx(np.pi/2,qb[j])
                        qc.rx(-np.pi/2,qb[a])

                        qc.cx(qb[j],qb[b])
                        qc.cx(qb[a],qb[b])

                        #qc.rz(-theta[i,j,a-n_fermi,b-n_fermi]/8,qb[b]) #(xyxx)
                        qc.rz(-theta[rho][t]/8,qb[b]) #(xyxx)
                        
                        qc.cx(qb[j],qb[b])
                        qc.cx(qb[i],qb[b])

                        qc.rx(np.pi/2,qb[i])
                        qc.rx(-np.pi/2,qb[j])
                        
                        qc.cx(qb[i],qb[b])
                        qc.cx(qb[j],qb[b])

                        #qc.rz(-theta[i,j,a-n_fermi,b-n_fermi]/8,qb[b]) #(yxxx)
                        qc.rz(-theta[rho][t]/8,qb[b]) #(yxxx)

                        qc.cx(qb[a],qb[b])
                        qc.cx(qb[j],qb[b])

                        qc.rx(np.pi/2,qb[j])
                        qc.rx(np.pi/2,qb[a])

                        qc.cx(qb[j],qb[b])
                        qc.cx(qb[a],qb[b])

                        #qc.rz(-theta[i,j,a-n_fermi,b-n_fermi]/8,qb[b]) #(yyyx)
                        qc.rz(-theta[rho][t]/8,qb[b]) #(yyyx)
                        
                        qc.cx(qb[a],qb[b])
                        qc.cx(qb[j],qb[b])
                        qc.cx(qb[i],qb[b])

                        qc.rx(-np.pi/2,qb[i])
                        qc.rx(np.pi/2,qb[b])

                        qc.cx(qb[i],qb[b])
                        qc.cx(qb[j],qb[b])
                        qc.cx(qb[a],qb[b])

                        #qc.rz(theta[i,j,a-n_fermi,b-n_fermi]/8,qb[b]) #(xyyy)
                        qc.rz(theta[rho][t]/8,qb[b]) #(xyyy)

                        qc.cx(qb[j],qb[b])
                        qc.cx(qb[i],qb[b])

                        qc.rx(np.pi/2,qb[i])
                        qc.rx(-np.pi/2,qb[j])
                        
                        qc.cx(qb[i],qb[b])
                        qc.cx(qb[j],qb[b])

                        #qc.rz(theta[i,j,a-n_fermi,b-n_fermi]/8,qb[b]) #(yxyy)
                        qc.rz(theta[rho][t]/8,qb[b]) #(yxyy)

                        qc.cx(qb[a],qb[b])
                        qc.cx(qb[j],qb[b])

                        qc.rx(np.pi/2,qb[j])
                        qc.rx(-np.pi/2,qb[a])

                        qc.cx(qb[j],qb[b])
                        qc.cx(qb[a],qb[b])

                        #qc.rz(-theta[i,j,a-n_fermi,b-n_fermi]/8,qb[b]) #(yyxy)
                        qc.rz(-theta[rho][t]/8,qb[b]) #(yyxy)

                        qc.cx(qb[j],qb[b])
                        qc.cx(qb[i],qb[b])

                        qc.rx(-np.pi/2,qb[i])
                        qc.rx(-np.pi/2,qb[j])
                        
                        qc.cx(qb[i],qb[b])
                        qc.cx(qb[j],qb[b])

                        #qc.rz(theta[i,j,a-n_fermi,b-n_fermi]/8,qb[b]) #(xxxy)
                        qc.rz(theta[rho][t]/8,qb[b]) #(xxxy)

                        qc.cx(qb[a],qb[b])
                        qc.cx(qb[j],qb[b])
                        qc.cx(qb[i],qb[b])

                        qc.h(qb[i])
                        qc.h(qb[j])
                        qc.h(qb[a])
                        qc.rx(-np.pi/2,qb[b])
                        qc.h(qb[b])
                        
                        t += 1

        return qc
