import numpy as np
import qiskit as qk

def UCCS(system,qc,qb,theta,depth=10):
    n_qubits = system[0]
    n_fermi = system[1]
    theta = theta.reshape(n_qubits,n_qubits-n_fermi)

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
                
                # Pauli strings
                qc.h(qb[i])
                qc.rx(np.pi/2,qb[i])
                qc.h(qb[a])

                qc.cx(qb[i],qb[a])
                qc.rz(theta[i,a-n_fermi]/2,qb[a])
                qc.cx(qb[i],qb[a])

                qc.rx(-np.pi/2,qb[i])
                qc.rx(np.pi/2,qb[a])

                qc.cx(qb[i],qb[a])
                qc.rz(-theta[i,a-n_fermi]/2,qb[a])
                qc.cx(qb[i],qb[a])

                qc.h(qb[i])
                qc.rx(-np.pi/2,qb[a])
                qc.h(qb[a])

    return qc

def UCCD(system,qc,qb,theta,depth=1):
    n_qubits = system[0]
    n_fermi = system[1]
    theta = theta.reshape(n_fermi,n_fermi,n_qubits-n_fermi,n_qubits-n_fermi)

    # Prepare reference state
    for k in range(n_fermi,n_qubits):
        qc.x(qb[k])

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

                        # Pauli strings
                        qc.h(qb[i])
                        qc.h(qb[j])
                        qc.h(qb[a])
                        qc.rx(np.pi/2,qb[a])
                        qc.h(qb[b])

                        qc.cx(qb[i],qb[b])
                        qc.cx(qb[j],qb[b])
                        qc.cx(qb[a],qb[b])

                        qc.rz(theta[i,j,a-n_fermi,b-n_fermi]/8,qb[b]) #(xxyx)
                        
                        qc.cx(qb[a],qb[b])
                        qc.cx(qb[j],qb[b])

                        qc.rx(np.pi/2,qb[j])
                        qc.rx(-np.pi/2,qb[a])

                        qc.cx(qb[j],qb[b])
                        qc.cx(qb[a],qb[b])

                        qc.rz(-theta[i,j,a-n_fermi,b-n_fermi]/8,qb[b]) #(xyxx)
                        
                        qc.cx(qb[j],qb[b])
                        qc.cx(qb[i],qb[b])

                        qc.rx(np.pi/2,qb[i])
                        qc.rx(-np.pi/2,qb[j])
                        
                        qc.cx(qb[i],qb[b])
                        qc.cx(qb[j],qb[b])

                        qc.rz(-theta[i,j,a-n_fermi,b-n_fermi]/8,qb[b]) #(yxxx)

                        qc.cx(qb[a],qb[b])
                        qc.cx(qb[j],qb[b])

                        qc.rx(np.pi/2,qb[j])
                        qc.rx(np.pi/2,qb[a])

                        qc.cx(qb[j],qb[b])
                        qc.cx(qb[a],qb[b])

                        qc.rz(-theta[i,j,a-n_fermi,b-n_fermi]/8,qb[b]) #(yyyx)
                        
                        qc.cx(qb[a],qb[b])
                        qc.cx(qb[j],qb[b])
                        qc.cx(qb[i],qb[b])

                        qc.rx(-np.pi/2,qb[i])
                        qc.rx(np.pi/2,qb[b])

                        qc.cx(qb[i],qb[b])
                        qc.cx(qb[j],qb[b])
                        qc.cx(qb[a],qb[b])

                        qc.rz(theta[i,j,a-n_fermi,b-n_fermi]/8,qb[b]) #(xyyy)

                        qc.cx(qb[j],qb[b])
                        qc.cx(qb[i],qb[b])

                        qc.rx(np.pi/2,qb[i])
                        qc.rx(-np.pi/2,qb[j])
                        
                        qc.cx(qb[i],qb[b])
                        qc.cx(qb[j],qb[b])

                        qc.rz(theta[i,j,a-n_fermi,b-n_fermi]/8,qb[b]) #(yxyy)

                        qc.cx(qb[a],qb[b])
                        qc.cx(qb[j],qb[b])

                        qc.rx(np.pi/2,qb[j])
                        qc.rx(-np.pi/2,qb[a])

                        qc.cx(qb[j],qb[b])
                        qc.cx(qb[a],qb[b])

                        qc.rz(-theta[i,j,a-n_fermi,b-n_fermi]/8,qb[b]) #(yyxy)

                        qc.cx(qb[j],qb[b])
                        qc.cx(qb[i],qb[b])

                        qc.rx(-np.pi/2,qb[i])
                        qc.rx(-np.pi/2,qb[j])
                        
                        qc.cx(qb[i],qb[b])
                        qc.cx(qb[j],qb[b])

                        qc.rz(theta[i,j,a-n_fermi,b-n_fermi]/8,qb[b]) #(xxxy)

                        qc.cx(qb[a],qb[b])
                        qc.cx(qb[j],qb[b])
                        qc.cx(qb[i],qb[b])

                        qc.h(qb[i])
                        qc.h(qb[j])
                        qc.h(qb[a])
                        qc.rx(-np.pi/2,qb[b])
                        qc.h(qb[b])

    return qc

def sort_theta(theta):
    """
    Return array with two elements corresponding to singles and doubles 
    parameters respectively.
    """
    return theta

def UCC(system,qc,qb,theta,terms='SD',depth=1):
    if terms == 'S':
        qc = UCCS(system,qc,qb,theta,depth=depth)
    elif terms == 'D':
        qc = UCCD(system,qc,qb,theta,depth=depth)
    elif terms == 'SD':
        theta = sort_theta(theta)
        qc = UCCS(system,qc,qb,theta[0],depth=depth)
        qc = UCCD(system,qc,qb,theta[1],depth=depth)
    else:
        raise ValueError('Need to spesify correct UCC type')
    return qc
