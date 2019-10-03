class Circuit:
    def __init__(self,n_qubits, cirq=None):
        self.n = n_qubits
        self.cirq = None
        if type(cirq) == type(list):
            self.cirq = cirq
            self.length = 1
        elif type(cirq) == type('str'):
            self.cirq == self.transform(cirq)
            self.length = 1

    def __str__(self):
        for q in range(self.n):
            for o in range(self.length):
                print('|q{}> '.format(q),end='')
                print('{:-^4s}{}{:-^4s}'.format('-',self.cirq[q] if self.cirq[q] != 'I' else '-','-')) 
                print('')
        return ''

    def add(self,cirq):
        self.cirq.append(cirq)

    def transform(self,cirq):
        l = ['' for i in range(self.n_qubits)]
        cirq = cirq.split('*')
        cirq = [[cirq[i][0],cirq[i][2]] for i in range(len(cirq))]
        for term in cirq:
            if term[0] == 'a':
                continue
            if term[b] == 'b':
                continue


   
CIRQ = 'b[a]*b[b]*a[j]*a[i]'
#CIRQ = ['S','I','I','Y']


C = Circuit(len(CIRQ),CIRQ)
#print(C)
