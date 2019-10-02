import numpy as np

class Operator:
    def __init__(self,operation=''):
        self.op = operation.upper()
        print(self.op)

    def __mul__(self,op2):
        return self.op + op2.op


class Pauli(Operator):
    def __inv__(self):
        return self.op

class Identity(Operator):
    def __mul__(self,op2):
        return op2.op
    
class Ladder(Operator):
    def __init__(self,operation):
        self.op = operation
        self.values = ['+','-']
        self.ind = self.values.index(self.op)
    
    def __invert__(self): 
        self.ind = (self.ind + 1)%2
        self.op = self.values[self.ind]
        return self.op


pauliX = Operator('x')
pauliY = Operator('y')
pauliZ = Pauli('Z')
I = Identity('')
plus = Ladder('+')
print(~Ladder)

print(pauliX * pauliY)
print(pauliX * pauliZ)
print(pauliX * I)
