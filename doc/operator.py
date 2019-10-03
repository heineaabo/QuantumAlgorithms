class Operator:
    def __init__(self,operation=''):
        self.op = operation.upper()
        print(self.op)

    def __mul__(self,other):
        return self.op + other.op

    def __str__(self):
        return self.op

class Pauli(Operator):
    def __invert__(self):
        return self.op

class Identity(Operator):
    def __mul__(self,other):
        return other.op

    def __invert__(self):
        return self.op
    
class Ladder(Operator):
    def __init__(self,operation):
        self.op = operation
        self.values = ['+','-']
        self.ind = self.values.index(self.op)
    
    def __invert__(self): 
        self.ind = (self.ind + 1)%2
        self.op = self.values[self.ind]
        return self.op

