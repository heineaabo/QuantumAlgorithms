class Operator:
    def __init__(self):
        self.op = 'U'
        self.unit = 1
        self.im = False

    def __mul__(self,other):
        self.op += other.op
        return self

    def __str__(self):
        if self.op == '':
            return 'I'
        else:
            return self.op
    
    def __invert__(self): # Assume unitary
        return self
    
class Identity(Operator):
    def __init__(self):
        self.op = ''

class Pauli(Operator):
    def __init__(self,operation):
        self.op = operation.upper()
        self.axis = ['X','Y','Z']
        self.ind = self.axis.index(self.op)
        self.unit = 1
        self.im = False
        
    def __mul__(self,other):
        if other.op == self.op:
            self.op = ''
            return 'I'
        elif other.op != self.op and other.op in self.axis:
            op_ind = 3 - (self.ind + other.ind)
            self.op = self.axis[op_ind]
            self.im = not self.im
            if 3*self.ind + other.ind not in [1,5,6]:
                self.unit *= -1
            return self
        else:
            self.op += other.op
            return self
        
    def __str__(self):
        if self.op == '':
            return 'I'
        else:
            return '{}{}{}'.format('-' if self.unit==-1 else '',
                                   'i' if self.im else '',
                                   self.op)
        
class Ladder(Operator):
    def __init__(self,operation):
        self.op = operation
        self.values = ['+','-']
        self.ind = self.values.index(self.op)
        
    def __mul__(self,other):
        if other.op == self.op:
            self.op = '0'
        return self
    
    def __invert__(self): 
        self.ind = (self.ind + 1)%2
        self.op = self.values[self.ind]
        return self.op


