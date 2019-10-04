class Operator:
    def __init__(self,operation):
        self.op = operation.lower()
        if operation == 'I':
            self.op = ''
        self.unit = 1
        self.im = False
        self.ladder = ['+','-']
        self.pauli= ['x','y','z']
        if self.op in self.pauli:
            self.ind = self.pauli.index(self.op)
        if self.op in self.ladder:
            self.ind = self.ladder.index(self.op)

    def __mul__(self,other):
        # PAULI
        if self.op in self.pauli:
            if other.op == self.op:
                self.op = ''
                return 'I'
            elif other.op != self.op and other.op in self.pauli:
                op_ind = 3 - (self.ind + other.ind)
                self.op = self.pauli[op_ind]
                self.im = not self.im
                if 3*self.ind + other.ind not in [1,5,6]:
                    self.unit *= -1
                return self
            else:
                self.op += other.op
                return self
        # LADDER
        elif self.op in self.ladder:
            # LADDER
            if other.op in self.ladder:
                if self.ind == other.ind:
                    return 0
                else: # != -> 1/2*(I +- Z)
                    """
                    last operator sign determines + or - for Z
                    """
                    I = Operator(I)
                    Z = Operator('z')
                    I.factor *= 0.5
                    Z.factor *= eval('{}0.5'.format(other.op))
                    return [I,Z]
            # PAULI
            if other.op in self.pauli:
                if other.op == 'x':
                    I = Operator(I)
                    Z = Operator('z')
                    I.im = not 

                    I.factor *= 0.5
                    Z.factor *= eval('{}0.5'.format(other.ladder[(self.ind+1)%2]))
                    return [I,Z]
                if other.op == 'y':
                    I = Operator(I)
                    Z = Operator('z')
                    Z.unit *= eval('{}1'.format(other.ladder[(self.ind+1)%2]))
                    I.factor *= 0.5*
                    Z.factor *= 0.5
                    return [I,Z]

                if other.op == 'z':
                    self.unit *= eval('{}1'.format(other.op))



        # OTHER/IDENTITY
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
    
    def __invert__(self): # Assume unitary
        # LADDER
        if self.op in self.ladder:
            self.ind = (self.ind + 1)%2
            self.op = self.ladder[self.ind]
            return self.op
        # OTHER/ IDENTITY PAULI
        else:
            return self
