class Operator:
    def __init__(self,operation):
        self.op = operation.lower()
        if operation == 'I':
            self.op = ''
        self.factor = 1
        self.im = False
        self.ladder = ['+','-']
        self.pauli= ['x','y','z']
        if self.op in self.pauli:
            self.ind = self.pauli.index(self.op)
        if self.op in self.ladder:
            self.ind = self.ladder.index(self.op)

    def __mul__(self,other):
        getSingle = self.getSingle
        getDouble = self.getDouble
        mulcopy = self.mulcopy
        # PAULI
        if self.op in self.pauli:
            if other.op == self.op:
                i,factor = mulcopy(other)
                self.op = ''
                self.factor = factor
                self.im = i
                return self
            elif other.op != self.op and other.op in self.pauli:
                op_ind = 3 - (self.ind + other.ind)
                i,factor = mulcopy(other)
                self.op = self.pauli[op_ind]
                self.check_ind()
                self.im = i
                self.factor = factor
                self.im = not self.im
                if 3*self.ind + other.ind not in [1,5,6]:
                    self.factor *= -1
                return self
            # LADDER
            elif other.op in self.ladder:
                if self.op == 'x':
                    z_fact = eval('{}1'.format(other.op))
                    return getDouble(other,'I','z',0.5,z_fact)
                if self.op == 'y':
                    z_fact = eval('{}1'.format(other.op))
                    return getDouble(other,'I','z',0.5,z_fact,imag=True) 
                if self.op == 'z':
                    fact = eval('{}1'.format(other.ladder[(other.ind+1)%2]))
                    return getSingle(other,other.op,fact)
            else: # If identity
                self.op += other.op
                self.check_ind()
                i,factor = mulcopy(other)
                self.im = i
                self.factor = factor
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
                    z_fact = eval('{}1'.format(other.op))
                    return getDouble(other,'I','z',0.5,z_fact)
            # PAULI
            elif other.op in self.pauli:
                if other.op == 'x':
                    z_fact = eval('{}1'.format(other.ladder[(self.ind+1)%2]))
                    return getDouble(other,'I','z',0.5,z_fact)
                elif other.op == 'y':
                    z_fact = eval('{}1'.format(other.ladder[(self.ind+1)%2]))
                    return getDouble(other,'I','z',0.5,z_fact,imag=True) 
                elif other.op == 'z':
                    fact = eval('{}1'.format(self.op))
                    return getSingle(other,self.op,fact)
            else: #If identity
                self.op += other.op
                self.check_ind()
                i,factor = mulcopy(other)
                self.im = i
                self.factor = factor
                return self

        # OTHER/IDENTITY
        else:
            self.op += other.op
            self.check_ind()
            i,factor = mulcopy(other)
            self.im = i
            self.factor = factor
            return self

    def __str__(self):
        sign_fact = ''
        if self.factor < 0:
            if self.factor == -1:
                sign_fact = '-'
            else:
                sign_fact = str(self.factor)
        else:
            if self.factor != 1:
                sign_fact = str(self.factor)

        return '{}{}{}{}'.format(sign_fact,
                                 'i' if self.im else '',
                                 '*' if sign_fact != '' else '',
                                 self.op if self.op != '' else 'I')
    
    def __invert__(self): # Assume unitary
        # LADDER
        if self.op in self.ladder:
            self.ind = (self.ind + 1)%2
            self.op = self.ladder[self.ind]
            if self.im:
                self.factor *= -1
            return self.op
        # OTHER/ IDENTITY PAULI
        else:
            if self.im:
                self.factor *= -1
            return self

    def check_ind(self):
        if self.op in self.pauli:
            self.ind = self.pauli.index(self.op)
        if self.op in self.ladder:
            self.ind = self.ladder.index(self.op)

    def getSingle(self,other,new,sign=1):
        """
        Handle factor and imaginary part of multiplication 
        of two operators when single operator is the product.
        """
        op1 = Operator(new)
        i,factor = self.mulcopy(other)
        op1.factor = factor*sign
        op1.im = i
        return op1

    def getDouble(self,other,new1,new2,factor,sign,imag=False):
        """
        When multiplication of two operators lead to summation of two new operators
        new1,new2 -> str of new operators
        factor    -> Factor of sum
        sign      -> Sign between operators
        """
        op1 = Operator(new1)
        op2 = Operator(new2)
        i,new_factor = self.mulcopy(other)
        op1.factor *= factor*new_factor
        op2.factor *= factor*new_factor*sign
        if imag:
            if i:
                op1.factor *= -1
                op2.factor *= -1
            else:
                op1.im = True
                op2.im = True
        else:
            op1.im = i
            op2.im = i
        print(op1.factor)
        return [op1,op2]

    def mulcopy(self,other):
        """
        When multiplying two operators. Will handle factor and imaginary part.
        """
        i1 = self.im
        i2 = other.im
        i = False
        factor = self.factor*other.factor
        if i1 and not i2 or i2 and not i1:
            i = True
        if i1 and i2:
            i = False
            factor *= -1
        if not i1 and not i2:
            i = False
        return i,factor
    
    def ladder2pauli(self):
        if self.op not in self.ladder:
            print('Operator is not ladder')
            return self
        else:
            op1 = Operator('x')
            op2 = Operator('y')
            op1.im = self.im
            if self.im:
                op2.factor *= -1
                op2.im = False
            else:
                op2.im = True
            op1.factor *= self.factor
            op2.factor *= self.factor
            if self.op == '-':
                op2.factor *= -1
            return [op1,op2]

