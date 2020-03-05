import time

def levi_cevita(i,j):
    """
    Input two coordinate-axes i and j. k is found from that.
    Output -1 or 1.
    """
    k = 3-i-j
    epsilon = int(((i-j)*(j-k)*(k-i))/2)
    return epsilon

t1 = time.process_time()
print('x*y =',levi_cevita(0,1))
t2 = time.process_time()
print(t2-t1)
print('y*x =',levi_cevita(1,0))
print('x*z =',levi_cevita(0,2))
print('z*x =',levi_cevita(2,0))
print('y*z =',levi_cevita(1,2))
print('z*y =',levi_cevita(2,1))


t1 = time.process_time()
if 4 not in [1,2,3]:
    print('non')
t2 = time.process_time()
print(t2-t1)
