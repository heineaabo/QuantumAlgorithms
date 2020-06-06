import numpy as np

def get_best(a,ac,b,bc,K):
    best = np.zeros_like(a)
    coef = np.zeros_like(ac)
    i = 0
    for a1,b1 in zip(a,b):
        ar = np.abs(a1-K[i])
        br = np.abs(b1-K[i])
        if ar < br:
            best[i] = a1
            coef[i] = ac[i]
        else:
            best[i] = b1
            coef[i] = bc[i]
        i += 1
    return best,coef

f = np.load('../cc.npy')
f = f[[i for i in range(0,len(f),5)]]
b = np.load('../bonds.npy')
b = b[[i for i in range(0,len(b),5)]]
print(b)


#ucc1 = np.load('E_ideal_UCCSDr_old.npy')
#coef1 = np.load('coeff_ideal_UCCSDr.npy')
#
#ucc2 = np.load('E_ideal_UCCSDr2.npy')
#coef2 = np.load('coeff_ideal_UCCSDr2.npy')
#
#newE,newC = get_best(ucc1,coef1,ucc2,coef2,f)
#
uccE = np.load('UCCSDr_good_E.npy')
uccC = np.load('UCCSDr_good_coeff.npy')
#
#
#ryrz1 = np.load('E_ideal_RYRZ.npy')
#coef1 = np.load('coeff_ideal_RYRZ.npy')
#
#ryrz2 = np.load('E_ideal_RYRZ2.npy')
#coef2 = np.load('coeff_ideal_RYRZ2.npy')
#
#newE,newC = get_best(ryrz1,coef1,ryrz2,coef2,f)

ryE = np.load('RYRZ_good_E.npy')
ryC = np.load('RYRZ_good_coeff.npy')

Un = np.load('single_UCCSDr_E.npy')
UnC = np.load('single_UCCSDr_coef.npy')
Rn = np.load('single_RYRZ_E.npy')
RnC = np.load('single_RYRZ_coef.npy')

#print([f[8],f[16]])
#print('')
#print([uccE[8],uccE[16]])
#print(Un)
#print('')
#print([ryE[8],ryE[16]])
#print(Rn)
print(uccE)
print(ryE)

uccE[16] = Un[1]
uccC[16] = UnC[1]

ryE[8] = Rn[0]
ryC[8] = RnC[0]

print(uccE)
print(ryE)

np.save('UCCSDr_good_E.npy',uccE)
np.save('UCCSDr_good_coeff.npy',uccC)

np.save('RYRZ_good_E.npy',ryE)
np.save('RYRZ_good_coeff.npy',ryC)
