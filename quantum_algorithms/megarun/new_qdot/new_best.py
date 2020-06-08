import numpy as np

def get_best(a,ac,b,bc,K):
    assert len(a) == len(b)
    assert len(a) == len(K)
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

f = np.load('fci.npy')
b = np.load('omegas.npy')


#ucc1 = np.load('E_ideal_UCCSDr_old.npy')
#coef1 = np.load('coeff_ideal_UCCSDr.npy')
#
#ucc2 = np.load('E_ideal_UCCSDr2.npy')
#coef2 = np.load('coeff_ideal_UCCSDr2.npy')
#
#newE,newC = get_best(ucc1,coef1,ucc2,coef2,f)
#
uccE = np.load('UCCSD_E_i.npy')
uccC = np.load('UCCSD_coeff_i.npy')
#
#
#ryrz1 = np.load('E_ideal_RYRZ.npy')
#coef1 = np.load('coeff_ideal_RYRZ.npy')
#
#ryrz2 = np.load('E_ideal_RYRZ2.npy')
#coef2 = np.load('coeff_ideal_RYRZ2.npy')
#

ryE = np.load('RYRZ_E_i.npy')
ryC = np.load('RYRZ_coeff_i.npy')

Un = np.load('new_UCCSD_E_i.npy')
UnC = np.load('new_UCCSD_coeff_i.npy')
Rn = np.load('new_RYRZ_E_i.npy')
RnC = np.load('new_RYRZ_coeff_i.npy')

newE,newC = get_best(ryE,ryC,Rn,RnC,f)
for i in range(len(f)):
    print('{:<.5f} -> {:<.5f}, {:<.5f} <= {:<.5f}'.format(f[i],ryE[i],Rn[i],newE[i]))


nuE,nuC = get_best(uccE,uccC,Un,UnC,f)
for i in range(len(f)):
    print('{:<.5f} -> {:<.5f}, {:<.5f} <= {:<.5f}'.format(f[i],uccE[i],Un[i],nuE[i]))

np.save('UCCSD_best_E.npy',nuE)
np.save('UCCSD_best_coeff.npy',nuC)

np.save('RYRZ_best_E.npy',newE)
np.save('RYRZ_best_coeff.npy',newC)
