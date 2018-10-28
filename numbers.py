from random import randint 


n = 1
L = (1/12.0)




def L_Multiples(L,n):
    x = 0.0
    multi = []
    for i in range(1,n):
        x += L*(1.0/i)
        multi.append(x)
    return multi

def L_Difference(lmultiples):
    newMulti = []
    for i in xrange(len(lmultiples)):
        x = (lmultiples[i] - lmultiples[i-1])
        newMulti.append(x)
    return newMulti

def L_sum(lmultiples):
    newMulti = []
    for i in xrange(len(lmultiples)):
        x = (lmultiples[i] + lmultiples[i-1])
        newMulti.append(x)
    return newMulti


def sums(lSome):
    x = 0.0
    for i in xrange(0,len(lSome)):
        x += lSome[i]
    return x
'''
n = 1
while n == 1:
    into = int(input("choose the amount:"))
    if into == "stop":
        n = 0
    else:
    
    '''
        
zeroNumbers = []
def dod(n):
    afterOp = {}
    for i in range(n):
        a = sums(L_sum(L_Multiples(L, i)))
        if a == 0.0:
            print str(i) + " is zero"
            zeroNumbers.append(i)
        
        else:
            afterOp[i] = a
    return afterOp

print L_Multiples(L,12)
a = sums(L_sum(L_Multiples(L, 12)))
print a
b = sums(L_Difference(L_Multiples(L, 12)))
print b



