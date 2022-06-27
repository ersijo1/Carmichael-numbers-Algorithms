"""
This program calculates Carmichael numbers following the algorithm created for my degree project.
Details are written in "Comparing characterizations of Carmichael numbers for computation" (Simon Johansson (2022)).
"""

import time
import math
import numpy as np

def loadPrimes(lim):
    #Loads primes from a .txt-file
    #param: lim(float): the sharp estimate (Theorem 10)
    #returns a list of the primes
    primes = []
    with open("primes.txt", "r") as file:
        for row in file:
            if int(row) < lim:
                primes.append(int(row))
    return primes

def compute(primes):
    """
    The function computing the Carmichael numbers (Algorithm 4)
    param: primes(int[]): list of prime numbers
    returns a list of Carmichael numbers
    """
    c = []
    for j in range(3, len(primes)):
        p = primes[j]
        upper_limit = min([int(((X/p)-1)/(p-1)), 1+ int((p**(d-1))/(p-1))])
        for f in range(2, upper_limit + 1):
            most_factors = []
            #From Theorem 7
            P = (f * (p-1)) +1
            #Factorization
            temp = P
            for i in range(1, primes.index(p)):
                if temp == 1:
                    break
                if temp%primes[i]==0:
                    if (P*p-1)%(primes[i]-1) == 0:
                        temp /= primes[i]
                        most_factors.append(primes[i])
                    else:
                        break
            if P/np.prod(most_factors) == 1:
                if len(most_factors) > 1:
                    c.append(p*np.prod(most_factors, dtype="int64"))
    return c

######################################################################################################################################

#START UP

#set d and X to the desired range
a_t = 1/(math.sqrt(2-(1/17)))
X = 10**8
lim = a_t * math.sqrt(X)
d = 5
primes = loadPrimes(lim)
start = time.process_time()
c = compute(primes)
end = time.process_time()
c.sort()
print(c)
print("No of numbers:", len(c))
print("Time:", end - start)

#controls the obtained numbers by comparing them to a list. List obtained from: https://www.numbersaplenty.com/set/Carmichael_number/more.php
a = """561, 1105, 1729, 2465, 2821, 6601, 8911, 10585, 15841, 29341, 41041, 46657, 52633, 62745, 63973, 75361, 101101, 115921, 126217, 162401, 172081, 188461, 252601, 278545, 294409, 314821, 334153, 340561, 399001, 410041, 449065, 488881, 512461, 530881, 552721, 656601, 658801, 670033, 748657, 825265, 838201, 852841, 997633, 1024651, 1033669, 1050985, 1082809, 1152271, 1193221, 1461241, 1569457, 1615681, 1773289, 1857241, 1909001, 2100901, 2113921, 2433601, 2455921, 2508013, 2531845, 2628073, 2704801, 3057601, 3146221, 3224065, 3581761, 3664585, 3828001, 4335241, 4463641, 4767841, 4903921, 4909177, 5031181, 5049001, 5148001, 5310721, 5444489, 5481451, 5632705, 5968873, 6049681, 6054985, 6189121, 6313681, 6733693, 6840001, 6868261, 7207201, 7519441, 7995169, 8134561, 8341201, 8355841, 8719309, 8719921, 8830801, 8927101, 9439201, 9494101, 9582145, 9585541, 9613297, 9890881, 10024561, 10267951, 10402561, 10606681, 10837321, 10877581, 11119105, 11205601, 11921001, 11972017, 12261061, 12262321, 12490201, 12945745, 13187665, 13696033, 13992265, 14469841, 14676481, 14913991, 15247621, 15403285, 15829633, 15888313, 16046641, 16778881, 17098369, 17236801, 17316001, 17586361, 17812081, 18162001, 18307381, 18900973, 19384289, 19683001, 20964961, 21584305, 22665505, 23382529, 25603201, 26280073, 26474581, 26719701, 26921089, 26932081, 27062101, 27336673, 27402481, 28787185, 29020321, 29111881, 31146661, 31405501, 31692805, 32914441, 33302401, 33596641, 34196401, 34657141, 34901461, 35571601, 35703361, 36121345, 36765901, 37167361, 37280881, 37354465, 37964809, 38151361, 38624041, 38637361, 39353665, 40160737, 40280065, 40430401, 40622401, 40917241, 41298985, 41341321, 41471521, 42490801, 43286881, 43331401, 43584481, 43620409, 44238481, 45318561, 45877861, 45890209, 46483633, 47006785, 48321001, 48628801, 49333201, 50201089, 53245921, 53711113, 54767881, 55462177, 56052361, 58489201, 60112885, 60957361, 62756641, 64377991, 64774081, 65037817, 65241793, 67371265, 67653433, 67902031, 67994641, 68154001, 69331969, 70561921, 72108421, 72286501, 74165065, 75151441, 75681541, 75765313, 76595761, 77826001, 78091201, 78120001, 79411201, 79624621, 80282161, 80927821, 81638401, 81926461, 82929001, 83099521, 83966401, 84311569, 84350561, 84417985, 87318001, 88689601, 90698401, 92625121, 93030145, 93614521, 93869665, 94536001, 96895441, 99036001, 99830641, 99861985"""
b = a.split(", ")
q = []
for num in b:
    q.append(int(num))
    if int(num) not in c and int(num) < X:
        print(num)
for num in c:
   if num not in q:
        print(num)

