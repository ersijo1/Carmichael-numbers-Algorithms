"""
This program calculates Carmichael numbers following the algorithm given by Richard G. E. Pinch in "The Carmichael Numbers up to 10^15" (1993).
It was created in conjunction with my degree project: "Comparing characterizations of Carmichael numbers for computation" (Simon Johansson (2022)).
"""
import time
import math
import numpy as np

#USEFUL FUNCTIONS

def loadPrimes(lim):
    #Loads primes from a .txt-file
    #param: lim(float): the sharp estimate (Theorem 10)
    #returns a list of the primes
    primes = []
    with open("primes.txt", "r") as file:
        for row in file:
            if int(row) <= lim:
                primes.append(int(row))
    return primes

def check_primality(p, primes):
    #checks if the number is part of the prime list
    #param: p(int): the number to be checked, primes(int[]): the list of primes
    #returns True if prime
    if p in primes:
        return True
    else:
        return False

def korselt(factors):
    #Checks if the factors satisfies the divisibility of Korselt's criterion
    #param: factors(int[]): the factors to be checked
    #returns True if all satisfy the criterion
    N = np.prod(factors)
    for factor in factors:
        if ((N-1) % (factor-1) != 0):
            return False
    return True

def lcm(factors):
    #Calculates the lcm of the given numbers
    #param: factors(int[]): the numbers to examine
    #returns the lcm of the numbers
    return np.lcm.reduce(factors)

def cong_modulo(num1, num2, divisor):
    #Checks if num1 is congruent num2 (mod divisor)
    #returns True if num1 is congruent
    if (num1 - num2) % divisor == 0:
        return True
    else:
        return False


def get_factors(primes, X, d):
    #obtains successive lists of primes
    #param: primes(int[]): the list of primes
    #X(int): the upper limit
    #d(int): the number of factors
    #returns a list of possible prime combinations
    most_factors = []
    for i in range(1,len(primes)):
        most_factors.append([primes[i]])
    starting_point=0
    for a in range(4,d+1):
        end_point = len(most_factors)
        for j in range(starting_point, end_point):
            for k in range(primes.index(most_factors[j][-1])+1, len(primes)):
                if primes[k] < (X/np.prod(most_factors[j]))**(1/(a-len(most_factors[j]))):
                    most_factors.append(most_factors[j] + [primes[k]])
        starting_point = end_point
    return most_factors


############################################################################################################################

#CALCULATION

def low_P(primes, P, most_factors):
    """
    Calculates Carmichael numbers from a sequence of primes according to Proposition 12. (Algorithm 2)
    param: primes(int[]): list of primes
    P(int): the product of the sequence of primes
    most_factors(int[]): the sequence of possible prime factors up to d-2
    returns a list of the Carmichael numbers for that sequence of primes
    """
    p_d2 = most_factors[-1]
    c_num = []
    for D in range(2, P):
        #Proposition 12(3)
        lower_limit = int((P**2)/D)            
        upper_limit = int((P**2)*(p_d2+3)/(D*(p_d2+1)))
        for C in range(lower_limit+1, upper_limit+1):
            delta = C*D-(P**2)
            #Proposition 12(1)
            p_d1 = ((P-1)*(P+D)/(delta))+1
            ##Proposition 12(2)
            p_d = ((P-1)*(P+C)/(delta))+1
            #checks if the resulting factors make a Carmichael number
            isPrime = check_primality(p_d1, primes) and check_primality(p_d, primes)
            factors = most_factors + [p_d1, p_d]
            new_P = P*p_d1
            if (isPrime and p_d < (X/new_P) and p_d1> most_factors[-1]):
                if korselt(factors):
                    c_num.append(int(np.prod(most_factors) * p_d1 * p_d))
    return c_num

def high_P(primes, X, P, most_factors):
    """
    Calculates Carmichael numbers from a sequence of primes according to Proposition 11(ii). (Algorithm 3)
    param: primes(int[]): list of primes
    X(int): the upper limit
    P(int): the product of the sequence of primes
    most_factors(int[]): the sequence of possible prime factors up to d-2
    returns a list of the Carmichael numbers for that sequence of primes
    """
    cNum = []
    #Proposition 11(i) and Proposition 13
    upper_limit = min([2*(P**2), (X/P)**(1/2)])
    for i in range(primes.index(most_factors[-1]), len(primes)):
        if primes[i] > upper_limit:
            break
        if primes[i] not in most_factors:
            new_P = P*primes[i]
            factors = most_factors + [primes[i]]
            L = lcm(np.subtract(factors, 1))
            #The algorithm is here split in two, "reach" is the index of the prime at the point of the split. 
            reach = len(primes)//6
            alg = high_P_low(i, primes, new_P, L, reach, X)
            if (new_P*primes[reach] < X):
                alg += high_P_high(primes, new_P, L, reach, X)
            if alg != []:
                for p_d in alg:
                    cNum.append(new_P * p_d)
    return cNum

def high_P_low(i, primes, new_P, L, reach, X):
    """
    Calculates for lower values of the largest factor. I found skipping the calculation of P' to be more efficient.
    param: i(int): index of the previous prime factor
    primes(int[]): list of primes
    new_P(int): the product of the prime factors, except the last one
    L(int): the lcm according to Proposition 11(ii)
    reach(int): the largest index checked in this function
    X(int): the upper limit
    returns a list of viable largest factors
    """
    p_d = []
    for j in range(i+1, reach):
        #makes sure too large Carmichael numbers aren't calculated
        if primes[j]*new_P > X:
            break
        if cong_modulo(new_P*primes[j], 1, L):
            if (new_P-1)%(primes[j]-1)==0:
                if primes[j]*new_P > X:
                    print(primes[j])
                p_d.append(primes[j])
    return p_d


def high_P_high(primes, new_P, L, reach, X):
    """
    Calculates for lower values of the largest factor. I found skipping the calculation of P' to be more efficient.
    param: i(int): primes(int[]): list of primes
    new_P(int): the product of the prime factors, except the last one
    L(int): the lcm according to Proposition 11(ii)
    reach(int): the smallest index checked in this function
    X(int): the upper limit
    returns a list of viable largest factors
    """
    p_d = []
    #makes sure too large Carmichael numbers aren't calculated
    start = max([round(new_P*(new_P-1)/(X-1)),2])
    end = (new_P-1)//(primes[reach]-1)
    for j in range(start, end+1):
        if (new_P-1) % j == 0:
            p = int(((new_P-1)/j)+1)
            if new_P*p < X and cong_modulo(p*new_P, 1, L) and p in primes:
                p_d.append(p)
    return p_d


def compute(primes, X, d):
    """
    main function for the computation. Obtains lists of primes and then splits into the two different methods depending on the size of P. (Algorithm 1)
    param: primes(int[]): list of primes
    X(int): upper limit
    d(int): number of factors
    returns a list of Carmichael numbers
    """
    carmichael = []
    most_factors = []
    #obtains lists of primes with a maximum size of d-2.
    most_factors += get_factors(primes,X,d)
    for i in range(len(most_factors)):
        P = np.prod(most_factors[i], dtype="int64")
        #splits into two algorithms. I found P=53 to be the most efficient place to split.
        if P < 53:
            #for lower values of P
            c_num = low_P(primes, P, most_factors[i])
            if c_num != []:
                for c in c_num:
                    carmichael.append(c)

        else:
            #for higher values
            c_num = high_P(primes, X, P, most_factors[i])
            if c_num != []:
                for c in c_num:
                    carmichael.append(c)
    return carmichael

##############################################################################################################################

#START UP

#set d and X to the desired range
d = 5
X = 10**7
a_t = 1/(math.sqrt(2-(1/17)))
lim = a_t * math.sqrt(X)
primes = loadPrimes(lim)
start = time.process_time()
carmichael = compute(primes, X, d)
end = time.process_time()
carmichael.sort()
print(carmichael)
print(end - start)
print("No of numbers:", len(carmichael))

#controls the obtained numbers by comparing them to a list. List obtained from: https://www.numbersaplenty.com/set/Carmichael_number/more.php
a = """561, 1105, 1729, 2465, 2821, 6601, 8911, 10585, 15841, 29341, 41041, 46657, 52633, 62745, 63973, 75361, 101101, 115921, 126217, 162401, 172081, 188461, 252601, 278545, 294409, 314821, 334153, 340561, 399001, 410041, 449065, 488881, 512461, 530881, 552721, 656601, 658801, 670033, 748657, 825265, 838201, 852841, 997633, 1024651, 1033669, 1050985, 1082809, 1152271, 1193221, 1461241, 1569457, 1615681, 1773289, 1857241, 1909001, 2100901, 2113921, 2433601, 2455921, 2508013, 2531845, 2628073, 2704801, 3057601, 3146221, 3224065, 3581761, 3664585, 3828001, 4335241, 4463641, 4767841, 4903921, 4909177, 5031181, 5049001, 5148001, 5310721, 5444489, 5481451, 5632705, 5968873, 6049681, 6054985, 6189121, 6313681, 6733693, 6840001, 6868261, 7207201, 7519441, 7995169, 8134561, 8341201, 8355841, 8719309, 8719921, 8830801, 8927101, 9439201, 9494101, 9582145, 9585541, 9613297, 9890881, 10024561, 10267951, 10402561, 10606681, 10837321, 10877581, 11119105, 11205601, 11921001, 11972017, 12261061, 12262321, 12490201, 12945745, 13187665, 13696033, 13992265, 14469841, 14676481, 14913991, 15247621, 15403285, 15829633, 15888313, 16046641, 16778881, 17098369, 17236801, 17316001, 17586361, 17812081, 18162001, 18307381, 18900973, 19384289, 19683001, 20964961, 21584305, 22665505, 23382529, 25603201, 26280073, 26474581, 26719701, 26921089, 26932081, 27062101, 27336673, 27402481, 28787185, 29020321, 29111881, 31146661, 31405501, 31692805, 32914441, 33302401, 33596641, 34196401, 34657141, 34901461, 35571601, 35703361, 36121345, 36765901, 37167361, 37280881, 37354465, 37964809, 38151361, 38624041, 38637361, 39353665, 40160737, 40280065, 40430401, 40622401, 40917241, 41298985, 41341321, 41471521, 42490801, 43286881, 43331401, 43584481, 43620409, 44238481, 45318561, 45877861, 45890209, 46483633, 47006785, 48321001, 48628801, 49333201, 50201089, 53245921, 53711113, 54767881, 55462177, 56052361, 58489201, 60112885, 60957361, 62756641, 64377991, 64774081, 65037817, 65241793, 67371265, 67653433, 67902031, 67994641, 68154001, 69331969, 70561921, 72108421, 72286501, 74165065, 75151441, 75681541, 75765313, 76595761, 77826001, 78091201, 78120001, 79411201, 79624621, 80282161, 80927821, 81638401, 81926461, 82929001, 83099521, 83966401, 84311569, 84350561, 84417985, 87318001, 88689601, 90698401, 92625121, 93030145, 93614521, 93869665, 94536001, 96895441, 99036001, 99830641, 99861985"""
b = a.split(", ")
q = []
for num in b:
    q.append(int(num))
    if int(num) not in carmichael and int(num) < X:
        print(num)
for num in carmichael:
   if num not in q:
        print(num)