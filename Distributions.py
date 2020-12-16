### PROVIDE DESCRIPTIONS FOR DISTRIBUTIONS.

### TO-DO LIST:
# VET THE CODE.
# STANDARDIZE ERROR-PRINTOUT.
# STANDARDIZE DISTRIBUTION INFO PRINTOUT.
# ENSURE ALL DISTRIBUTIONS HAVE MEAN, VAR, SD, INF AND SUP.


import operator
import math
from math import *
from pyerf import erfinv
from random import random
from copy import deepcopy

# Global controls.
class DistEnv:
    differentiation_accuracy = 10e-4
    integration_slice_count = int(10e4)
    quantile_precision = 10e-4


# Error messages.
ERR_REAL_SAMPLE = "Your sample must be a list of real numbers."
ERR_SAMPLE_SIZE = "Sample size must be a positive integer."
ERR_RATE_PARAM = "Rate parameter must be a non-negative number."
ERR_LOC_PARAM = "Location parameter must be a numerical value."
ERR_SCALE_PARAM = "Scale parameter must be a positive number."
ERR_TIME = "Time must be a non-negative integer."
ERR_EMPTY_DATA = "No valid data detected."
ERR_DISTRIBUTION = "No parametric distribution exists for the result of this operation."
ERR_PROB = "Probability input must be a numerical value between zero and one exclusive."
ERR_VALUE = "Value supplied is beyond supported range."


# Functions for convenience.
Inf = math.inf

def ln(x):
    return log(x, e)

def beta(a, b):
    return gamma(a) * gamma(b) / gamma(a + b)

def ifelse(predicate, value_true, value_false):
    if predicate:
        return value_true
    else:
        return value_false



# Calculus-related functions.
def differentiate(function, x):
    dx = 1
    previous_slope = (function(x + dx) - function(x)) / dx
    while True:
        dx /= 10
        slope = (function(x + dx) - function(x)) / dx
        print("dx = " + str(dx) + ", slope = " + str(slope))
        if abs(slope - previous_slope) < DistEnv.differentiation_accuracy:
            return slope
        else:
            previous_slope = slope



def integrate(function, lower, upper):
    SLICE_COUNT = int(10e4)
    subtotal = 0
    dx = (upper - lower) / DistEnv.integration_slice_count
    for i in range(SLICE_COUNT):
        a = lower + i * dx
        b = lower + (i + 1) * dx
        subtotal += dx * (function(b) + function(a)) / 2
    return subtotal



# Basic functions for type-checking.
# This package implements rigorous type-checking before any arguments are passed into the functions or constructors.
# Note that for real numbers, there is no substantial difference between <IsPosNum> and <IsNonNegNum>, as well as
# between <IsNegNum> and <IsNonPosNum>, for generated random data, as zero is a point value. However, the distinction
# matters when arguments are passed explicitly into a function.
# <IsProb> checks if the value passed is a probability, which must be a real number between zero and one inclusive.
### THIS SECTION IS VALIDATED.
def isInt(x):
    return type(x) == int

def isPosInt(x):
    return type(x) == int and x > 0

def isNonNegInt(x):
    return type(x) == int and x >= 0

def isNegInt(x):
    return type(x) == int and x < 0

def isNonPosInt(x):
    return type(x) == int and x <= 0

def isNum(x):
    return type(x) == int or type(x) == float

def isPosNum(x):
    return (type(x) == int or type(x) == float) and x > 0

def isNonNegNum(x):
    return (type(x) == int or type(x) == float) and x >= 0

def isNegNum(x):
    return (type(x) == int or type(x) == float) and x < 0

def isNonPosNum(x):
    return (type(x) == int or type(x) == float) and x <= 0

# For distributions, we only consider probabilities that are not degenerate.
def isProb(x):
    return isNum(x) and x > 0 and x < 1

# Determines if the observations form a non-empty sample where all elements satisfy the predicate, usually for data type.
def isSample(observations, predicate):
    if (type(observations) is not list) or (len(observations) == 0):
        return False
    for x in observations:
        if not predicate(x):
            return False
    return True



### THIS SECTION IS VALIDATED.
class Sample:

    # Sample mean.
    @staticmethod
    def mean(elements):
        if not isSample(elements, isNum):
            raise Exception(ERR_REAL_SAMPLE)
        return sum(elements) / len(elements)

    # Sample variance.
    @staticmethod
    def var(elements):
        if not isSample(elements, isNum):
            raise Exception(ERR_REAL_SAMPLE)
        if len(elements) < 2:
            raise Exception("Your sample must have at least 2 elements for sample variance to be meaningful.")
        sample_mean = Sample.mean(elements)
        sum_of_squared_errors = 0
        for element in elements:
            sum_of_squared_errors += (element - sample_mean) ** 2
        return sum_of_squared_errors / (len(elements) - 1)

    # Sample standard deviation.
    @staticmethod
    def sd(elements):
        return sqrt(Sample.var(elements))



class Bundle:

    @staticmethod
    def getPrintout(name, parameters):
        output = "==== " + str(name) + " ===="
        for [key, value] in parameters.items():
            output += "\n# " + str(key) + " = " + str(value)
        return output


# << GENERAL FEATURES FOR PROBABILITY DISTRIBUTIONS >>
# A probability distribution is simply a rule to realize random variables.
# A family of distributions is realized as a class.
# A distribution of specific parameters is instantiated by the corresponding class constructor.
# Once instantiated, the distribution may be used repeatedly to generate random variables, or to compute
# probability masses, probability densities, or cumulative probabilities.
# < Compulsory features that all probability distributions must implement to interface with the current package >
# Regardless of the type of the distribution, the infimum and supremum of the support must be specified as the
# <inf> and <sup> members, because when the distribution is passed into some processes, a check will be performed
# to ensure that the support of the distribution is compatible with the processes. <self.inf = None> implies
# that the infimum is negative infinity, while <self.sup = None> implies that the supremum is positive infinity.
# The operator overload for <__call__> returns a random value generated from the distribution.
# The <sample> method returns an array of independent elements based on the sample size.
# The <p> method would be used to compute the probability mass for discrete distributions.
# The <f> method would be used to compute the probability density for continuous distributions.
# The <F> method would be used to compute the cumulative probability.
# The <Fc> method would be used to compute the tail probability.



# These are generic functions that makes statistical calculations look similar to their mathematical symbols.
### THIS SECTION IS VALIDATED.
class pClass:
    def __getitem__(self, rv):
        return rv.p
p = pClass()

class fClass:
    def __getitem__(self, rv):
        return rv.f
f = fClass()

class FClass:
    def __getitem__(self, rv):
        return rv.F
F = FClass()

class FcClass:
    def __getitem__(self, rv):
        return rv.Fc
Fc = FcClass()

class QClass:
    def __getitem__(self, rv):
        return rv.Q
Q = QClass()




# Computes the probability distribution for a sum of INDEPENDENT random variables.
# Currently available only for discrete rvs with support [0, +Infinity).
# VALIDATED.
class ConvolutionBase:
    def __getitem__(self, rvs):
        for rv in rvs:
            if not isinstance(rv, Distribution):
                raise Exception("You must supply a list of random variables only.")
        if len(rvs) < 2:
            raise Exception("You must supply at least two random variables.")
        return ConvolutionClass(rvs)

class ConvolutionClass:

    def __init__(self, rvs):
        self.rvs = deepcopy(rvs)

    def p(self, t):
        if not isNonNegInt(t):
            raise Exception("To compute convolution, you must supply a non-negative integer.")
        subtotal = 0
        rvs = self.rvs
        for i in range(t + 1):
            # We have already checked that <rvs> has at least 2 elements.
            if len(rvs) == 2:
                subtotal += p[rvs[0]](i) * p[rvs[1]](t - i)
            else:
                subtotal += p[self.rvs[0]](i) * Convolution[self.rvs[1:]].p(t - i)
        return subtotal

    def F(self, t):
        if not isNonNegInt(t):
            raise Exception("To compute convolution, you must supply a non-negative integer.")
        subtotal = 0
        for i in range(t + 1):
            subtotal += p[self](i)
        return subtotal

Convolution = ConvolutionBase()



### VALIDATED.
class MLEStat:

    # The Fisher information refers to estimated Fisher information.
    def __init__(self, name, estimate, sample_size, fisher_information):
        self.name = name
        self.est = estimate
        self.n = sample_size
        self.fisher = fisher_information
        if fisher_information is None:
            self.asymptotic_var = None
        else:
            self.asymptotic_var = 1 / (sample_size * fisher_information)

    def __repr__(self):
        return Bundle.getPrintout("MLE for: " + str(self.name), {
            "estimate": self.est,
            "sample size": self.n,
            "Fisher information": self.fisher,
            "asymptotic variance": self.asymptotic_var
        })

    # Confidence interval is only defined for MLEs with defined Fisher information and asymptotic variance,
    # where the distribution's support must be independent of value of parameter.
    def getCI(self, significance_level):
        if not isProb(significance_level):
            raise Exception("Significance level must be a real number between zero and one.")
        if self.asymptotic_var is None:
            return None
        else:
            z_cutoff = Q[Z](1 - significance_level / 2)
            est_sd = sqrt(self.asymptotic_var)
            return (self.est - est_sd * z_cutoff, self.est + est_sd * z_cutoff)



# Parent class for distribution, which defines methods whose implementation is the same for all distributions.
# Some distributions, due to their peculiar algorithms for sampling or probability calculation, would need to
# override the parent method with a homonymous child method.
# Kurtosis: Refers to excess kurtosis.
### VALIDATED.
class Distribution:

    def Fc(self, x):
        return 1 - self.F(x)

    def Q(self, prob):
        if not isProb(prob):
            raise Exception(ERR_PROB)
        shrinkage = 2
        # First obtain a starting point based on Chebyshev's inequality assuming the worse-case scenario.
        offset = sqrt(1 - prob)
        x = offset * self.sd + self.mean
        if abs(self.F(x) - prob) < DistEnv.quantile_precision:
            return x
        else:
            direction = ifelse(self.F(x) < prob, 1, -1)
            dx = self.sd
            while True:
                if abs(self.F(x) - prob) < DistEnv.quantile_precision:
                    return x
                if direction == 1 and self.F(x) > prob:
                    direction = -1
                    dx /= shrinkage
                elif direction == -1 and self.F(x) < prob:
                    direction = 1
                    dx /= shrinkage
                x += direction * dx

    def sample(self, sample_size):
        if not isPosInt(sample_size):
            raise Exception(ERR_SAMPLE_SIZE)
        elements = [0] * sample_size
        for i in range(sample_size):
            elements[i] = self()
        return elements



# Uniform distribution.
class Uniform(Distribution):

    @staticmethod
    def getMLE(observations):
        if not isSample(observations, isNum):
            raise Exception("Your sample must be a non-empty list where all elements are numerical.")
        sample_size = len(observations)
        lower_mle = MLEStat("lower", min(observations), sample_size, None)
        upper_mle = MLEStat("upper", max(observations), sample_size, None)
        return [lower_mle, upper_mle]

    def __init__(self, lower, upper):
        if (not isNum(lower)) or (not isNum(upper)) or lower >= upper:
            raise Exception("Lower and upper bounds must be two numerical values where lower < upper.")
        self.lower = lower
        self.upper = upper
        self.inf = lower
        self.sup = upper
        self.mean = (lower + upper) / 2
        self.var = pow(upper - lower, 2) / 12
        self.sd = sqrt(self.var)
        #self.skew = 0
        #self.kurt = - 1.2

    def __repr__(self):
        return Bundle.getPrintout("Uniform distribution", {
            "lower": self.lower,
            "upper": self.upper
        })

    def __call__(self):
        return self.lower + random() * (self.upper - self.lower)

    def f(self, x):
        if not isNum(x):
            raise Exception("Value must be a number.")
        if x < self.inf or x > self.sup:
            return 0
        else:
            return 1 / (self.upper - self.lower)

    def F(self, x):
        if not isNum(x):
            raise Exception("Value must be a number.")
        if x < self.lower:
            return 0
        elif x > self.upper:
            return 1
        else:
            return (x - self.lower) / (self.upper - self.lower)

    def Q(self, prob):
        if not isProb(prob):
            raise Exception(ERR_PROB)
        return self.lower + prob * (self.upper - self.lower)



# Generator for Bernoulli distribution.
class Ber(Distribution):

    @staticmethod
    def getMLE(observations):
        if not isSample(observations, lambda x: x == 0 or x == 1):
            raise Exception("Your sample must be a non-empty list where each element is either zero or one.")
        sample_size = len(observations)
        prob1_mle = Sample.mean(observations)
        fisher_information = 1 / (prob1_mle * (1 - prob1_mle))
        return MLEStat("prob1", prob1_mle, sample_size, fisher_information)

    def __init__(self, prob1):
        if not isProb(prob1):
            raise Exception("Probability of success must be between zero and one.")
        prob0 = 1 - prob1
        self.prob1 = prob1
        self.prob0 = prob0
        self.sup = 1
        self.inf = 0
        self.mean = prob1
        self.var = prob1 * prob0
        self.sd = sqrt(self.var)

    def __repr__(self):
        return Bundle.getPrintout("Bernoulli distribution", {
            "prob1": self.prob1
        })

    def __call__(self):
        return ifelse(random() < self.prob1, 1, 0)

    def p(self, x):
        if not isInt(x):
            raise Exception("Value must be an integer.")
        if x == 1:
            return self.prob1
        elif x == 0:
            return self.prob0
        else:
            return 0

    def F(self, x):
        if not isInt(x):
            raise Exception("Value must be an integer.")
        if x >= 1:
            return 1
        elif x >= 0:
            return self.prob0
        else:
            return 0

    def Q(self, prob):
        if not isProb(prob):
            raise Exception(ERR_PROB)
        if prob <= 1 - self.prob1:
            return 0
        else:
            return 1



class Beta(Distribution):

    def __init__(self, shape1, shape2):
        if (not isPosNum(shape1)) or (not isPosNum(shape2)):
            raise Exception("Shape parameters must both be positive numbers.")
        self.shape1 = shape1
        self.shape2 = shape2
        self.inf = 0
        self.sup = 1
        self.mean = shape1 / (shape1 + shape2)
        self.var = shape1 * shape2 / (pow(shape1 + shape2, 2) * (shape1 + shape2 + 1))
        self.sd = sqrt(self.var)

    def __repr__(self):
        return Bundle.getPrintout("Beta distribution", {
            "shape1": self.shape1,
            "shape2": self.shape2
        })

    def __call__(self):
        x = Gamma(self.shape1, 1)()
        y = Gamma(self.shape2, 1)()
        return x / (x + y)

    def f(self, x):
        if not isNum(x):
            raise Exception("Value must be a number.")
        if x < 0 or x > 1:
            return 0
        else:
            return pow(x, self.shape1 - 1) * pow(1 - x, self.shape2 - 1) / beta(self.shape1, self.shape2)

    def F(self, x):
        if not isNum(x):
            raise Exception("Value must be a number.")
        if x >= 1:
            return 1
        elif x < 0:
            return 0
        else:
            return integrate(self.f, 0, x)



# Generator for binomial distribution.
class Bin(Distribution):

    @staticmethod
    def getMLE(n, observations):
        if not isSample(observations, isNonNegInt):
            raise Exception("Your sample must be a list of non-negative integers.")
        sample_size = len(observations)
        prob1_mle = Sample.mean(observations) / n
        fisher_information = n / (prob1_mle * (1 - prob1_mle))
        return MLEStat("prob1", prob1_mle, sample_size, fisher_information)

    def __init__(self, n, prob1):
        if not isPosInt(n):
            raise Exception("Trial count must be a positive integer.")
        if not isProb(prob1):
            raise Exception("Probability of success must be between zero and one.")
        prob0 = 1 - prob1
        self.n = n
        self.prob1 = prob1
        self.prob0 = prob0
        self.sup = n
        self.inf = 0
        self.mean = n * prob1
        self.var = n * prob1 * prob0
        self.sd = sqrt(self.var)

    def __repr__(self):
        return Bundle.getPrintout("Binomial distribution", {
            "n": self.n,
            "prob1": self.prob1
        })

    def __call__(self):
        # Draws <n> Bernoulli trials and computes the sum.
        return sum(Ber(self.prob1).sample(self.n))

    def p(self, x):
        if not isInt(x):
            raise Exception("Value must be an integer.")
        if x < 0 or x > self.n:
            return 0
        else:
            # Calculates the binomial coefficient.
            coefficient = 1
            for i in range(self.n - x + 1, self.n + 1):
                coefficient *= i
            for i in range(1, x + 1):
                coefficient /= i
            return coefficient * (self.prob1 ** (x)) * (self.prob0 ** (self.n - x))

    def F(self, x):
        if not isInt(x):
            raise Exception("Value must be an integer.")
        if x < 0:
            return 0
        elif x >= self.n:
            return 1
        else:
            cumulative = 0
            for i in range(0, x + 1):
                cumulative += self.p(i)
            return cumulative

    def Q(self, prob):
        if not isProb(prob):
            raise Exception(ERR_PROB)
        for i in range(0, self.n + 2):
            if self.F(i) >= prob:
                return i



# Generator for Cauchy distribution.
class Cauchy(Distribution):

    def __init__(self, x0, scale):
        if not isNum(x0):
            raise Exception("Location parameter must be a numerical value.")
        if not isPosNum(scale):
            raise Exception("Scale parameter must be a positive number.")
        self.x0 = x0
        self.scale = scale
        self.inf = -Inf
        self.sup = Inf
        self.mean = None
        self.var = None
        self.sd = None

    def __repr__(self):
        return Bundle.getPrintout("Cauchy distribution", {
            "x0": self.x0,
            "scale": self.scale
        })

    def __call__(self):
        u = random()
        return self.x0 + self.scale * tan(pi * (u - 0.5))

    def f(self, x):
        if not isNum(x):
            raise Exception("Value must be a number.")
        return 1 / (pi * self.scale * (1 + pow((x - self.x0) / self.scale, 2)))

    def F(self, x):
        if not isNum(x):
            raise Exception("Value must be a number.")
        return (1 / pi) * atan((x - self.x0) / self.scale) + 0.5

    def Q(self, prob):
        if not isProb(prob):
            raise Exception(ERR_PROB)
        return self.x0 + self.scale * tan(pi * (prob - 0.5))



class Chi2(Distribution):

    def __init__(self, df):
        if not isPosInt(df):
            raise Exception("Degree of freedom must be a positive integer.")
        self.df = df
        self.inf = 0
        self.sup = Inf
        self.mean = df
        self.var = 2 * df
        self.sd = sqrt(self.var)

    def __repr__(self):
        return Bundle.getPrintout("Chi-square distribution", {
            "df": self.df
        })

    def __call__(self):
        subtotal = 0
        for i in range(self.df):
            subtotal += pow(Z(), 2)
        return subtotal

    def f(self, x):
        if not isNum(x):
            raise Exception("Value must be a number.")
        if x < 0:
            return 0
        else:
            return pow(x, 0.5 * self.df - 1) * exp(- x / 2) / (pow(2, self.df / 2) * gamma(self.df / 2))

    # Instead of involving the incomplete gamma function, we compute CDF by direct definite integral.
    def F(self, x):
        if not isNum(x):
            raise Exception("Value must be a number.")
        if x < 0:
            return 0
        else:
            return integrate(self.f, 0, x)



# Generator for Erlang distribution.
class Erlang(Distribution):

    @staticmethod
    def getMLE(k, observations):
        if not isSample(observations, isNonNegNum):
            raise Exception("Your sample must be a list of non-negative numbers.")
        sample_size = len(observations)
        rate_mle = k / Sample.mean(observations)
        fisher_information = 1 / (rate_mle ** 2)
        return MLEStat("rate", rate_mle, sample_size, fisher_information)

    def __init__(self, k, rate):
        if not isPosInt(k):
            raise Exception("K must be a positive integer.")
        if not isNonNegNum(rate):
            raise Exception(ERR_RATE_PARAM)
        self.k = k
        self.rate = rate
        self.inf = 0
        self.sup = Inf
        self.mean = k / rate
        self.var = k / (rate ** 2)
        self.sd = sqrt(self.var)

    def __repr__(self):
        return Bundle.getPrintout("Erlang distribution", {
            "k": self.k,
            "rate": self.rate
        })

    # Generation of Erlang-distributed random variables is based on the fact that
    # the sum of independent exponentially distributed random variables are Erlang-distributed.
    def __call__(self):
        return sum(Exp(self.rate).sample(self.k))

    def f(self, x):
        if not isNum(x):
            raise Exception("Value must be a number.")
        if x < 0:
            return 0
        else:
            return pow(self.rate, self.k) * pow(x, self.k - 1) * exp(- self.rate * x) / factorial(self.k - 1)

    # Note that for the Poisson process, Nt = Pois(lambda * t) and Sn = Erlang(n, lambda).
    # We can use the equivalence of Nt >= n <=> Sn <= t to obtain the CDF based on the associated Poisson distribution.
    def F(self, x):
        if not isNum(x):
            raise Exception("Value must be a number.")
        if x < 0:
            return 0
        else:
            poisson = Pois(self.rate * x)
            return 1 - poisson.F(self.k - 1)



# Generator for exponential distribution.
# Do not confuse the notation with the <exp> function.
class Exp(Distribution):

    @staticmethod
    def getMLE(observations):
        if not isSample(observations, isNonNegNum):
            raise Exception("Your sample must be a list of non-negative numbers.")
        sample_size = len(observations)
        rate_mle = 1 / Sample.mean(observations)
        fisher_information = 1 / (rate_mle ** 2)
        return MLEStat("rate", rate_mle, sample_size, fisher_information)

    def __init__(self, rate):
        if not isNonNegNum(rate):
            raise Exception("Rate must be a non-negative number.")
        self.rate = rate
        self.inf = 0
        self.sup = Inf
        self.mean = 1 / rate
        self.var = 1 / (rate ** 2)
        self.sd = 1 / rate

    def __repr__(self):
        return Bundle.getPrintout("Exponential distribution", {
            "rate": self.rate
        })

    # Generation of exponentially distributed random variables directly uses the quantile function.
    def __call__(self):
        u = random()
        return -ln(u) / self.rate

    def f(self, x):
        if not isNum(x):
            raise Exception("Value must be a number.")
        if x < 0:
            return 0
        else:
            return self.rate * exp(- self.rate * x)

    def F(self, x):
        if not isNum(x):
            raise Exception("Value must be a number.")
        if x < 0:
            return 0
        else:
            return 1 - exp(- self.rate * x)

    def Q(self, prob):
        if not isProb(prob):
            raise Exception(ERR_PROB)
        return - ln(1 - prob) / self.rate



# Generator for Gamma distribution:
class Gamma(Distribution):

    def __init__(self, shape, rate):
        if not isPosNum(shape):
            raise Exception("Shape parameter must be a positive number.")
        if not isPosNum(rate):
            raise Exception("Scale parameter must be a positive number.")
        self.shape = shape
        self.rate = rate
        self.inf = 0
        self.sup = Inf
        self.mean = shape / rate
        self.var = shape / (rate ** 2)
        self.sd = sqrt(self.var)
        #self.skew = 2 / sqrt(shape)
        #self.kurt = 6 / shape

    def __repr__(self):
        return Bundle.getPrintout("Gamma distribution", {
            "shape": self.shape,
            "rate": self.rate
        })

    # <delta_value> is generated by Ahrens-Dieter method.
    def __call__(self):
        n = int(self.shape)
        delta = self.shape % 1
        # Generates a random variable for Gamma(n, 1).
        if n > 0:
            n_value = sum(Exp(1).sample(n))
        else:
            n_value = 0
        # Generates a random variable for Gamma(delta, 1) if delta is non-zero.
        if delta == 0:
            delta_value = 0
        else:
            while True:
                u = random()
                v = random()
                w = random()
                if u <= e / (e + delta):
                    xi = pow(v, 1 / delta)
                    eta = w * pow(xi, delta - 1)
                else:
                    xi = 1 - ln(v)
                    eta = w * exp(-xi)
                if eta <= pow(xi, delta - 1) * exp(-xi):
                    delta_value = xi
                    break
        # (Gamma(n, 1) + Gamma(delta, 1)) / rate
        return (n_value + delta_value) / self.rate

    def f(self, x):
        if not isNum(x):
            raise Exception("Value must be a number.")
        if x < 0:
            return 0
        else:
            return (pow(self.rate, self.shape) / gamma(self.shape)) * pow(x, self.shape - 1) * exp(- self.rate * x)

    # Instead of involving the incomplete gamma function, we compute CDF by direct definite integral.
    def F(self, x):
        if not isNum(x):
            raise Exception("Value must be a number.")
        if x < 0:
            return 0
        else:
            return integrate(self.f, 0, x)



class Geo(Distribution):

    @staticmethod
    def getMLE(observations):
        if not isSample(observations, isPosInt):
            raise Exception("Your sample must be a list of positive integers.")
        sample_size = len(observations)
        prob1_mle = 1 / Sample.mean(observations)
        fisher_information = 1 / (prob1_mle ** 2) + 1 / ((1 - prob1_mle) * prob1_mle)
        return MLEStat("prob1", prob1_mle, sample_size, fisher_information)

    def __init__(self, prob1):
        if not isProb(prob1):
            raise Exception(ERR_PROB)
        self.prob1 = prob1
        self.prob0 = 1 - prob1
        self.inf = 1
        self.sup = Inf
        self.mean = 1 / prob1
        self.var = (1 - prob1) / pow(prob1, 2)
        self.sd = sqrt(self.var)

    def __repr__(self):
        return Bundle.getPrintout("Geometric distribution", {
            "prob1": self.prob1
        })

    def __call__(self):
        trial = Ber(self.prob1)
        i = 0
        while True:
            i += 1
            if trial() == 1:
                return i

    def p(self, x):
        if not isInt(x):
            raise Exception("Trial count must be an integer.")
        if x <= 0:
            return 0
        else:
            return pow(self.prob0, x - 1) * self.prob1

    def F(self, x):
        if not isInt(x):
            raise Exception("Trial count must be an integer.")
        if x <= 0:
            return 0
        else:
            return 1 - self.prob0 ** x

    def Q(self, prob):
        if not isProb(prob):
            raise Exception(ERR_PROB)
        i = 1
        while True:
            if self.F(i) >= prob:
                return i
            i += 1



# Levy distribution.
class Levy(Distribution):

    def __init__(self, loc, scale):
        if not isNum(loc):
            raise Exception(ERR_LOC_PARAM)
        if not isPosNum(scale):
            raise Exception(ERR_SCALE_PARAM)
        self.loc = loc
        self.scale = scale
        self.inf = loc
        self.sup = Inf
        self.mean = Inf
        self.var = Inf
        self.sd = Inf

    def __repr__(self):
        return Bundle.getPrintout("Levy distribution", {
            "loc": self.loc,
            "scale": self.scale
        })

    # Generation of Levy-distributed random variables is based on transformations of normal random variables.
    def __call__(self):
        std_levy = 1 / pow(Z(), 2)
        return std_levy * self.scale + self.loc

    def f(self, x):
        if not isNum(x):
            raise Exception("Value must be a number.")
        if x <= self.inf:
            return 0
        else:
            return sqrt(self.scale / (2 * pi)) * exp(-0.5 * self.scale / (x - self.loc)) / pow(x - self.loc, 1.5)

    def F(self, x):
        if not isNum(x):
            raise Exception("Value must be a number.")
        if x <= self.inf:
            return 0
        else:
            return erfc(sqrt(0.5 * self.scale / (x - self.loc)))

    def Q(self, prob):
        if not isProb(prob):
            raise Exception(ERR_PROB)
        return self.loc + self.scale / pow(Q[Z](1 - prob * 0.5), 2)



class LogNorm(Distribution):

    @staticmethod
    def getMLE(observations):
        if not isSample(observations, isPosNum):
            raise Exception("Your sample must be a list of positive numbers.")
        ln_observations = [ln(x) for x in observations]
        # Ln(X) for LogNorm is normally distributed.
        return N.getMLE(ln_observations)

    def __init__(self, norm_mean, norm_var):
        if (not isNum(norm_mean)) or (not isNum(norm_var)):
            raise Exception("Log mean and log variance must be numerical values.")
        self.norm_mean = norm_mean
        self.norm_var = norm_var
        self.norm_sd = sqrt(norm_var)
        self.inf = 0
        self.sup = Inf
        self.mean = exp(norm_mean + norm_var / 2)
        self.var = (exp(norm_var) - 1) * exp(2 * norm_mean + norm_var)
        self.sd = sqrt(self.var)

    def __repr__(self):
        return Bundle.getPrintout("Lognormal distribution", {
            "norm_mean": self.norm_mean,
            "norm_var": self.norm_var
        })

    def __call__(self):
        return exp(N(self.norm_mean, self.norm_var)())

    def f(self, x):
        if not isNum(x):
            raise Exception("Value must be a number.")
        if x <= 0:
            return 0
        else:
            return (1 / (x * self.norm_sd * sqrt(2 * pi))) * exp(- pow(ln(x) - self.norm_mean, 2) / (2 * self.norm_var))

    def F(self, x):
        if not isNum(x):
            raise Exception("Value must be a number.")
        if x <= 0:
            return 0
        else:
            return 0.5 + 0.5 * erf((ln(x) - self.norm_mean) / sqrt(2 * self.norm_var))

    def Q(self, prob):
        if not isProb(prob):
            raise Exception(ERR_PROB)
        return exp(Q[N(self.norm_mean, self.norm_var)](prob))



# Generator for normal distribution.
### VALIDATED.
class N(Distribution):

    @staticmethod
    def getMLE(observations):
        if not isSample(observations, isNum):
            raise Exception("Your sample must be a list of numbers.")
        sample_size = len(observations)
        mean_mle = Sample.mean(observations)
        var_mle = Sample.var(observations) * (sample_size - 1) / sample_size
        mean_fisher_information = 1 / var_mle
        var_fisher_information = 1 / (2 * pow(var_mle, 2))
        return [
            MLEStat("mean", mean_mle, sample_size, mean_fisher_information),
            MLEStat("var", var_mle, sample_size, var_fisher_information)
        ]

    def __init__(self, mean, var):
        if not isNum(mean):
            raise Exception("Distribution mean must be a numerical value.")
        if not isPosNum(var):
            raise Exception("Distribution variance must be a positive number.")
        self.inf = -Inf
        self.sup = Inf
        self.mean = mean
        self.var = var
        self.sd = sqrt(var)

    def __repr__(self):
        return Bundle.getPrintout("Normal distribution", {
            "mean": self.mean,
            "var": self.var
        })

    def __call__(self):
        u = random()
        v = random()
        return (sqrt(- 2 * ln(u)) * cos(2 * pi * v)) * self.sd + self.mean

    # Generation of normally distributed random variables are based on Box-Muller method.
    # Reference: https://en.wikipedia.org/wiki/Normal_distribution.
    # Note that Box-Muller method generates two random variables at each time, and thus is not compatible with
    # the sampling method of the parent class. Therefore, overriding is required.
    # Although we can run <__call__> to generate observations one by one, this is clearly not efficient.
    def sample(self, sample_size):
        if not isPosInt(sample_size):
            raise Exception(ERR_SAMPLE_SIZE)
        elements = []
        for i in range(sample_size // 2):
            u = random()
            v = random()
            z1 = sqrt(- 2 * ln(u)) * cos(2 * pi * v)
            z2 = sqrt(- 2 * ln(u)) * sin(2 * pi * v)
            elements += [z1 * self.sd + self.mean, z2 * self.sd + self.mean]
        if sample_size % 2 == 1:
            u = random()
            v = random()
            elements += [(sqrt(- 2 * ln(u)) * cos(2 * pi * v)) * self.sd + self.mean]
        return elements

    def f(self, x):
        if not isNum(x):
            raise Exception("Value must be a number.")
        return (1 / (self.sd * sqrt(2 * pi))) * exp(- 0.5 * pow(x - self.mean, 2) / self.var)

    def F(self, x):
        return 0.5 * (1 + erf((x - self.mean) / (self.sd * sqrt(2))))

    def Q(self, prob):
        if not isProb(prob):
            raise Exception(ERR_PROB)
        return self.mean + self.sd * sqrt(2) * erfinv(2 * prob - 1)



### VALIDATED.
class NegBin(Distribution):

    @staticmethod
    def getMLE(k, observations):
        if not isSample(observations, isPosInt):
            raise Exception("Your sample must be a list of positive integers.")
        sample_size = len(observations)
        prob1_mle = k / Sample.mean(observations)
        fisher_information = k / (prob1_mle ** 2) + k / (prob1_mle * (1 - prob1_mle))
        return MLEStat("prob1", prob1_mle, sample_size, fisher_information)

    def __init__(self, k, prob1):
        self.k = k
        self.prob1 = prob1
        self.prob0 = 1 - prob1
        self.inf = k
        self.sup = Inf
        self.mean = k / prob1
        self.var = k * (1 - prob1) / pow(prob1, 2)
        self.sd = sqrt(self.var)

    def __repr__(self):
        return Bundle.getPrintout("Negative binomial distribution", {
            "k": self.k,
            "prob1": self.prob1
        })

    def __call__(self):
        return sum(Geo(self.prob1).sample(self.k))

    def p(self, x):
        if not isInt(x):
            raise Exception("Trial count must be an integer.")
        if x < self.inf:
            return 0
        # Computes (x - 1)C(k - 1) = (x - 1)! / [(k - 1)! * (x - k)!
        coeff = 1
        for i in range(1, self.k):
            coeff *= x - i
        coeff /= factorial(self.k - 1)
        return coeff * pow(self.prob1, self.k) * pow(self.prob0, x - self.k)

    def F(self, x):
        if not isInt(x):
            raise Exception("Trial count must be an integer.")
        if x < self.inf:
            return 0
        else:
            cumulative = 0
            for i in range(self.k, x + 1):
                cumulative += self.p(i)
            return cumulative

    def Q(self, prob):
        if not isProb(prob):
            raise Exception(ERR_PROB)
        i = 1
        while True:
            if self.F(i) >= prob:
                return i
            i += 1



# Generator for Poisson distribution.
class Pois(Distribution):

    @staticmethod
    def getMLE(observations):
        if not isSample(observations, isNonNegInt):
            raise Exception("Your sample must be a list of non-negative integers.")
        sample_size = len(observations)
        rate_mle = Sample.mean(observations)
        fisher_information = 1 / rate_mle
        return MLEStat("rate", rate_mle, sample_size, fisher_information)

    def __init__(self, rate):
        if not isNonNegNum(rate):
            raise Exception("Rate must be a non-negative number.")
        self.rate = rate
        self.inf = 0
        self.sup = Inf
        self.mean = rate
        self.var = rate
        self.sd = sqrt(rate)

    def __repr__(self):
        return Bundle.getPrintout("Poisson distribution", {
            "rate": self.rate
        })

    # Generation of Poisson-distributed random variables are based on inverse transform sampling.
    # Reference: https://en.wikipedia.org/wiki/Poisson_distribution.
    def __call__(self):
        rate = self.rate
        x = 0
        p = exp(-rate)
        s = p
        u = random()
        while u > s:
            x += 1
            p *= (rate / x)
            s += p
        return x

    def p(self, x):
        if not isInt(x):
            raise Exception("Value must be an integer.")
        if x < 0:
            return 0
        else:
            return (self.rate ** x) * exp(-self.rate) / factorial(x)

    def F(self, x):
        if not isInt(x):
            raise Exception("Value must be an integer.")
        if x < 0:
            return 0
        else:
            cumulative = 0
            for i in range(0, x + 1):
                cumulative += self.p(i)
            return cumulative

    def Q(self, prob):
        if not isProb(prob):
            raise Exception(ERR_PROB)
        i = 0
        while True:
            if self.F(i) >= prob:
                return i
            i += 1



# To avoid confusion with the CDF notation F, we use the full name of the F-distribution.
class SnedecorF(Distribution):

    def __init__(self, df1, df2):
        if (not isPosNum(df1)) or (not isPosNum(df2)):
            raise Exception("Degrees of freedom must be positive numbers.")
        self.df1 = df1
        self.df2 = df2
        self.inf = 0
        self.sup = Inf
        if df2 > 2:
            self.mean = df2 / (df2 - 2)
        else:
            self.mean = None
        if df2 > 4:
            self.var = 2 * pow(df2, 2) * (df1 + df2 - 2) / (df1 * pow(df2 - 2, 2) * (df2 - 4))
        else:
            self.var = None
        if self.var is None:
            self.sd = None
        else:
            self.sd = sqrt(self.var)

    def __repr__(self):
        return Bundle.getPrintout("F-distribution", {
            "df1": self.df1,
            "df2": self.df2
        })

    def __call__(self):
        c1 = Chi2(self.df1)()
        c2 = Chi2(self.df2)()
        return (c1 / self.df1) / (c2 / self.df2)

    def f(self, x):
        if not isNum(x):
            raise Exception("Value must be a number.")
        if x <= 0:
            return 0
        df1 = self.df1
        df2 = self.df2
        numerator = pow(df1 / df2, 0.5 * df1) * pow(x, 0.5 * df1 - 1) * pow(1 + x * df1 / df2, -0.5 * (df1 + df2))
        denominator = beta(df1 / 2, df2 / 2)
        return numerator / denominator

    def F(self, x):
        if not isNum(x):
            raise Exception("Value must be a number.")
        if x <= 0:
            return 0
        else:
            return integrate(self.f, 0, x)



class t(Distribution):

    def __init__(self, df):
        if not isPosNum(df):
            raise Exception("Degree of freedom must be a positive number.")
        self.df = df
        self.inf = -Inf
        self.sup = Inf
        if df > 1:
            self.mean = 0
        else:
            self.mean = None
        if df > 2:
            self.var = df / (df - 2)
        else:
            self.var = None
        if self.var is None:
            self.sd = None
        else:
            self.sd = sqrt(self.var)

    def __repr__(self):
        return Bundle.getPrintout("t-distribution", {
            "df": self.df
        })

    def __call__(self):
        return Z() / sqrt(Chi2(self.df)() / self.df)
        #n = self.df + 1
        #elements = Z.sample(n)
        #return Sample.mean(elements) * sqrt(n) / Sample.sd(elements)

    def f(self, x):
        if not isNum(x):
            raise Exception("Value must be a number.")
        coeff_term = gamma((self.df + 1) / 2) / (sqrt(self.df * pi) * gamma(self.df / 2))
        power_term = pow(1 + pow(x, 2) / self.df, -(self.df + 1) / 2)
        return coeff_term * power_term

    def F(self, x):
        if not isNum(x):
            raise Exception("Value must be a number.")
        if x >= 0:
            return 0.5 + integrate(self.f, 0, x)
        else:
            return 0.5 - integrate(self.f, x, 0)



### < SECTION ON HYPOTHESIS TESTING >
# For hypothesis testing, p-values would be returned so that interpretations can be made on different significance levels.
# Each hypothesis testing returns a TestStat object that contains the statistic and p-value.
class TestStat:

    def __init__(self, type, hypothesized_value, observed_value, test_stat, p_value):
        self.type = type
        self.hypothesized_value = hypothesized_value
        self.observed_value = observed_value
        self.test_stat = test_stat
        self.p_value = p_value

    def __repr__(self):
        return Bundle.getPrintout(self.type, {
            "hypothesized": self.hypothesized_value,
            "observed": self.observed_value,
            "statistic": self.test_stat,
            "p-value": self.p_value
        })



def zTest1s1t(mean_h0, var, observations):
    if not isNum(mean_h0):
        raise Exception("Null hypothesis must have a numerial mean value.")
    if not isPosNum(var):
        raise Exception("Population variance must be a positive number.")
    if not isSample(observations, isNum):
        raise Exception("Your sample must be a list of numbers.")
    sample_mean = Sample.mean(observations)
    z_value = (sample_mean - mean_h0) / sqrt(var / len(observations))
    p_value = ifelse(sample_mean > mean_h0, Fc[Z](z_value), F[Z](z_value))
    return TestStat("z-test, one-tailed", mean_h0, sample_mean, z_value, p_value)



def zTest1s2t(mean_h0, var, observations):
    if not isNum(mean_h0):
        raise Exception("Null hypothesis must have a numerial mean value.")
    if not isPosNum(var):
        raise Exception("Population variance must be a positive number.")
    if not isSample(observations, isNum):
        raise Exception("Your sample must be a list of numbers.")
    sample_mean = Sample.mean(observations)
    z_value = (sample_mean - mean_h0) / sqrt(var / len(observations))
    p_value = ifelse(sample_mean > mean_h0, 2 * Fc[Z](z_value), 2 * F[Z](z_value))
    return TestStat("z-test, two-tailed", mean_h0, sample_mean, z_value, p_value)



def zTest2s1t(var1, var2, observations1, observations2):
    if (not isPosNum(var1)) or (not isPosNum(var2)):
        raise Exception("Population variances for both samples must be positive numbers.")
    if (not isSample(observations1, isNum)) or (not isSample(observations2, isNum)):
        raise Exception("Both samples must be lists of numerical values.")
    sample_size_1 = len(observations1)
    sample_mean_1 = Sample.mean(observations1)
    sample_size_2 = len(observations2)
    sample_mean_2 = Sample.mean(observations2)
    mean_diff = sample_mean_1 - sample_mean_2
    z_value = mean_diff / sqrt(var1 / sample_size_1 + var2 / sample_size_2)
    p_value = ifelse(mean_diff > 0, Fc[Z](z_value), F[Z](z_value))
    return TestStat("two sample z-test, one-tailed", 0, mean_diff, z_value, p_value)



def zTest2s2t(var1, var2, observations1, observations2):
    if (not isPosNum(var1)) or (not isPosNum(var2)):
        raise Exception("Population variances for both samples must be positive numbers.")
    if (not isSample(observations1, isNum)) or (not isSample(observations2, isNum)):
        raise Exception("Both samples must be lists of numerical values.")
    sample_size_1 = len(observations1)
    sample_mean_1 = Sample.mean(observations1)
    sample_size_2 = len(observations2)
    sample_mean_2 = Sample.mean(observations2)
    mean_diff = sample_mean_1 - sample_mean_2
    z_value = mean_diff / sqrt(var1 / sample_size_1 + var2 / sample_size_2)
    p_value = ifelse(mean_diff > 0, 2 * Fc[Z](z_value), 2 * F[Z](z_value))
    return TestStat("two sample z-test, two-tailed", 0, mean_diff, z_value, p_value)



def tTest1s1t(mean_h0, observations):
    if not isNum(mean_h0):
        raise Exception("Null hypothesis must have a numerial mean value.")
    if not isSample(observations, isNum):
        raise Exception("Your sample must be a list of numbers.")
    sample_size = len(observations)
    df = sample_size - 1
    sample_mean = Sample.mean(observations)
    sample_var = Sample.var(observations)
    t_value = (sample_mean - mean_h0) / sqrt(sample_var / sample_size)
    p_value = ifelse(sample_mean > mean_h0, Fc[t(df)](t_value), F[t(df)](t_value))
    return TestStat("Student's t-test, one-tailed", mean_h0, sample_mean, t_value, p_value)



def tTest1s2t(mean_h0, observations):
    if not isNum(mean_h0):
        raise Exception("Null hypothesis must have a numerial mean value.")
    if not isSample(observations, isNum):
        raise Exception("Your sample must be a list of numbers.")
    sample_size = len(observations)
    df = sample_size - 1
    sample_mean = Sample.mean(observations)
    sample_var = Sample.var(observations)
    t_value = (sample_mean - mean_h0) / sqrt(sample_var / sample_size)
    p_value = ifelse(sample_mean > mean_h0, 2 * Fc[t(df)](t_value), 2 * F[t(df)](t_value))
    return TestStat("Student's t-test, two-tailed", mean_h0, sample_mean, t_value, p_value)



def tTest2s1t(observations1, observations2, same_var = False):
    if (not isSample(observations1, isNum)) or (not isSample(observations2, isNum)):
        raise Exception("Both samples must be lists of numerical values.")
    sample_size_1 = len(observations1)
    sample_mean_1 = Sample.mean(observations1)
    sample_var_1 = Sample.var(observations1)
    sample_size_2 = len(observations2)
    sample_mean_2 = Sample.mean(observations2)
    sample_var_2 = Sample.var(observations2)
    mean_diff = sample_mean_1 - sample_mean_2
    if same_var:
        overall_var = ((sample_size_1 - 1) * sample_var_1 + (sample_size_2 - 1) * sample_var_2) / (sample_size_1 + sample_size_2 - 2)
        t_value = mean_diff / sqrt(overall_var / sample_size_1 + overall_var / sample_size_2)
        df = sample_size_1 + sample_size_2 - 2
    else:
        t_value = mean_diff / sqrt(sample_var_1 / sample_size_1 + sample_var_2 / sample_size_2)
        df_numerator = pow(sample_var_1 / sample_size_1 + sample_var_2 / sample_size_2, 2)
        df_denominator = pow(sample_var_1 / sample_size_1, 2) / (sample_size_1 - 1) + \
            pow(sample_var_2 / sample_size_2, 2) / (sample_size_2 - 1)
        df = df_numerator / df_denominator
    p_value = ifelse(mean_diff > 0, Fc[t(df)](t_value), F[t(df)](t_value))
    return TestStat("two-sample t-test, one-tailed", 0, mean_diff, t_value, p_value)



def tTest2s2t(observations1, observations2, same_var = False):
    if (not isSample(observations1, isNum)) or (not isSample(observations2, isNum)):
        raise Exception("Both samples must be lists of numerical values.")
    sample_size_1 = len(observations1)
    sample_mean_1 = Sample.mean(observations1)
    sample_var_1 = Sample.var(observations1)
    sample_size_2 = len(observations2)
    sample_mean_2 = Sample.mean(observations2)
    sample_var_2 = Sample.var(observations2)
    mean_diff = sample_mean_1 - sample_mean_2
    if same_var:
        overall_var = ((sample_size_1 - 1) * sample_var_1 + (sample_size_2 - 1) * sample_var_2) / (sample_size_1 + sample_size_2 - 2)
        t_value = mean_diff / sqrt(overall_var / sample_size_1 + overall_var / sample_size_2)
        df = sample_size_1 + sample_size_2 - 2
    else:
        t_value = mean_diff / sqrt(sample_var_1 / sample_size_1 + sample_var_2 / sample_size_2)
        df_numerator = pow(sample_var_1 / sample_size_1 + sample_var_2 / sample_size_2, 2)
        df_denominator = pow(sample_var_1 / sample_size_1, 2) / (sample_size_1 - 1) + \
            pow(sample_var_2 / sample_size_2, 2) / (sample_size_2 - 1)
        df = df_numerator / df_denominator
    p_value = ifelse(mean_diff > 0, 2 * Fc[t(df)](t_value), 2 * F[t(df)](t_value))
    return TestStat("two-sample t-test, two-tailed", 0, mean_diff, t_value, p_value)



def chi2TestVar1t(var_h0, observations):
    if not isPosNum(var_h0):
        raise Exception("Null hypothesis must have a positive number for variance.")
    if (not isSample(observations, isNum)) or len(observations) < 2:
        raise Exception("Your sample must be a list of at least two numbers.")
    df = len(observations) - 1
    sample_var = Sample.var(observations)
    chi2_value = df * sample_var / var_h0
    p_value = ifelse(sample_var > var_h0, Fc[Chi2(df)](chi2_value), F[Chi2(df)](chi2_value))
    return TestStat("Chi-square test of variance, one-tailed", var_h0, sample_var, chi2_value, p_value)



def chi2TestVar2t(var_h0, observations):
    if not isPosNum(var_h0):
        raise Exception("Null hypothesis must have a positive number for variance.")
    if (not isSample(observations, isNum)) or len(observations) < 2:
        raise Exception("Your sample must be a list of at least two numbers.")
    df = len(observations) - 1
    sample_var = Sample.var(observations)
    chi2_value = df * sample_var / var_h0
    p_value = ifelse(sample_var > var_h0, 2 * Fc[Chi2(df)](chi2_value), 2 * F[Chi2(df)](chi2_value))
    return TestStat("Chi-square test of variance, two-tailed", var_h0, sample_var, chi2_value, p_value)



def fTestVar1t(observations1, observations2):
    if (not isSample(observations1, isNum)) or (not isSample(observations2, isNum)):
        raise Exception("Both samples must be lists of numerical values.")
    if len(observations1) < 2 or len(observations2) < 2:
        raise Exception("Both samples must have sample size of at least two.")
    f_value = Sample.var(observations1) / Sample.var(observations2)
    df1 = len(observations1) - 1
    df2 = len(observations2) - 1
    f_dist = SnedecorF(df1, df2)
    p_value = ifelse(f_value > 1, Fc[f_dist](f_value), F[f_dist](f_value))
    return TestStat("F-test of variance, one-tailed", 1, f_value, f_value, p_value)



def fTestVar2t(observations1, observations2):
    if (not isSample(observations1, isNum)) or (not isSample(observations2, isNum)):
        raise Exception("Both samples must be lists of numerical values.")
    if len(observations1) < 2 or len(observations2) < 2:
        raise Exception("Both samples must have sample size of at least two.")
    f_value = Sample.var(observations1) / Sample.var(observations2)
    df1 = len(observations1) - 1
    df2 = len(observations2) - 1
    f_dist = SnedecorF(df1, df2)
    p_value = ifelse(f_value > 1, 2 * Fc[f_dist](f_value), 2 * F[f_dist](f_value))
    return TestStat("F-test of variance, two-tailed", 1, f_value, f_value, p_value)



# Standard normal distribution.
Z = N(0, 1)