import pandas as pd
import matplotlib.pyplot as plt
import math


# Change this directory to the location of the data file
fileLocation = 'C:\\Users\\Daniel F\\Downloads\\Consolidated Clean Data.xlsx'


def lin_uncertainty(arr_y, arr_x):  # this works

    N = len(arr_y)  # arr_y is of the same length as arr_y

    term_00 = 0  # first thing in the delta eqn (just the sums)
    term_01 = 0  # second thing in the delta eqn (just the sums)
    term_10 = 0  # product of x_i and y_i, gets summed
    term_12 = 0  # sum of y vals

    for i in range(0, N):
        term_00 += arr_x[i] ** 2
        term_01 += arr_x[i]
        term_10 += arr_x[i] * arr_y[i]
        term_11 = term_01  # they're the same thing, just for consistency
        term_12 += arr_y[i]

    delta = N * term_00 - term_01 ** 2

    m = (N * term_10 - term_11 * term_12) / delta
    b = (term_12 - m * term_01)/N

    return [m, b]


def linApprox(attx, atty):

    plt.plot(attx, atty)
    plt.show()

    print("For slope computations, pick data points that do not have repeating x or y values.")

    lb = float(input("Input linear region lower bound on x data: "))
    ub = float(input("Input linear region upper bound on x data: "))

    lbx = min(attx, key=lambda x: abs(x - lb))
    ubx = min(attx, key=lambda x: abs(x - ub))

    num = attx.index(ubx) - attx.index(lbx)
    print("Number of samples used: " + str(num))

    return lin_uncertainty(atty[attx.index(lbx) : attx.index(ubx) + 1], attx[attx.index(lbx) : attx.index(ubx) + 1])


class sample:

    def __init__(self, MTSForce, laserDisplacement, strainGauge1, strainGauge2, area):

        self.MTSForce = MTSForce.tolist()
        self.laserDisplacement = laserDisplacement.tolist()
        self.strainGauge1 = strainGauge1.tolist()
        self.strainGauge2 = strainGauge2.tolist()
        self.area = area

    def stress(self):

        sigma = []

        for i in range(0, len(self.MTSForce)):

            sigma.append(self.MTSForce[i] / self.area)

        return sigma

    def yieldStress(self):

        epsilonOffset = []
        sigmaExtrap = []
        res = []
        minimum = 10e15
        minInd = -1

        for i in range(0, len(self.strainGauge1)):

            epsilonOffset.append(self.strainGauge1[i] + 0.002)

        [m, b] = linApprox(epsilonOffset, self.stress())

        indexOffset = self.strainGauge1.index(min(self.strainGauge1, key=lambda x: abs(x - epsilonOffset[0])))

        for j in range(0, len(self.strainGauge1) - indexOffset):

            sigmaExtrap.append(m * epsilonOffset[j] + b)
            res.append(abs(self.stress()[j + indexOffset] - sigmaExtrap[j]))

            if res[j] < minimum:

                minimum = res[j]
                minInd = j

        plt.plot(self.strainGauge1, self.stress(), epsilonOffset[:len(epsilonOffset) - indexOffset], sigmaExtrap)
        plt.show()

        return self.stress()[minInd]

    def elasticModulus(self):

        return linApprox(self.strainGauge1, self.stress())[0]

    def ultimate(self):

        return max(self.stress())

    def elongation(self):

        return max(self.strainGauge1)

    def rupture(self):

        return self.stress()[self.strainGauge1.index(self.elongation())]


# Collect data frames of each sample
sample1df = pd.read_excel(fileLocation, sheet_name = 0)
sample2df = pd.read_excel(fileLocation, sheet_name = 1)
sample3df = pd.read_excel(fileLocation, sheet_name = 2)
sample4df = pd.read_excel(fileLocation, sheet_name = 3)
sample5df = pd.read_excel(fileLocation, sheet_name = 4)


# Convert data frames to sample instances
sample1 = sample(sample1df.loc[:,"MTS force"], sample1df.loc[:,"laser displacement"], sample1df.loc[:,"strain gauge 1"],
                 sample1df.loc[:,"strain gauge 2"], 0.000045)
sample2 = sample(sample2df.loc[:,"MTS force"], sample2df.loc[:,"laser displacement"], sample2df.loc[:,"strain gauge 1"],
                 sample2df.loc[:,"strain gauge 2"], 0.000045)
sample3 = sample(sample3df.loc[:,"MTS force"], sample3df.loc[:,"laser displacement"], sample3df.loc[:,"strain gauge 1"],
                 sample3df.loc[:,"strain gauge 2"], 0.000045)
sample4 = sample(sample4df.loc[:,"MTS force"], sample4df.loc[:,"laser displacement"], sample4df.loc[:,"strain gauge 1"],
                 sample4df.loc[:,"strain gauge 2"], 0.000045)
sample5 = sample(sample5df.loc[:,"MTS force"], sample5df.loc[:,"laser displacement"], sample5df.loc[:,"strain gauge 1"],
                 sample5df.loc[:,"strain gauge 2"], 0.000045)

if __name__ == "__main__":

    print(sample5.ultimate())