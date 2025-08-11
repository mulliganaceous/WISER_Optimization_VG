import dimod
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from dimod import BinaryQuadraticModel

from types import NoneType
from typing import Union

# Static variables for guardrails
K_LOW = 0
K_TARGET = 1
K_HIGH = 2
B_LOW = 3
B_HIGH = 4

# Static variables for security tuples
CID = 0
MINTRADE = 1
MAXTRADE = 2
INVENTORY = 3
PMV = 4
DXS = 5
DELTA = 6
RATING = 7
BETAS = -1

# Static variables for outputs
ENDCOLOR = '\033[0m'
HEADER = '\033[5m\033[1;32m'
YELLOW = '\033[0;33m'
IN = "\033[0;32m"
OUT = "\033[0;31m"
OUT_LOW = '\033[0;31m'
OUT_HIGH = '\033[0;35m'

# Print methods for outputs
def _rangecolor_binary(lb, val, ub):
    if val < lb or val > ub:
        return OUT
    return IN

def _rangecolor(lb, val, ub):
    if np.isnan(val):
        return "\033[0;30m"
    if val < lb:
        return OUT_LOW
    if val > ub:
        return OUT_HIGH
    return IN

def _printresult(a_idx, s):
    if a_idx < 5:
        print(s)

# Utility methods
def x(c):
    """ Provided a Security, obtain the max as determined by the simplified OneOpt formula. """
    if not c[DELTA]:
        raise Exception("A bond has zero interval.")
    return (c[MINTRADE] + min(c[MAXTRADE], c[INVENTORY]))/(2*c[DELTA])

class SimplifiedOneOpto():
    """ Ad-hoc object for performing simplified OneOpto optimization, without automatic data processing. """
    def __init__(self, 
                 C: list[tuple],
                 KC: list[list[tuple]], 
                 indexmap: dict,
                 L: list[str], 
                 J: list[str], 
                 KB: list[list[list[float]]], 
                 N: int, 
                 R: list[float],
                 beta: Union[NoneType, list[list[float]]] = None, 
                 rho: Union[NoneType, list[float]] = None):
        
        # Checks for sanity
        if len(KB) != len(J):
            raise Exception("K/B matrix requires {0:d} set of guardrail sets; there are {1:d}".format(len(J), len(KB)))
        for jIdx in range(len(KB)):
            if len(KB[jIdx]) != len(L):
                raise Exception("K/B[{2:s}] by characteristic requires {0:d} set of guardrails; there are {1:d}".format(len(L), len(KB[jIdx]), J[jIdx]))
            for lIdx in range(len(KB[jIdx])):
                if len(KB[jIdx][lIdx]) < 5:
                    raise Exception("Guardrail information of K/B[{0:s},{1:s}] incomplete.".format(J[jIdx], L[lIdx]))
                
        # Set as instance variables
        self.C, self.KC, self.indexmap = C, KC, indexmap
        self.L = L
        self.J = J
        self.KB = KB
        self.beta = beta
        self.N = N
        self.R = R
        if not beta:
            self.beta = np.ones((len(J), len(C)))
        if not rho:
            self.rho = [1]*len(J)

        # Generate QUBO matrix
        matrix = np.zeros((len(C), len(C)))
        for lIdx in range(len(KC)):
            for jIdx in range(len(J)):
                # Obtain quadratic coefficients, including diagonal terms
                KCl = KC[lIdx]
                for c1 in KCl:
                    for c2 in KCl:
                        i1 = indexmap[c1[CID]]
                        i2 = indexmap[c2[CID]]
                        matrix[i1][i2] += self.rho[jIdx]*self.beta[jIdx][i1]*c1[PMV]*x(c1)*self.beta[jIdx][i2]*c2[PMV]*x(c2)
                # Obtain linear coefficients
                for c in KCl:
                    i0 = indexmap[c[CID]]
                    matrix[indexmap[c[CID]]][indexmap[c[CID]]] -= self.rho[jIdx]*2*KB[jIdx][lIdx][K_TARGET]*self.beta[jIdx][i0]*c1[PMV]*x(c)

        # Consider constant term
        constantterm = np.zeros((len(J), len(L)))
        for lIdx in range(len(L)):
            for jIdx in range(len(J)):
                constantterm[jIdx][lIdx] = KB[jIdx][lIdx][K_TARGET]*KB[jIdx][lIdx][K_TARGET]

        self.matrix = matrix
        self.constantterm = constantterm
        self.results = None

    def run(self, num_reads):
        # Convert to binary quadratic model
        bqm = BinaryQuadraticModel(vartype=dimod.BINARY)
        for i in range(len(C)):
            bqm.add_variable(C[i][CID])
            bqm.add_linear(C[i][CID], self.matrix[i][i])
        for i1 in range(len(C)):
            for i2 in range(i1 + 1, len(C)):
                if self.matrix[i1][i2]:
                    bqm.add_quadratic(C[i1][CID], C[i2][CID], 2*self.matrix[i1][i2])
        # Add master linear inequality constraints
        sy = []
        sx = []
        for i in range(len(self.C)):
            sy.append((self.C[i][CID], 1))
            sx.append((self.C[i][CID], self.C[i][PMV]*self.C[i][DELTA]*x(self.C[i])))

        lm = np.sum(abs(self.matrix), axis=1) # Hard restraint
        bqm.add_linear_inequality_constraint(sy, ub = self.N, lagrange_multiplier=np.dot(lm,lm), label="Sy") # Bonds in basket
        bqm.add_linear_inequality_constraint(sx, lb = self.R[0], ub = self.R[1], lagrange_multiplier=np.dot(lm,lm), label="Sx") # Residual cash flow

        # # Add bucket and characteristic based inequality constraints
        for j in range(len(self.J)):
            for l in range(len(self.L)):
                if len(self.KC[l]) == 0:
                    continue 
                syl = []
                sxl = []
                for i in range(len(self.KC[l])):
                    c = self.KC[l][i]
                    i0 = self.indexmap[c[CID]]
                    sxl.append((c[CID], c[PMV]*c[DELTA]*self.beta[j][i0]*x(c)))
                print(l, j, self.KB[j][l])
                # print(syl)
                print(sxl)
                # bqm.add_linear_inequality_constraint(syl, lb = KB[j][l][K_LOW], ub = KB[j][l][K_HIGH], lagrange_multiplier=np.sqrt(np.dot(lm,lm)), label="Sy_{0:s}_{1:s}".format(J[j], L[l]))
                bqm.add_linear_inequality_constraint(sxl, lb = self.KB[j][l][K_LOW], ub = self.KB[j][l][K_HIGH], lagrange_multiplier=np.dot(lm,lm), label="Sx_{0:s}_{1:s}".format(self.J[j], self.L[l]))

        from dwave.samplers import PathIntegralAnnealingSampler
        sampler = PathIntegralAnnealingSampler()
        sampleset = sampler.sample(bqm, num_reads=num_reads)
        df_s = sampleset.to_pandas_dataframe().sort_values(by="energy")
        df_s.to_csv("dwave_result.csv")
        self.results = df_s
        return df_s
    
    def display_results(self, CUTOFF, prefix1="results", prefix2="resultsK", suffix="OMG"):
        df_s = self.results
        answers = df_s[0:CUTOFF].to_numpy()[:,:len(C)]
        totalenergies = df_s[0:CUTOFF].to_numpy()[:,-2]
        energies = np.zeros(CUTOFF)
        cashflow = np.zeros(CUTOFF)
        characteristic = np.zeros((len(L), len(J), CUTOFF))
        c_ids = np.array(list(indexmap.keys()))

        # Gather and display text data
        for a_idx in range(min(len(answers),CUTOFF)):
            answer = answers[a_idx]
            totalenergy = np.trunc(totalenergies[a_idx]).astype(int)
            ids = c_ids[answer.astype(bool)]
            _printresult(a_idx, f"{HEADER}Solution #1 ({totalenergy}): {ids}{ENDCOLOR}")
            energy = (answer @ self.matrix @ answer) + np.sum(self.constantterm)
        
            energies[a_idx] = energy
            bonds = np.dot(np.ones(len(C)), answer).astype(int)

            for c_id in ids:
                c = C[indexmap[c_id]]
                xh = x(c)
                coeff = c[PMV]*c[DELTA]*x(c)
                cashflow[a_idx] += coeff
                _printresult(a_idx, "\t{0:.3f}\t{1:.3f}\t{2:.3f}\t{3:s}".format(xh, c[PMV], coeff, str(c)))
            _printresult(a_idx, YELLOW + "\tEnergy:\t{0:.3f}".format(energy) + ENDCOLOR)
            _printresult(a_idx, _rangecolor_binary(0,bonds,N) + "\tBonds:\t{0:d}".format(bonds) + ENDCOLOR)
            _printresult(a_idx, _rangecolor_binary(R[0],cashflow[a_idx],R[1]) + "\tFlow:\t{0:.3f}".format(cashflow[a_idx]) + ENDCOLOR)

            for lIdx in range(len(L)):
                for jIdx in range(len(J)):
                    KCl = KC[lIdx]
                    for c in KCl:
                        if not c[CID] in ids:
                            continue
                        xh = x(c)*(c[CID] in ids)
                        kcoeff = c[PMV]*c[DELTA]*self.beta[jIdx][self.indexmap[c[CID]]]*x(c)
                        characteristic[lIdx, jIdx, a_idx] += kcoeff
                        # bcoeff = c[DXS]/100*c[DELTA]*c[-1][jIdx]*how_much_bond(c)
                        # bucketflow += bcoeff
                        _printresult(a_idx, "\t\t{0:.3f}\t{1:.3f}\t{2:.3f}\t{3:s}".format(xh, c[PMV], kcoeff, str(c)))
                    k_diff = characteristic[lIdx, jIdx, a_idx] - KB[jIdx][lIdx][K_TARGET]
                    _printresult(a_idx, _rangecolor(KB[jIdx][lIdx][K_LOW], characteristic[lIdx, jIdx, a_idx], KB[jIdx][lIdx][K_HIGH]) 
                                + "\t\tBonds[{1:s}][{2:s}]: Σ={0:.3f} Δ={6:.3f} Δ²={7:.3f} ({3:.3f}, {4:.3f}, {5:.3f})"
                                .format(characteristic[lIdx, jIdx, a_idx], L[lIdx], J[jIdx], 
                                        KB[jIdx][lIdx][K_LOW], KB[jIdx][lIdx][K_TARGET], KB[jIdx][lIdx][K_HIGH],
                                        k_diff, k_diff*k_diff) + ENDCOLOR)
                    # _printresult(a_idx, _rangecolor(KB[jIdx][lIdx][B_LOW], bucketflow, KB[jIdx][lIdx][B_HIGH]) 
                    #       + "\tFlow[{1:s}][{2:s}]:  {0:.3f} ({3:.3f}, {4:.3f})".format(
                    #           bucketflow, L[lIdx], J[jIdx], KB[jIdx][lIdx][B_LOW], KB[jIdx][lIdx][B_HIGH]) + ENDCOLOR)
            _printresult(a_idx, "")

        # Obtain graphical charts
        T = np.linspace(0, CUTOFF, CUTOFF)
        fig = plt.subplots(3,1, figsize=(8, 9))

        plt.subplot(3,1,1)
        plt.title("N={0:d}, R={1:s}".format(N, str(R)))
        plt.plot(np.linspace(0, CUTOFF, CUTOFF), totalenergies)
        plt.plot(np.linspace(0, CUTOFF, CUTOFF), energies)
        plt.gca().set(xlim=(0,CUTOFF), xlabel="rank")
        plt.gca().set(ylabel="energy")

        plt.subplot(3,1,2)
        plt.plot(np.linspace(0, CUTOFF, CUTOFF), cashflow)
        plt.gca().fill_between(np.linspace(0, CUTOFF, CUTOFF), R[1], R[0], alpha=0.36)
        plt.gca().set(xlim=(0,CUTOFF), xlabel="rank")
        plt.gca().set(ylabel="cashflow")

        plt.subplot(3,1,3)
        plt.pcolormesh(answers.transpose(), cmap=mpl.colormaps["Greys"])
        plt.gca().plot(T, N*np.ones(CUTOFF), c='r', linewidth=2)
        plt.gca().fill_between(np.linspace(0, CUTOFF, CUTOFF), answers.sum(axis=1), alpha=0.36)
        plt.gca().set(xlim=(0,CUTOFF), xlabel="rank")
        plt.gca().set(ylabel="count")

        plt.savefig("{0:s}@{1:s}.png".format(prefix1, suffix))

        fig = plt.subplots(len(L), len(J), figsize=(8, 9))
        for lIdx in range(len(L)):
            for jIdx in range(len(J)):
                plt.subplot(len(L), len(J), 1 + len(J)*lIdx + jIdx)
                plt.title("{0:s},{1:s}".format(J[jIdx], L[lIdx]))
                plt.plot(T, characteristic[lIdx, jIdx])
                plt.gca().fill_between(T, KB[jIdx][lIdx][K_HIGH], KB[jIdx][lIdx][K_LOW], alpha=.36)
                plt.gca().set(xlim=(0,CUTOFF), xlabel="rank")
                plt.gca().set(ylabel="value")

        plt.savefig("{0:s}@{1:s}.png".format(prefix2, suffix))

if __name__ == '__main__':
    # Define sets
    L = ["I_Capital_Goods", "I_Insurance", "I_Transportation"]
    J = ["fund_enriched.pmv"]

    # Extract JSON
    _c_json = None
    indexmap = {}
    with open('./implementation/bond_data1.json', 'r') as f:
        import json
        _c_json = json.loads(f.read())
        del json

    C = []
    KC = [[],[],[],[]]
    count = 0
    for _series in _c_json["bonds"]:
        rating = _series["sector"]
        c = (_series['id'], _series["min"], _series["max"], _series["max"], _series["pmv"], _series["dxs"], 1, rating, None)
        
        C.append(c)
        lIdx = 4
        if rating == 'Capital_Goods':
            lIdx = 0
        elif rating == 'Insurance':
            lIdx = 1
        elif rating == 'Transportation':
            lIdx = 2
        else:
            lIdx = 3
        KC[lIdx].append(c)
        count += 1

    for lIdx in range(len(L)):
        for c in KC[lIdx]:
            indexmap.update({c[CID]: len(indexmap)})

    KB_pmv = [[4.758736462913, 5.803304406912, 6.847872350911, 0, 1000],
            [4.134392251914, 5.178960195913, 6.223528139912, 0, 1000], 
            [0.240014083139, 1.284582027138, 2.329149971137, 0, 1000]] # Each must be size L
    KB = [KB_pmv]

    N = 10
    R = [2.600974180557, 10.435233760548]

    for lIdx in range(len(L)):
        for c in KC[lIdx]:
            print("{0:s}:\t{1:8.3f}\t{2:8.3f} (=)\t{4:8.3f}\t{5:8.3f}\t{6:s}\t".format(c[CID], c[MINTRADE], c[MAXTRADE], c[INVENTORY], c[PMV], c[DXS], c[RATING]), c[DELTA], c[BETAS])
    indexmap

    # N range
    for N in range(0, 32):
        soo = SimplifiedOneOpto(C, KC, indexmap, L, J, KB, N, R)
        dfs = soo.run(65536)
        soo.display_results(1024, suffix="{0:d},{1:s}".format(N, str(R)))

    # K range
    for D in range(-9, 41):
        KB_pmv = [[4.758736462913, 5.803304406912, 6.847872350911, 0, 1000],
            [4.134392251914, 5.178960195913, 6.223528139912, 0, 1000], 
            [0.240014083139, 1.284582027138, 2.329149971137, 0, 1000]] # Each must be size L
        KB = [KB_pmv]
        for jIdx in range(len(J)):
            for lIdx in range(len(L)):
                radius = (KB[jIdx][lIdx][K_HIGH] - KB[jIdx][lIdx][K_LOW])/2
                KB[jIdx][lIdx][K_HIGH] += D*radius/10
                KB[jIdx][lIdx][K_LOW] -= D*radius/10
        try:
            soo = SimplifiedOneOpto(C, KC, indexmap, L, J, KB, N, R)
            dfs = soo.run(65536)
            soo.display_results(1024, suffix="{0:d}".format(D))
        except ValueError as e:
            print(D, e)
            continue