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

def x(c):
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
        self.L = L
        self.J = J
        self.KB = KB
        self.beta = beta
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
                    matrix[indexmap[c[CID]]][indexmap[c[CID]]] -= rho[jIdx]*2*KB[jIdx][lIdx][K_TARGET]*self.beta[jIdx][i0]*c1[PMV]*x(c)

        # Consider constant term
        constantterm = np.zeros((len(J), len(L)))
        for lIdx in range(len(L)):
            for jIdx in range(len(J)):
                constantterm[jIdx][lIdx] = KB[jIdx][lIdx][K_TARGET]*KB[jIdx][lIdx][K_TARGET]

        self.matrix = matrix
        self.constantterm = constantterm

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
                    sxl.append((self.c[CID], self.c[PMV]*self.c[DELTA]*self.beta[j][i0]*x(self.KC[l][i])))
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
        df_s

def generate_from_C(df):
    """ Generate C from the dataframe """
    C = []
    indexmap = {}
    for c in df.iterrows():
        _series = c[1]
        C.append((_series['price'], _series['minTradeSize'], _series['minTradeSize'], _series['minTradeSize'], _series['minTradeIncrement'], _series['assetId'], [BETA]*len(J)))
        indexmap.update({_series['assetId']: len(indexmap)})
    for c in C:
        print("{5:s}:\t{0:6.2f}%\t{1:8.3f}\t{1:8.3f}\t{1:8.3f}\t{1:8.3f}\t".format(c[0], c[1], c[2], c[3], c[4], c[5]), [BETA]*2)
    return C, indexmap

def generate_from_sample():
    """ Generate C from the sample data """
    # Extract JSON
    _c_json = None
    indexmap = {}
    with open('./bond_data1.json', 'r') as f:
        import json
        _c_json = json.loads(f.read())
        del json

    C = []
    KC = [[],[],[],[]]
    count = 0
    for _series in _c_json["bonds"]:
        rating = _series["sector"]
        c = (_series['id'], _series["min"], _series["max"], _series["max"], _series["pmv"], _series["dxs"], 1, rating, [BETA]*len(J))
        
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

    return C, KC, indexmap

if __name__ == '__main__':
    with open('./run_aer.txt', 'r') as f:
        print(f.read())

    C, KC, indexmap = generate_from_sample()
    L = ["I_Capital_Goods", "I_Insurance", "I_Transportation"]
    J = ["fund_enriched.pmv"]
    NL, NU = -2000000000, 2000000000
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

    soo = SimplifiedOneOpto(C, KC, indexmap, L, J, KB, N, R)
    soo.run(100)
