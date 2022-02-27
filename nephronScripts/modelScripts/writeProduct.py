import itertools
import pandas as pd
import feather

def main(): ## Runs differential equation for time span and outputs results to
    ## a csv file and a feather file.
    w = [0.25, 0.5, 0.75, 1]
    a = []
    for i in itertools.product(w, w, w, w):
        a.append(list(i))
    iterProd = pd.DataFrame(a)
    #feather.write_dataframe(iterProd, "../results/iterProd.feather")
    iterProd.to_csv("../results/iterProd.csv")
    w = [0.25, 0.5, 0.75, 1]
    h = [1., 1.15, 1.5, 2, 5, 10]
    p = [1., 1.15, 2]
    k = 0
    b = list()
    for i in itertools.product(w, w, w, w):
        for j in itertools.product(h, p):
            k+=1
            b.append([k]+list(i)+list(j))
    iterProdDrugSim = pd.DataFrame(b)
    iterProdDrugSim.to_csv(
        "../results/iterProdDrugSim.csv")

    w = [0.25, 0.5, 0.75, 1.0]
    h = [1., 1.15, 5, 10]
    o2 = [0.25, 0.5, 1]
    glyc = [0.0, 1.21e-6, 9.8e-6, 1.83e-5]
    b = list()

    k = 0
    for i in range(4):
        for l in w:
            for j in itertools.product(h, glyc, o2):
                weight = [1., 1., 1., 1.]
                weight[i] = l
                k+=1
                b.append([k]+weight+list(j))
    iterProdBroadSim = pd.DataFrame(b)
    iterProdBroadSim.to_csv(
        "../results/iterProdBroadSim.csv")

    w1 = [0.67, 0.75, 1]
    w3 = [0.75, 1]
    w4 = [0.25, 0.5, 0.75, 1]
    w5 = [1]
    h = [1., 1.15, 1.5, 5, 10]
    pO2 = [0.1, 0.5, 1]
    k = 0
    b = list()
    for i in itertools.product(w1, w3, w4, w5):
        for j in itertools.product(h, pO2):
            k+=1
            b.append([k]+list(i)+list(j))
    iterProdDiabetes = pd.DataFrame(b)
    iterProdDiabetes.to_csv(
    "../results/iterProdDiabetes.csv")




main()
