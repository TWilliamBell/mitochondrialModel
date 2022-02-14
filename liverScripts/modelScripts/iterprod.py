import itertools
import pandas as pd

def main():
    w = [0.25, 0.5, 0.75, 1.0]
    o2 = [0.25, 0.5, 1]
    glyc = [0.0, 2.3e-4, 3.4e-4]

    k = 0
    iter = list()
    for i in itertools.product(w, w, w, w):
        for j in itertools.product(glyc, o2):
            k += 1
            nums = [k]+list(i)+list(j)
            iter.append(nums)
    df = pd.DataFrame(iter, columns = ['nums', 'CI', 'CIII',
                                       'CIV', 'ATP',
                                       'Glycolysis', 'PO2'])
    df.to_csv("../results/iterprod.csv")

    w1 = [0.5, 0.75, 1]
    w3 = [0.75, 1]
    w4 = [0.25, 0.5, 0.75, 1, 1.75]
    w5 = [1]
    h = [1., 1.15, 1.5, 5, 10]
    pO2 = [0.1, 0.5, 1]

    k = 0
    iter = list()
    for i in itertools.product(w1, w3, w4, w5):
        for j in itertools.product(h, pO2):
            k += 1
            nums = [k] + list(i) + list(j)
            iter.append(nums)
    df = pd.DataFrame(iter, columns=['nums', 'CI', 'CIII',
                                         'CIV', 'ATP',
                                         'Leak', 'PO2'])
    df.to_csv("../results/iterprodDiabetes.csv")

main()