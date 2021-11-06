import numpy as np
import math as math

def main():
    ## Magic numbers taken from Anita's nephron model
    href = 1e-5
    Cref = 1e-3
    diam_male = 0.0020
    diam_female = 0.0017
    pi = math.pi

    convert_male = href*Cref/10.0 #*diam_male*pi/10.0
    convert_female = href*Cref/10.0 #*diam_female*pi/10.0
    ## Converts to mmol/(s*mm^2)    ## mmol/(s*mm)

    convert_male = convert_male/1e-3
    convert_female = convert_female/1e-3

    male4 = list()
    male5 = male4
    female4 = male4
    female5 = male4

    ending = ["sup.txt",
              "jux1.txt",
              "jux2.txt",
              "jux3.txt",
              "jux4.txt",
              "jux5.txt"]

    for i in range(len(ending)):
        male4.append(np.genfromtxt(
            "../results/mTAL_ATPc/male_rat_mTAL_NaKATPase_Na14_"+ending[i]
        ))
        male5.append(np.genfromtxt(
            "../results/mTAL_ATPc/male_rat_mTAL_NaKATPase_Na15_"+ending[i]
        ))
        female4.append(np.genfromtxt(
            "../results/mTAL_ATPc/female_rat_mTAL_NaKATPase_Na14_"+ending[i]
        ))
        female5.append(np.genfromtxt(
            "../results/mTAL_ATPc/female_rat_mTAL_NaKATPase_Na15_" + ending[i]
        ))

    male = list()
    female = list()

    for i in range(len(ending)):
        male.append((male4[i]+male5[i])*convert_male)
        female.append((female4[i]+female5[i])*convert_female)
    print(male)
    print(np.sum(male[0]))
    return male
    maximum = -1.0
    minimum = 100000.0

    all = male+female
    mean = np.array([])

    for i in range(len(all)):
        oldMaximum = maximum
        oldMinimum = minimum
        maximum = np.max(np.concatenate((all[i], [maximum])))
        minimum = np.min(np.concatenate((all[i], [minimum])))
        mean = np.append(mean, np.mean(all[i]))
        if not oldMaximum == maximum:
            j = i+1
        if not oldMinimum == minimum:
            k = i+1
    print(j)
    print(k)
    print("Maximum ATP consumption:")
    print((maximum/3.0))
    print("Minimum ATP consumption:")
    print((minimum/3.0))
    print("Mean in each simulated case:")
    print(mean/3.0)
    print("Mean of all:")
    print(np.mean(mean)/3.0)

main()