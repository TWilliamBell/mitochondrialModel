import numpy as np
import math as math

def main():
    ## Magic numbers taken from Anita's nephron model
    href = 1e-5
    Cref = 1e-3
    diam_male = 0.0020
    diam_female = 0.0017
    pi = math.pi

    convert_male = 2.*href*Cref*diam_male*pi/10.0
    convert_female = 2.*href*Cref*diam_female*pi/10.0
    ## Converts to mmol/(s*tubule)

    male_tubulevol = (1/1e-3)*(1/2.)*pi*diam_male
    female_tubulevol = (1/1e-3)*(1/2.)*pi*diam_female

    convert_male = convert_male*male_tubulevol
    convert_female = convert_female*female_tubulevol

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
        male.append(male4[i]+male5[i])
        female.append(female4[i]+female5[i])

    male = [convert_male*np.sum(male[i]) for i in range(len(male))]
    female = [convert_female*np.sum(female[i]) for i in range(len(female))]

    print(male)
    print(female)

main()