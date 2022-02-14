def cost(RCRmTAL, POmTAL):
    return ((RCRmTAL-10.0)/1.6)**2+((POmTAL-1.92)/0.13)**2

## Default
print(cost(9.9, 0.98))

## Half Leak, Norm kO2
print(cost(12.7, 1.03))

## Tenth Leak, Norm kO2
print(cost(16.72, 1.06))

## Norm Leak, Half kO2 - best so far
print("Norm Leak, Half kO2")
print(cost(10.44, 1.77))

## Half Leak, Half kO2
print(cost(13.4, 1.84))

## Tenth Leak, Half kO2
print(cost(17.61, 1.91))

## Norm Leak, Tenth kO2
print(cost(10.98, 5.1))

## Half Leak, Tenth kO2
print(cost(14.06, 5.23))

## Tenth Leak, Tenth kO2
print(cost(18.46, 5.51))

## 0.6 K, 1.1 H
print(cost(9.9, 1.51))

## 0.6 K, 1 H
print(cost(10.32, 1.52))

## 0.6 K, 0.9 H
print(cost(10.79, 1.54))

## 0.5 K, 1.1 H
print(cost(10.01, 1.75))

## 0.5 K, 0.9 H
print(cost(10.01, 1.75))

## 0.4 K, 1.1 H
print(cost(10.13, 2.09))

## 0.4 K, 1 H
print(cost(10.56, 2.11))

## 0.4 K, 0.9 H
print(cost(11.04, 2.12))

## 0.5 k_O2, 1.0 H Leak still the best case