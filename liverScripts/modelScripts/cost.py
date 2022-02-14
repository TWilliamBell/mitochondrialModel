def cost(S2, S3):
    return ((S2-0.3)/0.3)**2+((S3-2.0)/2.0)**2

# Parameters: protein density, k_O2, XHle

## 0.5 kO2, 1.5 XHle, 10 prot density
print(cost(0.32, 2.2))

## 0.5, 1.5, 8
print(cost(0.29, 2.1))

## 0.5, 1, 10
print(cost(0.26, 2.1))

## 0.75, 1.5, 8
print(cost(0.28, 2.1))

## 0.25, 1.5, 8
print(cost(0.29, 2.1))

## 0.5, 1.75, 6
print(cost(0.4, 2.8))

## 0.5, 1.25, 6
print(cost(0.34, 2.78))

## 1., 2.1, 8
print(cost(0.3, 2.0))

## 0.3 S2, 2.4 S3
## Parameters: XHle, kO2, J_AtCmax
## 0.25-1.75 for both of the first two, 0.01-1 for the latter

# # Leak 0.25, kO2 0.25, AtC 0.01
# print(cost(0.22, 2.1))
#
# # Leak 1.0, kO2 0.25, AtC 0.01
# print(cost(0.33, 2.1))
#
# # Leak 0.5, kO2 0.25, AtC 0.01
# print(cost(0.26, 2.1))
#
# # Leak 0.75, kO2 0.25, AtC 0.01
# print(cost(0.29, 2.1))
#
# # Leak 1.25, kO2 0.25, AtC 0.01
# print(cost(0.36, 2.1))
#
# # Leak 1.5, kO2 0.25, AtC 0.01
# print(cost(0.26, 2.1))
#
# # Leak 1.75, kO2 0.25, AtC 0.01
# print(cost(0.42, 2.1))
#
# # Leak 0.25, kO2 0.5, AtC 0.01
# print(cost(0.22, 2.0))
#
# # Leak 0.5, kO2 0.5, AtC 0.01
# print(cost(0.25, 2.0))
#
# # Leak 0.75, kO2 0.5, AtC 0.01
# print(cost(0.29, 2.0)) ## Second best due to conservatism
#
# # Leak 1, kO2 0.5, AtC 0.01
# print(cost(0.32, 2.0))
#
# # Leak 1.25, kO2 0.5, AtC 0.01
# print(cost(0.38, 2.1))
#
# # Leak 1.5, kO2 0.5, AtC 0.01
# print(cost(0.41, 2.1))
#
# # Leak 1.75, kO2 0.5, AtC 0.01
# print(cost(0.40, 2.1))
#
# # Leak 0.25, kO2 1, AtC 0.01
# print(cost(0.21, 1.9))
#
# # Leak 0.75, kO2 1, AtC 0.01
# print(cost(0.28, 1.9))
#
# # Leak 1.0, kO2 1, AtC 0.01
# print(cost(0.31, 1.9))
#
# # Leak 1.25, kO2 1, AtC 0.01
# print(cost(0.34, 1.9))
#
# # Leak 1.5, kO2 1, AtC 0.01
# print(cost(0.37, 1.9))
#
# # Leak 1.75, kO2 1, AtC 0.01
# print(cost(0.4, 1.9))
#
# # Leak 0.25, kO2 0.75, AtC 0.01
# print(cost(0.21, 2.0))
#
# # Leak 1, kO2 0.75, AtC 0.01
# print(cost(0.31, 2.0)) ## Best so far
#
# # Leak 1., kO2 0.75, AtC 0.01
# print(cost(0.35, 2.0))
#
# # Leak 1.5, kO2 0.75, AtC 0.01
# print(cost(0.38, 2.0))
#
# # Leak 1.75, kO2 0.75, AtC 0.01
# print(cost(0.41, 2.0))