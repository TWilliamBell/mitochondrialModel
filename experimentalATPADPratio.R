## 1130 is the mean mmol ATP
## 1650 is the mean mmol ATP+ADP
## 175 is the standard deviation of ATP mmol
## Assuming conservation between the two quantities
## Assuming the concentration follows a truncated Student's t
## distribution with 3 degrees of freedom in order to give a large
## range for the possible ATP concentraation.
empty <- rep(0, 1e6)

for (i in 1:1e6) {
  nums <- (rt(4, 3)*(175/2)+1130)
  for (j in 1:4) {
    while (nums[j] >= 1605 | nums[j] <= 0) {
      nums[j] <- (rt(1, 3)*(175/2)+1130)
    }
  }
  
  nums <- nums/(1605-nums)
  
  empty[i] <- mean(nums)
}

summary(empty)
