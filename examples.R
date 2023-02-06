##### One-sample t-Test Example

######################################################
### Run the Test of Tests for given data
######################################################

# Generate a data set
n <- 30
data <- data.frame(Y = rnorm(n, mean = 1, sd = 1))

# Define the pub_test function
pub_test <- function(data){
  t.test(x = data)$p.value
}

# Run the test
test.of.tests(data, pub_test, epsilon = 1, m = 3, alpha0 = 0.2)


######################################################
### Compute Theoretical Power
######################################################

# Define pub_power function
pub_power <- function(effect, n, alpha){
  power.t.test(n = n, delta = effect, sd = 1, sig.level = alpha,
               type = "one.sample", alternative = "two.sided",
               strict = T)$power
}

# Compute power
power.test.of.tests(effect = 1, pub_power, epsilon = 1, n = 30, m = 3, alpha0 = 0.2)

# Confirm result empirically
nsims <- 1000

samps <- rep(NA, nsims)
for(i in 1:nsims){
  data <- data.frame(Y = rnorm(n, mean = 1, sd = 1))
  samps[i] <- test.of.tests(data, pub_test, epsilon = 1, m = 3, alpha0 = 0.2)$p.value
}
mean(samps <= 0.05)




######################################################
### Find optimal m and alpha0 for unknown effect size
######################################################

# Find the parameters with target power rho = 0.8
pars <- practical.m.alpha0(pub_power, epsilon = 1, n = 30, rho = 0.8, 
                           effect_grid = seq(0.1, 2, 0.1))
pars

# Find the power with these parameters for effect = 1
power.test.of.tests(effect = 1, pub_power, epsilon = 1, n = 30, m = pars$m, 
                    alpha0 = pars$alpha0)

# Compare to optimized power when the effect is known
optimal.m.alpha0(effect = 1, pub_power, epsilon = 1, n = 30)$power

