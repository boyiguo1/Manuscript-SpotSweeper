# Original function with scaling applied twice
outliers <- function(x, s = 1.4826) {
  e <- (length(x) - 1) / sqrt(length(x)) 
  mad <- function (x, center = stats::median(x), constant = s,
                   low = FALSE, high = FALSE) {
    n <- length(x)
    constant * if ((low || high) && n %% 2 == 0) {
      if (low && high) 
        stop("'low' and 'high' cannot be both TRUE")
      n2 <- n %/% 2 + as.integer(high)
      sort(abs(x - center), partial = n2)[n2]
    } else stats::median(abs(x - center))
  }
  return(((0.6745 * (x - stats::median(x))) / mad(x)))
}

# Correct version of outliers function (for comparison)
correct_outliers <- function(x, s = 1.4826) {
  e <- (length(x) - 1) / sqrt(length(x)) 
  mad <- function (x, center = stats::median(x), constant = s,
                   low = FALSE, high = FALSE) {
    n <- length(x)
    constant * if ((low || high) && n %% 2 == 0) {
      if (low && high) 
        stop("'low' and 'high' cannot be both TRUE")
      n2 <- n %/% 2 + as.integer(high)
      sort(abs(x - center), partial = n2)[n2]
    } else stats::median(abs(x - center))
  }
  return(((x - stats::median(x)) / mad(x)))
}

# Test data
set.seed(42)
test_data <- rnorm(100)  # Generate random test data

# Compute results
incorrect_results <- outliers(test_data)
correct_results <- correct_outliers(test_data)

# Test whether dividing by 0.6745 fixes the scaling issue
divide_fixed <- incorrect_results / 0.6745
divide_check <- all.equal(divide_fixed, correct_results)
divide_check
# [1] TRUE

# Test whether multiplying by 0.6745 fixes the scaling issue
multiply_fixed <- incorrect_results * 0.6745
multiply_check <- all.equal(multiply_fixed, correct_results)
multiply_check
# [1] "Mean relative difference: 1.198043"
