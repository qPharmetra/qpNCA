roundString = function (x, digits = 0) 
{# from PKNCA
  if (length(digits) == 1) {
    if (digits < 0) {
      formatC(round(x, digits), format = "f", digits = 0)
    }
    else {
      formatC(round(x, digits), format = "f", digits = digits)
    }
  }
  else if (length(x) == length(digits)) {
    mapply(roundString, x, digits)
  }
  else {
    stop("digits must either be a scalar or the same length as x")
  }
}