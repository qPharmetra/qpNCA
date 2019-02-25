signifString = function (x, digits = 6) 
{ # from PKNCA
  toplog <- bottomlog <- rep(NA, length(x))
  bottomlog[x %in% 0] <- digits
  bottomlog[x %in% c(NA, NaN) | is.infinite(x)] <- 0
  toplog <- log10(abs(x))
  mask.exact.log <- (toplog%%1) == 0
  toplog[mask.exact.log] <- toplog[mask.exact.log] + 1
  toplog <- ceiling(toplog)
  bottomlog[is.na(bottomlog)] <- digits - toplog[is.na(bottomlog)]
  roundString(x, digits = bottomlog)
}