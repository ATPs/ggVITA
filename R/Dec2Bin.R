#real algerbra: 10->2 
Dec2Bin <- function(x) {
  i <- 0
  string <- numeric(32)
  while(x > 0) {
    string[32 - i] <- x %% 2
    x <- x %/% 2
    i <- i + 1 
  }
  first <- match(1, string)
  paste(as.character(string[first:32]),collapse="")
}