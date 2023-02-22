dominant <- function( x ) {
  return ( ifelse( additive(x) > 0, 1, 0 ) )
}

recessive <- function( x ) {
  return ( ifelse( additive(x) > 1, 1, 0 ) ) 
}

additive <- function( x ) {
  return(as.numeric(factor(x))-1)
}

genotype <- function( x ) {
  return(as.factor(x))
}

recode <- function( x ) {
  return(factor(x))
}

diallelic <- function(x) {
  n <- nlevels(recode(x))
  if ( n == 3 ) return(TRUE)
  else return(FALSE)
}

recode <- function(x) {
  return (as.numeric(x-1))
}
