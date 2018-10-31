# both functions adapted from https://gist.github.com/traversc/1446ebe1dcc2d84bccdca781d4f1fa2a

#' Fast calculation of ROC values
#'
#' calculate roc, giving both true positive rate and false positive rate
#' 
fast_roc <- function(probs, class) {
  class_sorted <- class[order(probs, decreasing=T)]
  TPR <- cumsum(class_sorted) / sum(class)
  FPR <- cumsum(class_sorted == 0) / sum(class == 0)
  return(list(tpr=TPR, fpr=FPR))
}

# Helpful function adapted from: https://stat.ethz.ch/pipermail/r-help/2005-September/079872.html
fast_auc <- function(probs, class) {
  x <- probs
  y <- class
  x1 = x[y==1]; n1 = length(x1); 
  x2 = x[y==0]; n2 = length(x2);
  r = rank(c(x1,x2))  
  auc = (sum(r[1:n1]) - n1*(n1+1)/2) / n1 / n2
  return(auc)
}
