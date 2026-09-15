# Shared predicate for deciding whether a matrix of expression values is already
# on a log scale and must therefore NOT be rescaled to CPM.
#
# Every normalization in this workflow is put on a common CPM scale before
# plotting or correlating, so that what is compared between methods is the part
# of the normalization that is not sequencing depth. The one exception is
# ratio_correction, which returns log2 ratios: CPM-rescaling those is
# meaningless (and divides by a column sum that can be near zero or negative).
#
# The test is value based rather than name based so that every consumer agrees
# without having to thread the normalization method name through the scripts.
# Caveat: a log-scale matrix whose values all happen to be positive would be
# treated as linear. If that ever becomes reachable, pass the method name in
# from normalization.R instead.
is_log_scale <- function(x) {
  values <- suppressWarnings(as.numeric(as.matrix(x)))
  any(is.finite(values) & values < 0)
}
