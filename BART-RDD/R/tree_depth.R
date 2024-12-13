nleaf <- function(tree)
{
  if (!is.list(tree$left)) return(1)
  nleaf(tree$left) + nleaf(tree$right)
}