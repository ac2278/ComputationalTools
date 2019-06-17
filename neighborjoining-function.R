# This function implements the neighbor joining algorithm for an arbitrary number of sequences.
# It takes as input a matrix representing the distance between all pair of species (the matrix 
# is symmetric with zero on the diagonal). The output is a tree in Newick format (described on 
# Joe Felsenstein's webpage at http://evolution.genetics.washington.edu/ phylip/newicktree.html).

neighbor.joining = function(distance.matrix){
  d = distance.matrix
  n = nrow(d)
  tree = list()
  # Step 1: Based on the current distance matrix calculate the matrix Q.
  while(n>3){
    Q = matrix(nrow=n, ncol=n, dimnames=dimnames(d))
    for (i in rownames(d)){
      for (j in setdiff(rownames(d), i)){
        # Based on the current distance matrix relating the n taxa, calculate the
        # corresponding Q matrix as follows:
        Q[i,j] = (n-2)*d[i,j] - sum(d[i,]) - sum(d[j,]) #, where d(i,j) is the distance between taxa i and j.
      }
    }
    
    # Step 2: Find the pair of distinct taxa i and j (i.e. with i \neq j) for which Q(i,j) has its lowest value.
    inside = rownames(which(Q == min(Q, na.rm=T), arr.ind = TRUE))[1:2]
    u = paste(inside, collapse=",")
    
    # These taxa are joined to a newly created node, which is connected to 
    # the central node. Let u denote the new node. 
    
    # Step 3: Calculate the distance from each of the taxa in the pair to this new node.
    # Taxa f and g are the paired taxa and u is the newly created node.
    
    # For each of the taxa in the pair being joined, use the following formula 
    # to calculate the distance to the new node:
    delta = structure(c(NA,NA), names=inside)
    delta[1] = d[inside[1],inside[2]]/2 + (sum(d[inside[1],])-sum(d[inside[2],]))/(2*(n-2))
    delta[2] = d[inside[1],inside[2]] - delta[1]
    delta[delta<0] = 0
    
    tree[[u]] = delta
    
    # Step 4: Calculate the distance from each of the taxa outside of this pair to the new node
    # (i.e. update the distance matrix). 
    # For each taxon not considered in the previous step, we calculate the distance to the new node as follows:
    outside = setdiff(rownames(d),inside)
    d.updated = matrix(nrow=nrow(d)-1, ncol=ncol(d)-1, dimnames=list(c(u,outside), c(u,outside)))
    for (k in outside){
      d.updated[u,k] = d.updated[k,u] = (d[inside[1],k] + d[inside[2],k] - d[inside[1],inside[2]])/2
    }
    d.updated[outside,outside] = d[outside,outside]
    diag(d.updated) = 0
    d = d.updated
    n = nrow(d)
  }
  
  b = as.matrix(setdiff(unique(as.vector(d)), 0))
  
  copy = structure(d, dimnames=list(letters[1:3], letters[1:3]))
  A = matrix(nrow=3, ncol=3, dimnames=list(letters[1:3], letters[1:3]))
  for (i in 1:length(b)){
    A[i,] = table(factor(rownames(which(copy == b[i], arr.ind = TRUE)), levels=letters[1:3]))
  }
  anchor_length = structure(round(solve(A, b),5), dimnames=list(rownames(d), NULL))
  anchor_length[anchor_length<0] = 0
  
  newick = sapply(1:length(tree), function(x){ names(tree[[x]]) %in% names(tree) })
  subtree_label = names(tree[apply(newick, 2, function(x) ifelse(length(grep(F, x))>1, TRUE, FALSE))])
  
  a = lapply(tree, function(x){ paste(names(x),":", x, collapse=",", sep="") })
  for (j in 1:length(subtree_label)){
    subtree_descendant = grep(subtree_label[j], names(tree))
    for (i in 1:length(subtree_descendant)){
      a[subtree_descendant[i+1]] = sub(names(a[subtree_descendant[i]]), paste("(",a[subtree_descendant[i]],")",sep=""), a[subtree_descendant[i+1]])
    }
  }
  
  x = paste("(", a[names(a) %in% rownames(anchor_length)], "):", anchor_length[rownames(anchor_length) %in% names(a)], collapse=",", sep="")
  y = paste(rownames(anchor_length)[grep(F,rownames(anchor_length) %in% names(a))], anchor_length[grep(F,rownames(anchor_length) %in% names(a))], sep=":")
  xy = ifelse(length(y)==0, paste(x,y,collapse=",",sep=""), paste(x,y,collapse=",",sep=","))
  newick = paste("(", xy, "):0.0;", sep="")
  return(newick)
}



