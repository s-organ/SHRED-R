library(dplyr)
library(glmnet)
library(ClustOfVar)

####### INPUTS ##############################

# x: a design matrix
# y: response vector
# q = 0.05 (or alternative desired FDR control level)


##### HIERARCHICAL CLUSTERING OF VARIABLES #######################

hc = hclustvar(x)

memb = hc[["clusmat"]]

# memb is the matrix detailing the clusters obtained from the hierarchical clustering

p = ncol(x)

# turn memb into a matrix with 2p-1 rows of 1's and 0's for which variables are in which set
mat2 = matrix(ncol = p)
for(v in 1:(nrow(memb))){
  testing = cbind.data.frame(c(1:p),memb[,v])
  num_sets = length(unique(c(memb[,v])))
  clusters = unique(memb[,v])
  mat = matrix(nrow = num_sets, ncol = p)
  for(i in 1:num_sets){
    set = which(testing[,2] == clusters[i])
    in_set = ifelse(c(1:p) %in% set, 1,0)
    mat[i,] = c(in_set)
  }
  mat2 = rbind(mat2,mat)
  mat2 = mat2 %>% as.data.frame() %>% distinct() %>% as.matrix()
}

mat2 = mat2[-1,]

############### HYPOTHESIS TESTING ##########################

# NOTE: this is Gaussian, change accordingly for logistic/poisson

reg = lm(y ~ x)
p_val = c()
for(i in 1:(nrow(mat2))){
  cat(paste0(round(i / (nrow(mat2)) * 100), '% Testing sets completed'))
  set = mat2[i,]
  in_set = which(set == 1)
  size = length(in_set)
  if(size < p){
    x2 = x[,-c(in_set)]
    nested = lm(y ~ x2)
    p_val[i] = anova(reg,nested)[2,6]
  } else{
    nested = lm(y ~ 1)
    p_val[i] = anova(reg,nested)[2,6]
  }
  if (i == (nrow(mat2)-1)) cat(': Done')
  else cat('\014')
}
results = cbind.data.frame(mat2,p_val)
colnames(results)[p+1] = "p"

size = c()
for(i in 1:(nrow(mat2))){
  size[i] = sum(mat2[i,])
}

hypothesis = cbind.data.frame(c(1:(nrow(mat2))),p_val,size)
colnames(hypothesis)[1] = "set"

# --- Inputs ---
# mat2: binary membership matrix (rows = sets, cols = variables)
# hypothesis: data.frame with columns: set (row index in mat2), p_val, size

# 1) Order hypotheses by ascending p-value
hypothesis <- hypothesis[order(hypothesis$p_val), ]
set_vars <- apply(mat2, 1, function(row) which(row != 0))
hypothesis$weight = 1/hypothesis$size

# Initialize tracking columns
hypothesis$nested_weight <- 0
hypothesis$nest <- "no"
minimal_sets_so_far <- list()
minimal_sets_indices <- integer(0)
cum_sigma <- numeric(nrow(hypothesis))

# 2) get minimal sum of weights (sigma) at each p-value

for (i in seq_len(nrow(hypothesis))) {
  current_set_idx <- hypothesis$set[i]
  current_vars <- set_vars[[current_set_idx]]
  
  # Check if current set contains any minimal set (if yes, current set is redundant)
  contains_minimal <- FALSE
  for (ms_vars in minimal_sets_so_far) {
    if (all(ms_vars %in% current_vars)) {
      contains_minimal <- TRUE
      break
    }
  }
  
  if (contains_minimal) {
    # Current set is redundant (larger containing a smaller minimal set)
    hypothesis$nested_weight[i] <- 0
    hypothesis$nest[i] <- "yes"
  } else {
    # Current set is minimal, add it
    hypothesis$nested_weight[i] <- 1 / length(current_vars)
    hypothesis$nest[i] <- "no"
    
    # Remove minimal sets contained in current set (they are now redundant)
    to_remove <- logical(length(minimal_sets_so_far))
    for (j in seq_along(minimal_sets_so_far)) {
      if (all(current_vars %in% minimal_sets_so_far[[j]])) {
        # current set is contained in minimal_sets_so_far[j], so the old minimal set is larger and should be removed
        to_remove[j] <- TRUE
      }
    }
    
    if (any(to_remove)) {
      minimal_sets_so_far <- minimal_sets_so_far[!to_remove]
      minimal_sets_indices <- minimal_sets_indices[!to_remove]
    }
    
    # Add current minimal set
    minimal_sets_so_far[[length(minimal_sets_so_far) + 1]] <- current_vars
    minimal_sets_indices <- c(minimal_sets_indices, i)
  }
  
  # Cumulative sum of nested weights (i.e., sigma)
  cum_sigma[i] <- sum(vapply(minimal_sets_so_far, function(s) 1 / length(s), numeric(1)))
}


hypothesis$sigma = cum_sigma

# NOTE: the last sigma (final row of hypothesis)  should be p

# 3) find c_max

total_weight <- sum(hypothesis$weight)
threshold <- hypothesis$p_val * (total_weight / q)

pass_index <- which(cum_sigma >= threshold)

if (length(pass_index) > 0) {
  k <- max(pass_index)
  selected_sets <- hypothesis[1:k, ]
} else {
  selected_sets <- hypothesis[0, ]
}

# 4) Final rejection: keep only minimal selected sets

set_vars <- apply(mat2, 1, function(row) which(row != 0))
selected_indices <- selected_sets$set
selected_vars_list <- set_vars[selected_indices]

# Keep only minimal (non-redundant) sets
keep <- rep(TRUE, length(selected_vars_list))

for (i in seq_along(selected_vars_list)) {
  for (j in seq_along(selected_vars_list)) {
    if (i == j) next
    vars_i <- selected_vars_list[[i]]
    vars_j <- selected_vars_list[[j]]
    if (all(vars_j %in% vars_i) && length(vars_j) < length(vars_i)) {
      keep[i] <- FALSE
      break
    }
  }
}

final_selected <- selected_sets[keep, ]

################ IF COMPUTING FDR AND POWER #######################

# NOTE:  truth = vector of true variables 

truth_in = c()
for(i in 1:nrow(final_selected)){
  set = final_selected$set[i]
  variables = which(mat2[set,] != 0)
  truth_in[i] = length(which(variables %in% truth))
}

final_selected$truth_in  = truth_in


final_selected = final_selected %>% dplyr::mutate(power_weight = ifelse(truth_in !=0, weight, 0))

power = sum(final_selected$power_weight)/length(truth)
FDR = sum(final_selected$weight[final_selected$truth_in == 0])/sum(final_selected$weight)

####### IF COMPUTING MSE: CHOOSE A VARIABLE AT RANDOM FROM EACH SET FOR MSE ######################

variables = c()
for(i in 1:nrow(final_selected)){
  tryCatch({
    set = final_selected$set[i]
    vars = which(mat2[set,] != 0)
    variables[i] = sample(vars,1)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}



