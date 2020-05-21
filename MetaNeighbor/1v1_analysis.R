get_variable_genes<-function(data) {
  pheno=colData(data)
  var.genes1=vector("list")
  experiment=unique(pheno$study_id)
  j=1
  for(exp in experiment){
    dat.sub=counts(data)[,pheno$study_id==exp]
      genes.list=vector("list")
      med.dat=apply(dat.sub,1,median)
      var.dat=apply(dat.sub,1,var)
      quant.med=unique(quantile(med.dat,prob=seq(0,1,length=11),type=5))
      genes.list=vector("list",length=length(quant.med))
      for(i in 1:(length(quant.med)-1)){
        if(i==1){
          filt1=med.dat<=quant.med[i]
          var.temp=var.dat[filt1]
          quant.var=quantile(var.temp,na.rm=T)
          filt2=var.temp>quant.var[4]
          genes.list[[i]]=names(var.temp)[filt2]
        }
        else {
          filt1=med.dat<=quant.med[i]&med.dat>quant.med[i-1]
          var.temp=var.dat[filt1]
          quant.var=quantile(var.temp,na.rm=T)
          filt2=var.temp>quant.var[4]
          genes.list[[i]]=names(var.temp)[filt2]
        }
      }
      var.genes1[[j]]=unlist(genes.list)
      j=j+1
  }
  var.genes=Reduce(intersect, var.genes1)
  return(var.genes)
  }

compute_best_hits = function(dataset, labels, one_vs_one = TRUE) {
  normalized_data = normalize_cols(assay(dataset))
  colnames(normalized_data) = paste(dataset$study_id, labels, sep = "|")
  voter_id = design_matrix(colnames(normalized_data))
  voters = normalized_data %*% voter_id
  result = c()
  for (study in unique(dataset$study_id)) {
    candidates = normalized_data[, dataset$study_id == study]
    votes = crossprod(candidates, voters)
    aurocs = compute_aurocs(votes, design_matrix(rownames(votes)))
    if (one_vs_one) {
        result = rbind(result, compute_1v1_aurocs(votes, aurocs))
    } else {
        result = rbind(result, aurocs)
    }
  }
  result = result[, rownames(result)]
  return(result)
}

normalize_cols <- function(M, ranked = TRUE) {
  M <- as.matrix(M)
  if (ranked) {
    M <- matrixStats::colRanks(M, ties.method = "average", preserveShape = TRUE)
  }
  return(scale_cols(M))
}

Rcpp::cppFunction('NumericMatrix scale_cols(NumericMatrix M) {
  NumericMatrix result(M.nrow(), M.ncol());
  for (int j = 0; j < M.ncol(); j++) {
    double m = 0;
    for (int i = 0; i < M.nrow(); i++) { m += M(i,j); }
    m /= M.nrow();
    for (int i = 0; i < M.nrow(); i++) { result(i,j) = M(i,j) - m; }
    double s = 0;
    for (int i = 0; i < M.nrow(); i++) { s += result(i,j) * result(i,j); }
    s = 1 / sqrt(s);
    for (int i = 0; i < M.nrow(); i++) { result(i,j) *= s; }
  }
  return result;
}')

design_matrix <- function(cell_type, scale = FALSE) {
  factors <- levels(as.factor(cell_type))
  if (length(factors) > 1) {
    result <- model.matrix(~cell_type-1)
  } else {
    result <- matrix(1, nrow = length(cell_type), ncol = 1)
  }
  colnames(result) <- factors
  if (scale) {
    result <- t(t(result) / colSums(result))
  }
  return(result)
}

compute_aurocs <- function(predictors, label_matrix) {
  n_positives <- colSums(label_matrix)
  n_negatives <- nrow(label_matrix) - n_positives
  sum_of_positive_ranks <- crossprod(
    label_matrix,
    matrixStats::colRanks(predictors, ties.method = "average", preserveShape=TRUE)
  )
  colnames(sum_of_positive_ranks) <- colnames(predictors)
  result <- (sum_of_positive_ranks / n_positives - (n_positives+1)/2) / n_negatives
  return(result)
}

compute_1v1_aurocs = function(votes, aurocs) {
  result = NA*aurocs
  for (i in seq_len(ncol(aurocs))) {
    top_candidates = find_top_candidates(votes[,i], aurocs[,i])
    result[top_candidates$best, i] = top_candidates$score
    result[top_candidates$second, i] = 1-top_candidates$score
  }
  return(result)
}

find_top_candidates = function(votes, aurocs) {
  candidates = extract_top_candidates(aurocs, 5)
  best = candidates[1]
  votes_best = votes[names(votes) == best]
  score = 1
  second_best = candidates[2]
  for (i in seq(2, length(candidates))) {
    contender = candidates[i]
    votes_contender = votes[names(votes) == contender]
    auroc = c(compute_aurocs(
      as.matrix(c(votes_best, votes_contender)),
      as.matrix(rep(c(1,0), c(length(votes_best), length(votes_contender))))
    ))
    if (auroc > 0.5) {
      if (auroc < score) {
        score = auroc
        second_best = contender
      }
    } else {
      second_best = best
      best = contender
      score = 1 - auroc
      votes_best = votes_contender
    }
  }
  return(list(score = score, best=best, second=second_best))
}

extract_top_candidates = function(aurocs, n = 10) {
  return(names(head(sort(aurocs, decreasing=TRUE), n = n)))
}
