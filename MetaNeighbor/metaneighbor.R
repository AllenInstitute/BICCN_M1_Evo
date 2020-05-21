
library(tidyverse)

compute_best_hits = function(dataset, labels, one_vs_all=TRUE) {
  normalized_data = normalize_cols(assay(dataset))
  colnames(normalized_data) = paste(dataset$study_id, labels, sep = "|")
  voter_id = design_matrix(colnames(normalized_data))
  voters = normalized_data %*% voter_id
  result = c()
  for (study in unique(dataset$study_id)) {
    candidates = normalized_data[, dataset$study_id == study]
    votes = crossprod(candidates, voters)
    aurocs = compute_aurocs(votes, design_matrix(rownames(votes)))
    if (one_vs_all) {
      result = rbind(result, aurocs)
    } else {
      result = rbind(result, compute_1v1_aurocs(votes, aurocs))
    }
  }
  result = result[, rownames(result)]
  return(result)
}

normalize_cols = function(M, ranked = TRUE) {
  M = as.matrix(M)
  if (ranked) {
    M = matrixStats::colRanks(M, ties.method = "average", preserveShape = TRUE)
  }
  return(scale_cols(M))
}

scale_cols = function(M) {
    apply(M, 2, function(x) {
        result = x - mean(x)
        result = result / sqrt(sum(result**2))
    })
}

compute_1v1_aurocs = function(votes, aurocs) {
  result = 0*aurocs
  for (i in seq_len(ncol(aurocs))) {
    top_candidates = find_top_candidate(votes[,i], aurocs[,i])
    result[top_candidates$best, i] <- top_candidates$score
    result[top_candidates$second, i] <- 1-top_candidates$score
  }
  return(result)
}

find_top_candidate = function(votes, aurocs) {
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
  return(list(score = score, best = best, second = second_best))
}

extract_top_candidates = function(aurocs, n = 10) {
  return(names(head(sort(aurocs, decreasing=TRUE), n = n)))
}

call_my_metaneighbor = function(dataset, gene_sets, label_matrix) {
  dat = assay(dataset)
    
  result = list()
  for(l in seq_along(gene_sets)){
    print(names(gene_sets)[l])
    subdat = dat[rownames(dat) %in% gene_sets[[l]],]
    nonzero_cells = colSums(subdat) > 0
  
    if (sum(nonzero_cells) > 0) {
        result[[names(gene_sets)[l]]] = average_auroc(
            subdat[, nonzero_cells],
            dataset$study_id[nonzero_cells],
            label_matrix[nonzero_cells,, drop = FALSE]
        )
    }
  }
  result = bind_rows(result, .id = "gene_set")
  return(result)
}

average_auroc = function(dat_sub, study_id, celltype_labels) {
  dat_sub = as.matrix(dat_sub)
  normalized_data = normalize_cols(dat_sub)
  aurocs = c()
  for (study in unique(study_id)) {
    is_study = study_id == study
    votes = crossprod(
      normalized_data[, is_study, drop = FALSE],
      normalized_data[, !is_study, drop = FALSE] %*% celltype_labels[!is_study,, drop = FALSE]
    )
    all_aurocs = compute_aurocs(votes, celltype_labels[is_study,, drop = FALSE])
    aurocs = cbind(aurocs, diag(all_aurocs))
  }
  colnames(aurocs) = unique(study_id)
  result = t(aurocs) %>%
    as_tibble(rownames = "test_dataset") %>%
    gather("cell_type", "auroc", na.rm = TRUE, -test_dataset)
  return(result)
}

compute_best_hits_from_label_sets = function(dataset, label_sets, one_vs_all = TRUE) {
  normalized_data = normalize_cols(assay(dataset))
  label_set_names = colnames(label_sets)
  
  result = c()
  for (voter_set in label_set_names) {
      voter_id = label_sets[, voter_set]
      is_not_na = !is.na(voter_id)
      voter_id = paste(dataset$study_id, voter_set, voter_id, sep = "|")
      voters = normalized_data[, is_not_na] %*% design_matrix(voter_id[is_not_na])
      
      subresult = c()
      for (candidate_set in label_set_names) {
          print(candidate_set)
          candidate_id = label_sets[, candidate_set]
          is_not_na = !is.na(candidate_id)
          candidate_id = paste(dataset$study_id, candidate_set, candidate_id, sep = "|")
          colnames(normalized_data) = candidate_id
          
          for (study in unique(dataset$study_id)) {
            if (sum(dataset$study_id == study & is_not_na) == 0) next
            candidates = normalized_data[, dataset$study_id == study & is_not_na]
            votes = crossprod(candidates, voters)
            aurocs = compute_aurocs(votes, design_matrix(rownames(votes)))
            if (one_vs_all) {
              subresult = rbind(subresult, aurocs)
            } else {
              subresult = rbind(subresult, compute_1v1_aurocs(votes, aurocs))
            }
          }
      }
      result = cbind(result, subresult)
  }
  return(result)
}

compute_cell_scores = function(dataset, labels, reference_study, target_study, scale = FALSE) {
  normalized_data = normalize_cols(assay(dataset))
  colnames(normalized_data) = paste(dataset$study_id, labels, sep = "|")
  voters = normalized_data[, dataset$study_id == reference_study]
  voters = voters %*% design_matrix(colnames(voters), scale = scale)
  candidates = normalized_data[, dataset$study_id == target_study]
  return(crossprod(candidates, voters))
}

compute_aurocs = function(predictors, label_matrix, return_tie_correction = FALSE) {
  n_positives = colSums(label_matrix)
  n_negatives = nrow(label_matrix) - n_positives
  ranks = matrixStats::colRanks(predictors, ties.method = "average", preserveShape=TRUE)
  sum_of_positive_ranks = crossprod(label_matrix, ranks)
  if (return_tie_correction) {
    tie_correction = compute_tie_correction(ranks)
  }
  colnames(sum_of_positive_ranks) = colnames(predictors)
  result = (sum_of_positive_ranks / n_positives - (n_positives+1)/2) / n_negatives
  if (return_tie_correction) {
    return(list(aurocs = result, tie_corrections = tie_correction))
  } else {
    return(result)
  }
}

compute_tie_correction = function(ranks) {
  ranks = as.matrix(ranks)
  observed_var = matrixStats::colVars(ranks)
  max_var = var(seq_len(nrow(ranks)))
  return((max_var-observed_var) * 12 / nrow(ranks))
}