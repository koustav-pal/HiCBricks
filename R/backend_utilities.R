return_chromosome_pairs <- function(chromosomes, type){
    Chromosome.pairs <- data.frame(chr1 = rep(seq_along(chromosomes), 
            each = length(chromosomes)), 
            chr2 = rep(seq_along(chromosomes), 
            times = length(chromosomes)))
    if(type == "both"){
        Chromosome.pairs <- Chromosome.pairs[
            Chromosome.pairs[,'chr2'] >= Chromosome.pairs[,'chr1'],]
    }else if(type == "cis"){
        Chromosome.pairs <- Chromosome.pairs[
            Chromosome.pairs[,'chr2'] == Chromosome.pairs[,'chr1'],]
    }else if(type == "trans"){
        Chromosome.pairs <- Chromosome.pairs[
            Chromosome.pairs[,'chr2'] > Chromosome.pairs[,'chr1'],]
    }
    Chromosome.pairs[,'chr1'] <- chromosomes[Chromosome.pairs[,'chr1']]
    Chromosome.pairs[,'chr2'] <- chromosomes[Chromosome.pairs[,'chr2']]
    Chromosome.pairs.split <- split(Chromosome.pairs[,'chr2'], 
        Chromosome.pairs[,'chr1'])
    return(Chromosome.pairs.split)
}

.normalize_by_distance_values <- function(a_matrix){
    max_dist <- nrow(a_matrix)
    distance_values <- vapply(seq(from = 1, to = max_dist), function(x){
        dist_val <- (x - 1)
        Rows <- seq(from = 1, to = (max_dist - dist_val))
        Cols <- Rows + dist_val
        mean(a_matrix[cbind(Rows, Cols)])
    }, 1)
    Rows <- row(a_matrix)
    Cols <- col(a_matrix)
    dist_mat <- Cols - Rows
    dist_mat <- abs(dist_mat) + 1
    norm_mat <- .remove_nas(a_matrix/distance_values[dist_mat])
    return(norm_mat)
}

.remove_nas <- function(a_vector){
    a_vector[is.na(a_vector) | is.null(a_vector) | is.nan(a_vector)] <- 0
    return(a_vector)
}