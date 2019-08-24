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