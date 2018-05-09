library(parallel)
if (!require(BiasedUrn, quietly = T)) {
  install.packages("BiasedUrn", quiet = T, repos='http://cran.us.r-project.org')
  library(BiasedUrn, quietly = T)
}
if (!require(argparser, quietly = T)) {
  install.packages("argparser", quiet = T, repos='http://cran.us.r-project.org')
  library(argparser, quietly = T)
}

get.p = function(gene, genome, genome.table.pa, num.genes, num.genomes, genome.size) {
  abund = genome.table.pa[gene, genome]
  if (abund==0) {
    return(1)
  }
  gene.abund = sum(genome.table.pa[gene,])
  if (num.genomes==gene.abund) {
    return(1)
  }
  gene.odds = gene.abund/num.genomes
  not.gene.odds = (num.genomes-gene.abund)/num.genomes
  return(dWNCHypergeo(abund, 1, num.genes-1, genome.size, gene.odds/not.gene.odds))
}

get.ko.p.table = function(genome.table.pa, module.genomes, cores) {
  num.genes = nrow(genome.table.pa)
  num.genomes = ncol(genome.table.pa)
  
  genome.table.p = sapply(module.genomes, function (genome) {
    print(genome)
    genome.size = sum(genome.table.pa[, genome])
    return(sapply(row.names(genome.table.pa),
                  function(gene) get.p(gene, genome, genome.table.pa,
                                       num.genes, num.genomes, genome.size)))
  })
  return(genome.table.p)
}

get.ko.p.table.3 = function(genome.table.pa, module.genomes, cores) {
  num.genes = nrow(genome.table.pa)
  num.genomes = ncol(genome.table.pa)
  
  pre.genome.table.p = mclapply(module.genomes, function (genome) {
    genome.size = sum(genome.table.pa[, genome])
    return(sapply(row.names(genome.table.pa), function(gene) get.p(gene, genome, genome.table.pa, num.genes, num.genomes, genome.size)))
  }, mc.cores = cores)
  genome.table.p <- data.frame(matrix(unlist(pre.genome.table.p), nrow=length(pre.genome.table.p[[1]]), byrow=T))
  row.names(genome.table.p) = row.names(genome.table.pa)
  colnames(genome.table.p) = module.genomes
  return(genome.table.p)
}

main = function(genome.list.loc, genome.table.loc, output.loc, cores=3) {
  genome.table.pa = read.csv(genome.table.loc, row.names = 1, check.names = F)
  module.genomes = readLines(genome.list.loc)
  genome.table.p = get.ko.p.table(genome.table.pa, module.genomes, cores)
  write.csv(genome.table.p, output.loc)
}

if (!interactive()) {
  # define argparser
  p <- arg_parser("Determine probabilities of KO genome pairs from a given list of genomes")
  p <- add_argument(p, "--genome_table", help="Table of presence abscense for a set of organisms and genes")
  p <- add_argument(p, "--genome_list", help="List of organisms which is a subset of organisms in the genome table")
  p <- add_argument(p, "--num_threads", help="Number of threads to use for multiprocessing", type="numeric", default=3)
  p <- add_argument(p, "--output", help="File name for output table of probabilities")
  
  # parse args
  argv <- parse_args(p)
  genome.table.loc <- argv$genome_table
  genome.list.loc <- argv$genome_list
  output.loc <- argv$output
  cores <- argv$num_threads
  main(genome.list.loc, genome.table.loc, output.loc, cores)
}
