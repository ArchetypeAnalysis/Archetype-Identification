#########################
#
# Compute optimal archetype suite(s) for given size and minimal richness
# for given binary data of cases and their attributes
#
# Klaus Eisenack, August 23, last update 30.4.25
#
# Main functions:
#  - atFCA(dataframe, verbose)
#  - atsel(ConceptSet, Nrange, rhorange, verbose, strict)
# Further functions:
#  - print_suites(candidates, conceptl)
#  - plot_pfrontier(nutable, filename)
#  - plot_configs(data, attrivec1, attrivec2) [experimental]

library(fcaR)

#########################
#
# Function to do the FCA and some standard analysis of the concepts
#
# atFCA(data, verbose)
# Arguments:
#   data (dataframe): contains the case/attribute data
#     - one row per case, rownames = casenames
#     - one column per attribute, colnames = attribute names
#     All attributes are considered in the FCA,
#     so it might make sense not to have outcome attributes here.
#   verbose (boolean): TRUE for more screen output, FALSE otherwise (default)
# Returns:
#   named list
#     $conceptlist: set of all concepts (ConceptSet)
#     $richness: list of concept intent sizes (vector)
#     $size: list of concept extent sizes (vector)
#     $n_configs: number of distinct configs in data
# usage makes only sense, if the FormalContext set is not needed later
# (currently, it is not re-used in this script)

atFCA <- function(data, verbose=FALSE) {
  n_attribs <- ncol(data)
  n_cases <- nrow(data)
  n_configs <- n_cases - sum(duplicated(data))

  # Do the FCA to determine the concepts
  context <- FormalContext$new(data)
  context$find_concepts()
  conceptl <- context$concepts
  n_concepts <- conceptl$size()

  size <- Matrix::colSums(conceptlist$extents()) 
  richness <- Matrix::colSums(conceptlist$intents()) 
  
  writeLines(paste('n_attribs=',  n_attribs, '  n_cases=', n_cases, 
              '  n_configs=', n_configs, '  n_concepts=', n_concepts))
  
  if (verbose) {
    context$plot()
    writeLines('Summary size:')
    print(summary(size))
    writeLines('Summary richness:')
    print(summary(richness))
    plot(richness, size)
  }
  
  return(list('context' = context,
              'conceptlist' = conceptl,
              'size' = size, 
              'richness' = richness,
              'n_configs' = n_configs))
}

#########################
#
# Function to compute optimal suites
#
# atsel(conceptl, Nrange, rhorange, verbose, strict)
#
# Function arguments:
#   conceptl (ConceptSet): contains all concepts to be considered for selection
#     (the names of the cases and the attributes are derived from conceptl)
#   Nrange (range of int, n:m): suites sizes to be considered for optimization
#   rhorange (range of int, k:l): richness to be considered for optimization
#     if k=l and n=m ("inspection mode"), joint size nu is maximized only for the specific size and richness,
#       but then, function returns *all* optima with their specific suite and configurations,
#     if otherwise ("selection mode"), the function returns joint size nu for all sizes and richness, but 
#       not the suites with their configurations
#   verbose (boolean): FALSE (default) only essential output, TRUE for screen output
#   strict (boolean): TRUE only accept concepts with exactly prescribed richness,
#                     FALSE accept concepts with prescribed or larger richness
#
# Returns a list
#   in selection mode:  $comparison (df: optimal joint size for each N, rho)
#   in inspection mode: $coverage (int: optimal joint size of all optima)
#                       $optima (matrix: all optima. One row per optimum. 
#                         Columns: concepts in suite, referenced by index in conceptl
#                       If no solution exists, returns coverage=0, optima=0.
#                       Is configured to determine not more than 100 optima. If more optima
#                           then only 100 optima a returned.
#
# The mathematics of the algorithm is explained in the appendix of:
#   Harmáčková, Eisenack, Yoshida, Sitas, O’Farrell (2025) Value Archetypes in
#   Future Scenarios: The Role of Scenario Co-Designers. Ecology and Society,
#   accepted.

library(lpSolve)

atsel <- function(conceptl, Nrange, rhorange, verbose = FALSE, strict = TRUE) {

  # clarify whether selection or inspection mode
  if (length(Nrange) == 1 && length(rhorange) == 1) 
    all_optima <- TRUE 
  else all_optima <- FALSE

  n_concepts <- conceptl$size()
  attnames <- conceptl[1]$sub(1)$get_intent()$get_attributes()
  casenames <- conceptl[1]$sub(1)$get_extent()$get_attributes()

  ####################
  # Reformat FCA results to be usable to set-up the integer program (IP)
  
  # Matrix A contains FCA configurations
  # (each row one concept, one column per attribute)

  A <- t(as.matrix(conceptl$intents()))
  rownames(A) <- c(1:n_concepts)
  colnames(A) <- attnames
  richness <- rowSums(A) # richness (cardinality of intents)

  # Matrix O contains FCA cases 
  # (each row one concept, one column per case)
  
  O <- t(as.matrix(conceptl$extents()))
  rownames(O) <- c(1:n_concepts)
  colnames(O) <- casenames

  #############
  # functions needed later to construct matrices which represent
  # the optimal suite selection problem as an integer program (IP)
  
  # generate the constraint row j with j \in 1..n_cases:
  #   nu[j] <= \sum_{i=1,..,nrow(Of)} O[i,j]*X[i]
  gen_part1 <- function(j) {
    c(Of[,j], rep(0, times=max(j-1,0)),
      -1,rep(0, times=n_cases-j)) %>% 
      as.integer()
  }
  # generate the constraint row (i,j) with i \in 1..nrow(Of), j \in 1..n_cases:
  #   nu[j] >= O[i,j]*X[i]
  gen_part2 <- function(i,j) {
    c(rep(0, times=max(i-1,0)), Of[i,j], rep(0,times=nrow(Of)-i),
      rep(0, times=max(j-1,0)),    -1,  rep(0, times=n_cases-j))
  }
  
  ###############
  # loop over N, rho (if inspection mode, loop over just one case)
  
  # table for comparison over multiple rho, N
  if (!all_optima) compare <- data.frame(rho=numeric(), N=numeric(), nu=numeric())
  
  if (verbose) start <- Sys.time()
  for (rho in rhorange) { for (N in Nrange) {
    writeLines(paste('Computing (rho, N)=', rho, N)) 
    
    ############
    # determine suite with (rho, N) which maximizes joint size

    # First setting matrices/vectors (Of, X, incl. f.obj, f.con, f.dir, f.rhs)
    # which define an equivalent integer program (IP)
    
    # Filter out all concepts with unfitting richness
    if (strict) Of <- O[richness == rho, , drop = FALSE]
      else Of <- O[richness >= rho, , drop = FALSE]

    # deal with the case of N larger than available concepts
    skip <- FALSE
    if (nrow(Of) < N) { 
      if (verbose) writeLines('Skip (rho, N) because not N concepts available')
      if (all_optima) {
        writeLines('Error: rho too large, not enough concepts with this richness')
        return(list('joint_size' = 0, 'optima' = 0))
      }
      else {
        compare[nrow(compare)+1,] <- c(rho, N, 0)
        skip <- TRUE # next rho, N
      }
    }
    if (skip) next # next rho, N

    # Vector for concept selection
    X <-  setNames(numeric(nrow(Of)), rownames(Of))

    # columns in the following matrices/vectors:
    # X[i]: column 1:nrow(Of),
    # nu[j]: column (nrow(Of)+1):(nrow(Of)+ON)
    
    # objective = \sum_{j \in 1..ON} nu[j]:
    f.obj <- c(rep(0, times=nrow(Of)), rep(1, times=n_cases)) %>% as.integer()

    # constraints 
    f.con <- matrix(nrow=(1+nrow(Of))*n_cases+2,
                     ncol=nrow(Of)+n_cases, byrow=T)
    # start building constraint matrix with rows 1:n_cases with nu[j] \le \sum...
    f.con[1:n_cases,] <- t(sapply(1:n_cases, 
                                   FUN = function(row) {gen_part1(row)}))
    # add further (nrow(Of)*ON) rows for the nu[j] >= O[i,j]*X[i] constraints:
    rowcount <- n_cases + 1
    for (i in 1:nrow(Of)) {for (j in 1:n_cases) {
      f.con[rowcount,] <- gen_part2(i,j) %>% as.integer()
      rowcount <- rowcount+1
    }}
    # add 2 further rows for N = sum X[i]
    f.con[rowcount,] <- c(rep(1, times = nrow(Of)), rep (0, times = n_cases))  %>% as.integer()
    f.con[rowcount+1,] <- c(rep(1, times = nrow(Of)), rep (0, times = n_cases))  %>% as.integer()

    f.rhs <- c(rep(0, times=n_cases),
               rep(0, times=n_cases*nrow(Of)), N, N)  %>% as.integer()
    
    f.dir <- c(rep('>=', times=n_cases),
               rep('<=', times=n_cases*nrow(Of)), '>=', '<=')
    
    # find a single solution of the IP for selection mode
    if (!all_optima) {
      sol <- lp('max', f.obj, f.con, f.dir, f.rhs, all.bin = TRUE)
      names(sol$solution) <- c(rownames(Of), casenames)
      # solution: the first n_concepts columns indicate concepts being part of
      # the optimal selection, the remaining columns indicate the
      # cases covered by at least one concept
      if (verbose) print(sol$solution)
      # sol$status # 0 = success, 2 = no feasible solution
      # sol$objval: the joint coverage
      if (sol$status == 0) compare[nrow(compare)+1,] <- c(rho, N, sol$objval)
      else compare[nrow(compare)+1,] <- c(rho, N, 0)
    }
    ## find all solutions of the IP for inspection mode
    else {
      sola <- lp('max', f.obj, f.con, f.dir, f.rhs, all.bin = TRUE, num.bin.solns=100)
      # Objective reached: joint size
      if (verbose) writeLines(paste('Joint size=', sola$objval))
      # the first colnames are the concept no., the remaining the cases' names
    }
  }} # end loop for rho, N
  
  if (verbose) {end <- Sys.time(); 
    writeLines(paste('computation time ', end-start))}
  
  if (!all_optima) {
    return (list('comparison' = compare))
    }
  else {
    numcols <- nrow(Of)+n_cases
    numsols <- sola$num.bin.solns
    # extract info on solutions in sall
    #   rows: solutions
    #   first N columns: the concepts contained in one solution (referenced by index)
    #   remaining columns: the cases contained
    sall <- matrix(head(sola$solution, numcols*numsols), nrow=numsols, byrow=TRUE)
    colnames(sall) <- c(rownames(Of), 1:n_cases)
    if (verbose) print(sall)
    # convert sall to new matrix, each row a suite, each column a concept by index
    opt <- sall[, 1:nrow(Of), drop=F]
    optima <- matrix(nrow=numsols, ncol=Nrange)
    for (j in 1:numsols) optima[j,] <- as.numeric(names(opt[j, opt[j,]==1]))
    
    return(list('joint_size' = sola$objval,
                'optima' = optima)) 
  }
}
### end of function atsel()


#########################
#
# Function to print the intents and the extents' cardinality of a candidate suite
#
# print_suites(candidates, conceptl)
#
# Argument:
#   candidates: A vector of computed optima, direct output from atsel in inspection mode
#   conceptl: ConceptSet where the candidates come from.
#     Note: conceptl and candidates need to consistent, i.e. same conceptl used
#     for atsel()
#   e.g. after result <- atsel(co_concepts, Nrange=i:i, rhorange=r:r)
#     use print_suites(result, co_concepts)
# Returns:
#   Prints the optima to the console

print_suites <- function(candidates, conceptl) {
  for (i in 1:nrow(candidates$optima)) {
    writeLines(paste0('== Candidate suite S', i, ':'))
    for (j in candidates$optima[i,]) {
      print(conceptl[j]$sub(1)$get_intent())
      writeLines(paste0('covers ',
                   conceptl[j]$sub(1)$get_extent()$cardinal(),
                   ' cases (concept C', j, ')'))
    }
  }
}
### end of function print_suites()

#########################
# Function to plot Pareto frontier from selection mode
#
# plot_pfrontier(nutable, filename)
#
# Function arguments:
#   nutable (df), as from atsel$comparison (df: optimal joint size nu for each rho, N)
#   filename (string) of png of graph 
#
# Returns:
#   plot of Pareto frontier to the screen
#   png file if graph written to working directory

library(reshape2)
library(plot3D)

plot_pfrontier <- function(nutable, filename) {
  mat <- dcast(nutable, N ~ rho, value.var = 'nu')
  rownames(mat) <- mat$N
  mat <- mat %>% subset(select = -N) %>% as.matrix()
  
  hist3D(x = as.numeric(rownames(mat)), 
         y = as.numeric(colnames(mat)), 
         z = mat, zlim = c(0, max(mat)),
         bty = "b", phi = 30,  theta = 220,
         xlab = "N", ylab = "rho", zlab = "nu", main = "optimal joint size",
         col = "#999999", border = "black", shade = 0.5,
         ticktype = "detailed", nticks = max(nutable$rho),
         space = 0.4, d = 2)
  dev.copy(png, file = paste0(filename, '.png'))
  dev.off()
}
### end of function plot_pfrontier

#########################
# Function to plot configuration frequencies and lift as heat map
#
# plot_configs(data, attrivec1, attrivec2)
#
# Function arguments:
#   data (df) orgininal binary data, possibly further binary rows to represent 
#     configurations, concepts or archetypes
#   attrivec1, attrivec2 (vectors of characters): names of the columns
#     to be used for the plot (attrilist1 on x-axis, attrilist2 on y-axis)
#
# prints
# - frequency and lift
# plots
# - numbers represents the frequency, color the lift
#   (pink above 1.0, blue belwo 1.0)

library(gplots)

plot_configs <- function(data, attrivec1, attrivec2) {
  
  n_cases <- nrow(data)
  
  freq <- matrix(ncol=length(attrivec1), nrow=length(attrivec2))
  freq <- outer(attrivec2, attrivec1,
                FUN = Vectorize(function(a2, a1) 
                  sum(data[[a1]] * data[[a2]], na.rm = TRUE)))
  dimnames(freq) <- list(attrivec2, attrivec1)
  
  lift <- matrix(ncol=length(attrivec1), nrow=length(attrivec2))
  col_sums <- vapply(attrivec1, function(a1) sum(data[[a1]], na.rm = TRUE), numeric(1))
  row_sums <- vapply(attrivec2, function(a2) sum(data[[a2]], na.rm = TRUE), numeric(1))
  lift <- outer(attrivec2, attrivec1,
    FUN = Vectorize(function(a2, a1) {
      denom <- col_sums[a1] * row_sums[a2]
      if (denom == 0) NA_real_ 
      else n_cases * freq[a2, a1] / denom
    })
  )  
  dimnames(lift) <- list(attrivec2, attrivec1)
  
  writeLines('Frequencies:')
  print(freq)
  writeLines('Lift:')
  print(lift)

  # x11() # plot on extra window, needs to be closed with dev.off()
  heatmap.2(lift,
            dendrogram = 'none', Rowv=NA, Colv=NA,
            col = cm.colors(10), 
            scale = 'none', breaks = c(seq(0,2, length=11)),
            cellnote = round(freq, digits = 2), notecol = 'black', 
            trace= 'none',
            density.info = 'none',
            key = TRUE,
            key.xlab = '', key.title = 'lift', key.par = list(cex=1.0), 
            # labRow = c('x', 'y', 'z'),
            # margins = c(5,13),
            cexCol = 2.0, cexRow = 2.0, notecex = 1.5,
            lmat=rbind( c(0, 3), c(2,1), c(0,4) ),
            lhei=c(0.1, 3, 1.1), lwid=c(0.1,0.9)
            # lmat=rbind( c(2,1), c(3,4) ) , lhei=c(3, 1), lwid=c(1, 5)
  )
}
### end of function plot_configs

