benchCoexp <- function(expPath, M, f, k, pref = "identity", prefret = FALSE) {
    
    # set functions to use
    f <- get(f) # association measure function
    pref <- get(pref) # preprocessing function

    # intialize a vector for timings
    timings <- c(totalSecs = 0, loadexpSecs = 0, prefSecs = 0, fSecs = 0, rankSecs = 0, readrankSecs = 0, mrSecs = 0, benchSecs = 0)
    
    # load expression data
    cat("\tLoading expression data\n"); flush.console()
    start.time <- Sys.time()
    E <- loadRData(expPath)
    end.time <- Sys.time()
    timings["loadexpSecs"] <- difftime(end.time, start.time, units = "secs")
    n <- nrow(E)

    # preprocess expression for association measure
    cat("\tPreprocessing expression data\n"); flush.console()
    start.time <- Sys.time()
    if(prefret) E <- pref(E) else pref(E)
    end.time <- Sys.time()
    timings["prefSecs"] <- difftime(end.time, start.time, units = "secs")

    # co-association matrix calculation
    cat("\tComputing associations\n"); flush.console()
    start.time <- Sys.time()
    A <- f(E)
    end.time <- Sys.time()
    timings["fSecs"] <- difftime(end.time, start.time, units = "secs")
    E <- NULL; gc()

    # compute rank matrix
    cat("\tComputing ranks matrix\n"); flush.console()
    start.time <- Sys.time()
    fileR <- getRankMatrix(A, n, ifelse(n > 5000, 5000, n))
    end.time <- Sys.time()
    timings["rankSecs"] <- difftime(end.time, start.time, units = "secs")
    A <- NULL; gc()

    # read ranks matrix
    cat("\tReading ranks matrix\n"); flush.console()
    start.time <- Sys.time()
    R <- bigmemory::read.big.matrix(fileR, sep = '\t', type = "integer") # 1.2 mins
    end.time <- Sys.time()
    timings["readrankSecs"] <- difftime(end.time, start.time, units = "secs")

    # compute mutual rank
    cat("\tComputing mutual rank\n"); flush.console()
    start.time <- Sys.time()
    mr <- mutualRank(R@address)
    end.time <- Sys.time()
    timings["mrSecs"] <- difftime(end.time, start.time, units = "secs")
    R <- NULL; unlink(fileR); gc()

    # compute confusion matrices
    cat("\tCalculating confusion matrices\n"); flush.console()
    opps <- seq(n * 1e+04, 0, -1e+04)
    start.time <- Sys.time()
    cms <- benchShort(mr, M, opps, k)
    end.time <- Sys.time()
    timings["benchSecs"] <- difftime(end.time, start.time, units = "secs")
    mr <- M <- NULL; gc()

    # organize output
    cat("\tDONE\n"); flush.console()
    cms <- lapply(1:length(opps), function(i) cms[(i + i - 1):(i + i), ])
    cms <- simplify2array(cms)
    dimnames(cms) <- list(c('pred true', 'pred false'), c('lab true', 'lab false'), paste0('MR<=', opps / 1e+04))

    timings["totalSecs"] <- sum(timings)
    return(list(cms = cms, timings = timings))
}