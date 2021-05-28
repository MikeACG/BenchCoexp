getRankMatrix <- function(assoc, n, nchunk) {

    file.create("temp.txt")
    R <- matrix(as.integer(0), nrow = nchunk, ncol = n)
    i <- 1
    last <- FALSE
    for(g in 1:n){

        if(g == n) last = TRUE
        glist <- getIndexPairs(g - 1, assoc, n)
        R[i, ] <- rank(-abs(glist), ties.method = "min") # highest is rank 1

        if(i %% nchunk == 0 || last){ # save progress and reset matrix if necessary
            data.table::fwrite(
                data.table::data.table(R),
                file = "temp.txt",
                append = TRUE,
                sep = '\t',
                col.names = FALSE,
                quote = FALSE
            )
            R <- NULL; gc()
            if(g < n){ # still not done, allocate another chunk
                if(n - g < nchunk) remaining <- n - g else remaining <- nchunk
                R <- matrix(as.integer(0), nrow = remaining, ncol = n)
                i <- 0
            }
        }
        i <- i + 1

    }

    return("temp.txt")
}
