## Wrapper script for simulator + mopsocd or other algosw

source('E:/GoForIT2/Optimizer/rvea/F_misc.R')
source('E:/GoForIT2/Optimizer/rvea/utils.R')
source('E:/GoForIT2/Optimizer/rvea/Rvea6.R')

# NOTE There must be only one RSM in parent_dir at a time
# so that script is not confused what RSM to read.

## User parameters
input_dir <- "deck/"
output_dir <- "output/"
parent_dir <- "./"
python_cmd <- "C:/Python27/python.exe converter.py" #"C:\python278\python.exe converter.py"
model_name <- "PUNQS3_SET1_HM"
simtype <- 1 # for type of simulator (1 or 2), to handle slight differences in formatting
varcnt <- 24 # number of parameters (24 in deck 1, 19 in deck 2)
nv <- 18 # number of production variables (18 in deck 1, 20 in deck 2)
weights <- rep(1,18) # wi for misfits
objfn <- 6 # number of objective functions
objf <- list(c(1,7,13), c(2,8,14), c(3,9,15), c(4,10,16), c(5,11,17), c(6,12,18)) # for sums of misfits (based on well)
npop <- 21
maxgen <- 110
p1 <- 2 # RVEA
p2 <- 0 # RVEA
algorithm <- rvea6 # mopsocd # mopsopsa
  
# NOTE Here iteration is in {1,... npop*maxgen+npop},
# i.e. it is the number for calls to simulator and all saved folders in output_dir

# current value of params
params <- NULL
# misfits values
misfits <- rep(0, nv)
# misfits against iterations for plot
stored_misfits <- matrix(0, nrow=npop*maxgen+npop, ncol=nv)
# parameters against iterations for plot
# (storedParam of mopsocd stores only npop*maxgen)
stored_param <- matrix(0, nrow=npop*maxgen+npop, ncol=varcnt)

## Variables that will be read from the files
# vector of strings, with parameter names
# populated by process_dat
paramnames <- NULL
# lower and upper bounds of parameters
# populated by process_dat
lbounds <- NULL
ubounds <- NULL
# vector of strings with production variables
# populated by process_obs and process_sig
prodnames <- NULL
# matrix with observation values from .OBS
# populated by process_obs
obs <- NULL
# matrix with sigma values from .SIG
# populated by process_sig
sig <- NULL
# matrix with sim values from .RSM with set colnames
# populated by process_rsm
rsm <- NULL

# names of files that must be deleted after successful run of algorithm
delfiles <- c()

library(stringi)
# in order to avoid troubles with RSM-type 2 reading
Sys.setlocale(locale="C")

type1 <- simtype == 1
type2 <- simtype == 2
# operators taken into account: +, *, pow(), ()
calc_formula <- function(line) {
    leq <- function(op1, op2) {
        op1 == op2 || (op1 == '+' && op2 == '*')
    }
    calc <- function(op, v1, v2) {
        if(op == '+') v1 + v2
        else if(op == '*') v1 * v2
        else if(op == 'p') v2 ^ v1
    }
    # no whitespaces
    line = stri_replace_all(line, "", fixed=" ")
    # borders between numbers
    borders <- stri_locate_all(line, regex="[+|*|p|(|)|,]")[[1]]
    i <- 1 # current position in line
    values <- c()
    operators <- c()
    b <- NULL
    while(length(borders) >= 2) {
        # the last border character in line or not
        b = ifelse(length(borders) == 2, borders[1], borders[1,1])
        if(i < b) {
            values = c(as.numeric(substr(line, i, b-1)), values)
            i = b
        } else if(i == b) {
            token <- substr(line, i, i)
            i = ifelse(token == 'p', b+4, b+1) # for "pow(" and one-symbol borders respectively
            if(token == 'p') {
                operators = c(token, operators)
                # there are 2 borders in 'pow(': p and (
                borders = borders[-1,]
            } else if(token == '(') {
                operators = c(token, operators)
            } else if(token == ')') {
                while(operators[1] != '(' && operators[1] != ',') {
                    v <- calc(operators[1], values[1], values[2])
                    operators = operators[-1]
                    values = c(v, values[-(1:2)])
                }
                if(operators[1] == ',') {
                    operators = operators[-1]
                    # now it is 'p' on top, so calculate it
                    v <- calc(operators[1], values[1], values[2])
                    values = c(v, values[-(1:2)])
                }
                # discard 'p' or '('
                operators = operators[-1]
            } else if(token == ',') {
                # calculate 1st argument for pow
                while(operators[1] != 'p') {
                    v <- calc(operators[1], values[1], values[2])
                    operators = operators[-1]
                    values = c(v, values[-(1:2)])
                }
                operators = c(token, operators)
            } else if(length(operators) > 0 && leq(token,operators[1])) {
                # can calculate + or * on top
                v <- calc(operators[1], values[1], values[2])
                operators = c(token, operators[-1])
                values = c(v, values[-(1:2)])
            } else if (token == '+' || token == '*') {
                operators = c(token, operators)
            }
            if(length(borders)==2) break;
            borders = borders[-1,]
        }
    }
    # check for last number
    if(i <= nchar(line)) {
        values = c(as.numeric(substring(line, i)), values)
    }
    while(length(operators) > 0) {
        v <- calc(operators[1], values[1], values[2])
        operators = operators[-1]
        values = c(v, values[-(1:2)])
    }
    values[1]
}

# write_formula_file('deck/PORO_PERM_SET1.PRP', './PORO_PERM_SET1.PRP')
write_formula_file <- function(f, outf){
    formula_line <- function(line) {
        line = stri_replace_all(line, 'pow', fixed='Math.pow')
        for(i in seq_along(paramnames)) {
            line = stri_replace_all(line, params[i], fixed=paramnames[i])
        }

        begin <- stri_locate_all(line, fixed='{')[[1]][,1]
        if(!any(is.na(begin))) {
            newline <- substr(line, 1, begin[1]-1)
            end <- stri_locate_all(line, fixed='}')[[1]][,1]
            for(i in seq_along(begin)) {
                expr <- substr(line, begin[i]+1, end[i]-1)
                newline = ifelse(i < length(begin),
                            paste(newline, calc_formula(expr), substr(line, end[i]+1, begin[i+1]-1), sep=""),
                            paste(newline, calc_formula(expr), substring(line, end[i]+1), sep=""))
            }
        } else newline <- line
        newline
    }
    l <- scan(f, what="", sep="\n")
    tocalc <- grepl('$', l, fixed=T)
    if(sum(tocalc) > 0) # there are parameter values
        l[tocalc] = unlist(lapply(l[tocalc], formula_line))
    fc <- file(outf)
    writeLines(l, fc)
    close(fc)
    cat("processed file ", f, "\n")
}
# process_dat(input_dir, 'Distributions_SET1.dat')
process_dat <- function(folder, f){
    f = paste(folder, f, sep="")
    file.copy(f, parent_dir, overwrite=TRUE)

    l <- scan(f, what=list("",""), sep="\t")
    paramnames <<- l[[1]]
    lbounds <<- as.numeric(lapply(sapply(l[[2]], strsplit, "[(|,|)]"), '[', 2))
    ubounds <<- as.numeric(lapply(sapply(l[[2]], strsplit, "[(|,|)]"), '[', 3))
    cat("processed file ", f, "\n")
}
# process_data(input_dir, 'PUNQS3_SET1_HM.DATA')
# NOTE params and paramnames must already be populated
process_data <- function(folder, f){
    dataf = paste(folder, f, sep="")

    l <- scan(dataf, what="", sep="\n")
    l = trimws(l[which(l=="INCLUDE")+1])
    tocalc = stri_trim(l, pattern="[:Letter:]")
    for(i in seq_along(tocalc)) {
        fname <- paste(folder,tocalc[i],sep="")
        outfname <- paste(parent_dir,tocalc[i],sep="")
        write_formula_file(fname, outfname)
    }
    # .DATA is also a formula file
    write_formula_file(dataf, paste(parent_dir,f,sep=""))
    cat("processed file ", dataf, "\n")
}
# parse numerical matrices as in .OBS, .SIG, .RSM files
parse_matrix <- function(lines) {
    v <- as.numeric(stri_split_coll(lines[1], "\t", simplify=T))
    n <- length(v)
    m <- matrix(v, ncol=1)
    for(i in 2:length(lines)) {
        v <- as.numeric(stri_split_coll(lines[i], "\t", simplify=T))
        # there is not enough fields, but there is a timestamp (1st field)
        if(length(v) < n) break
        m = cbind(m, v)
    }
    t(m)
}
# process_obs(input_dir, 'PUNQS3.OBS')
process_obs <- function(folder, f){
    f = paste(folder, f, sep="")
    file.copy(f, parent_dir, overwrite=TRUE)
    l <- scan(f, what="", sep="\n")
    prodnames <<- strsplit(l[1], '\t', fixed=TRUE)[[1]]
    i <- ifelse(type1, 3, 2)
    obs <<- parse_matrix(l[i:length(l)])
    cat("processed file ", f, "\n")
}
# process_sig(input_dir, 'PUNQS3.sig')
process_sig <- function(folder, f){
    f = paste(folder, f, sep="")
    file.copy(f, parent_dir, overwrite=TRUE)
    l <- scan(f, what="", sep="\n")
    prodnames <<- strsplit(l[1], '\t', fixed=TRUE)[[1]]
    i <- ifelse(type1, 3, 2)
    sig <<- parse_matrix(l[i:length(l)])
    cat("processed file ", f, "\n")
}
process_input_file <- function(folder, f){
    inf = paste(folder, f, sep="")
    outf = paste(parent_dir, f, sep="")
    # to avoid overwriting calculated files from .DATA
    if(!file.exists(outf))
        file.copy(inf, outf)
    cat("processed file ", f, "\n")
}
# returns name for RSM file in a given folder
locate_rsm <- function(folder) {
    files <- dir(folder)
    f <- files[stri_sub(tolower(files), -4)==".rsm"]
    if(!length(f)) stop("No RSM found.")
    paste(folder, f, sep="")
}
# process_rsm(locate_rsm(parent_dir))
process_rsm <- function(f){
    rsm <<- NULL
    empty <- function(tokens) all(apply(tokens, 2, `==`, ""))
    del_empty <- function(tokens) tokens[-which(tokens == "")]
    del_na <- function(m) m[,apply(m, 2, function(v) !any(is.na(v)))]
    read_block <- function(l) {
        h <- c()
        for(i in seq_along(l)) {
            tokens <- trimws(stri_split_coll(l[i], "\t", simplify=T))
                # there is an empty line between header and matrix
            if(empty(tokens)) break

            # first line of header
            if(!length(h)) h = tokens
            # for other lines glue the header
            else for(i in seq_along(h)) h[i] = paste(h[i], tokens[i])
        }
        h = del_empty(trimws(h))
        header <<- c(header, h)
        blockwidth <- length(h)
        # seek the matrix
        tokens <- trimws(stri_split_coll(l[i], "\t", simplify=T))
        while(empty(tokens)) {
            i = i+1
            tokens <- trimws(stri_split_coll(l[i], "\t", simplify=T))
        }
        m <- parse_matrix(l[i:length(l)])
        # empty cols turn into cols with NA, so delete them
        if(ncol(m) > length(h)) m = del_na(m)
        m
    }
    if(type1) l <- scan(f, what="", sep="\n")
    else l <- scan(f, what="", sep="\n", skipNul=TRUE)
    blockstarts <- which(grepl('SUMMARY', l, fixed=T))
    header <- c()
    for(i in seq_along(blockstarts)) {
        e <- blockstarts[i+1]-1
        b <- read_block(l[(blockstarts[i]+1):ifelse(is.na(e),length(l),e)])
        if(is.null(rsm)) rsm <<- b else rsm <<- cbind(rsm,b)
    }
    # there might be some empty fields in the end of header
    colnames(rsm) <<- header[1:ncol(rsm)]
    cat("processed file ", f, "\n")
}
# NOTE prodnames must already be populated
get_sim <- function() {
    cn <- colnames(rsm)
    # we'll need column with time
    ind <- c(1)
    # name of "TIME" is not needed
    pn <- strsplit(prodnames[-1], ':')
    for(i in seq_along(pn)) {
        w <- pn[[i]]
        pattern <- ifelse(length(w)==2, paste(w[1], ".*", w[2], "$", sep=""), paste("^", w, ".*", sep=""))
        # should always be just one match
        ind = c(ind, which(grepl(pattern, cn, fixed=F)))
    }
    rsm[,ind]
}
# NOTE prodnames, obs, sig, rsm must already be populated and have correct dimensions
misfit_calc <- function() {
    sim <- get_sim()
    for(i in 1:nv) {
        ind <- which(!is.na(obs[,(i+1)]))
        nts <- length(ind)
        times <- obs[ind, 1]
        simind <- unlist(lapply(lapply(times, `==`, sim[,1]), which, TRUE))
        misfits[i] <<- sum((obs[ind, (i+1)] - sim[simind, (i+1)])^2 / (2 * sig[ind, (i+1)]^2)) * weights[i]
    }
}
# dimension check
validate_input <- function() {
    length(prodnames)-1 == nv &&
        ncol(obs) == length(prodnames) &&
        ncol(obs) == ncol(sig) &&
        nrow(obs) == nrow(sig) &&
        prodnames[1] == "TIME" &&
        nrow(rsm) >= nrow(obs) &&
        ncol(rsm) >= ncol(obs) &&
        length(paramnames) == length(ubounds) &&
        length(lbounds) == length(ubounds) &&
        length(paramnames) == varcnt &&
        length(weights) == nv
}
call_simulator <- function() {
    system(paste("$eclipse ", model_name))
    #system(paste("E:/1_HWU/I_Installer/tNav423/tNavigator-con.exe -d ", model_name))
  
    if(type2) system(python_cmd)
}
# processes input files from input_dir and copies them to parent_dir
read_input_files <- function() {
    input <- dir(input_dir)
    for(i in seq_along(input)) {
        f <- input[i]
        # since f will be copied to parent_dir, store its name like this:
        delfiles <<- unique(c(paste(parent_dir, f, sep=""), delfiles))
        switch(tolower(tail(unlist(strsplit(f, '.', fixed=TRUE)), 1)),
            "dat" = process_dat(input_dir, f),
            "data" = process_data(input_dir, f),
            "obs" = process_obs(input_dir, f),
            "sig" = process_sig(input_dir, f),
            process_input_file(input_dir, f)
            )
    }
}
# copies files from parent_dir to output_dir
write_output_dir <- function(iter) {
    outd <- paste(output_dir, iter, sep="")
    dir.create(outd) # fails with warning if outd already exists
    output <- dir(parent_dir)
    for(i in seq_along(output)) {
        f <- output[i]
        if(f != input_dir && f != output_dir) {
            switch(tolower(tail(unlist(strsplit(f, '.', fixed=TRUE)), 1)),
                "r" = {}, # don't touch R files
                "py" = {}, # don't touch Python files
                {file.copy(f, outd, overwrite=TRUE)} # copy everything else
                )
        }
    }
    fc <- file(paste(outd, "ParameterValues", sep="/"))
    writeLines(mapply(paste, paramnames, params, sep="\t"), fc)
    close(fc)
}
# deletes files listed in delfiles array
cleanup <- function() {
    for(i in seq_along(delfiles)) {
        file.remove(delfiles[i])
        cat("deleted file ", delfiles[i], "\n")
    }
    delfiles <<- c()
}
plot_misfits <- function() {
    i1 <- c(1:maxgen)*npop+1
    i2 <- c((i1-1)[-1], nrow(stored_misfits))
    i <- mapply(c, i1, i2)
    mm <- matrix(0, ncol=ncol(i), nrow=nv)
    for(j in 1:ncol(i))
        mm[,j] = apply(stored_misfits[i[1]:i[2],], 2, min)
    png(paste(output_dir, 'misfits.png', sep=''))
    boxplot(mm, xlab="Generation", ylab="Generational Minimum")
    dev.off()
}
plot_params <- function() {
    i <- 1
    # get rid off $ sign
    pn <- substring(paramnames, 2)
    invisible(apply(stored_param, 2, function(v) {
        png(paste(output_dir, pn[i], '.png', sep=''))
        plot(v, xlab="misfit evaluation no", ylab=pn[i])
        dev.off()
        i <<- i+1
    }))
}
# progress_plot <- function(progress) {
#     plot(progress, xlab="Iteration", ylab="Sum of objectives")
# }

wrapper <- function() {
    iter <- 1
    # for progress_plot
    progress <- c()
    objectives <- function(x) {
        params <<- x
        stored_param[iter,] <<- x
        read_input_files()
        call_simulator()
        rsmf <- locate_rsm(parent_dir)
        process_rsm(rsmf)
        stopifnot(validate_input())
        misfit_calc()
        # misfits <- runif(nv)
        stored_misfits[iter,] <<- misfits
        if(is.null(colnames(stored_misfits))) colnames(stored_misfits) <<- prodnames[-1]
        write_output_dir(iter)
        iter <<- iter + 1
        of <- c()
        for(i in seq_along(objf))
            of = c(of, sum(misfits[objf[[i]]]))
#         progress <<- c(progress, sum(of))
#         progress_plot(progress)
        of
    }
    # NOTE hack? .dat file must be read before everything else
    # in order to have bounds for mopsocd call
    input <- dir(input_dir)
    process_dat(input_dir, input[stri_sub(tolower(input), -4)==".dat"])
    colnames(stored_param) <<- substring(paramnames, 2)

#   For MOPSO    
#    s <- algorithm(objectives, varcnt=varcnt, fncnt=objfn,
#	           lowerbound=lbounds, upperbound=ubounds, opt=0, popsize=npop, maxgen=maxgen)

# For RVEA    
    s <- algorithm(objectives, varcnt=varcnt, fncnt=objfn,
            lowerbound=lbounds, upperbound=ubounds, opt=0, popsize=npop, maxgen=maxgen,
            p1=p1, p2=p2, alpha=2, fr=0.1, FE=0)
    
# For NSGAII   
    # s <- algorithm(objectives, varNo=varcnt, objDim=objfn, lowerBounds=lbounds, upperBounds=ubounds,
    #            popSize=npop, tourSize=2, generations=maxgen, cprob=0.7, XoverDistIdx=5, mprob=0.2, MuDistIdx=10)
    
    # if we are here, execution was successful
    cleanup()

    write.table(stored_misfits, file=paste(output_dir,"stored_misfits_rvea.csv",sep=""), col.names=TRUE, row.names=FALSE, sep=',')
    write.table(stored_param, file=paste(output_dir,"stored_param_rvea.csv",sep=""), col.names=TRUE, row.names=FALSE, sep=',')

    plot_params()
}
