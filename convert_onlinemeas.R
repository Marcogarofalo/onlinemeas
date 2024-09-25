# args <- commandArgs(trailingOnly = TRUE)
#
# if (length(args) < 1) {
#   stop("usage: Rscript convert.R  infile", call. = FALSE)
# }

# source(args[1])
library("fs")
source("merging_info.R")

read_LTncorr <- function(df) {
  size <- dim(df)
  L <- df[size[1], 3]
  T <- L * 2
  ncorr <- size[1] / (L + 1)
  return(list("L" = L, "T" = T, "ncorr" = ncorr, "gammas" = unique(df[, 1]),
              "size" = ncorr * 2 * T))
}
compare_LTncorr <- function(head, head1) {
  stopifnot(
    head$L == head1$L, head$T == head1$T, head$ncorr == head1$ncorr,
    all(head$gamma == head1$gamma)
  )
}

files <- dir_ls(indir)
id_online <- grep("onlinemeas", files)
confs <- unique(gsub(paste0(indir, "onlinemeas\\.s...\\."), "", files[id_online]))

### init data from the first onlinemeas
df <- read.table(files[id_online][1])
head <- read_LTncorr(df)

############ printing header

con <- file(out, "wb")
Sint <- 4
writeBin(as.integer(length(confs)), con, size = Sint, endian = "little")
writeBin(as.integer(head$T), con, size = Sint, endian = "little")
writeBin(as.integer(head$L), con, size = Sint, endian = "little")
writeBin(as.integer(head$ncorr), con, size = Sint, endian = "little")

Sdoub <- 8
writeBin(as.numeric(beta), con, size = Sdoub, endian = "little")
writeBin(as.numeric(kappa), con, size = Sdoub, endian = "little")

writeBin(length(mus), con, size = Sint, endian = "little")
for (mu in mus) {
  writeBin(as.numeric(mu), con, size = Sdoub, endian = "little")
}
writeBin(length(rs), con, size = Sint, endian = "little")
for (r in rs) {
  writeBin(as.numeric(r), con, size = Sdoub, endian = "little")
}
writeBin(length(thetas), con, size = Sint, endian = "little")
for (r in thetas) {
  writeBin(as.numeric(r), con, size = Sdoub, endian = "little")
}
writeBin(length(head$gammas), con, size = Sint, endian = "little")
for (r in paste(head$gammas)) {
  writeChar(r, con, nchars = nchar(r, type = "chars"), eos = "", useBytes = FALSE)
}
writeBin(length(smearing), con, size = Sint, endian = "little")
for (r in smearing) {
  writeChar(r, con, nchars = nchar(r, type = "chars"), eos = "", useBytes = FALSE)
}
writeBin(length(bananas), con, size = Sint, endian = "little")
for (r in bananas) {
  writeBin(as.integer(r), con, size = Sint, endian = "little")
}
writeBin(length(oranges), con, size = Sint, endian = "little")
for (r in oranges) {
  writeBin(as.numeric(r), con, size = Sdoub, endian = "little")
}
writeBin(as.integer(head$size), con, size = Sint, endian = "little")


############## alocate data for average
data <- array(rep(0, head$T * head$ncorr), c(head$ncorr, head$T))
first <- c(1:(head$L+1))
second <- c((head$L + 2):head$T)



for (conf in confs) {
  cat(conf,"\t")
  ids_src <- grep(conf, files[id_online])

  ## source average
  data[, ] <- 0
  for (id_src in ids_src) {
    df <- read.table(files[id_src])
    head1 <- read_LTncorr(df)
    compare_LTncorr(head, head1)
    for (id_cor in seq_along(head$gammas)) {
      rows_forward <- c(1:(head$L+1)) + (id_cor - 1) * (head$L + 1)
      data[id_cor, first] <- data[id_cor, first] + df[rows_forward, 4]
      rows_backward <- c(head$L:2) + (id_cor - 1) * (head$L + 1)
      data[id_cor, second] <- data[id_cor, second] + df[rows_backward, 5]
    }
  }
  ## printing
  writeBin(as.integer(conf), con, size = Sint, endian = "little")
  for (id_cor in seq_along(head$gammas)) {
    for (t in c(1:T)) {
      
      writeBin(data[id_cor,t], con, size = Sdoub, endian = "little")
      writeBin(as.numeric(0), con, size = Sdoub, endian = "little")
    }
  }
}

cat("\n")
close(con)
