#!/usr/bin/env Rscript
require( graphics )

df <- read.table("1d.txt")
acf( df, lag.max = 50000, type = c("correlation"), plot = TRUE )

