.libPaths(c(".R-library",.libPaths()))
args <- commandArgs(trailingOnly=TRUE)
port <- if(length(args)) as.integer(args[1]) else 3876L
shiny::runApp("shiny",host="127.0.0.1",port=port,launch.browser=FALSE)
