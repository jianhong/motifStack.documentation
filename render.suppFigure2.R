library(rmarkdown)
render('suppFigure2.Rmd')
x <- readLines('suppFigure2.tex')
pos <- grep('begin\\{figure\\}', x)
x[pos] <- gsub('begin\\{figure\\}', 'begin{figure}[!h]', x[pos])
writeLines(x, 'suppFigure2.tex')
tools::texi2pdf('suppFigure2.tex', clean = TRUE)

file.remove('suppFigure2.tex')
unlink("figure", recursive = TRUE)
unlink("suppFigure2_files", recursive = TRUE)

