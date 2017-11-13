safeColor <- function(motifs){
  DNAcolors <- c(A="#009E73", C="#0072B2", G="#E69F00", 'T'="#D55E00")
  RNAcolors <- c(A="#009E73", C="#0072B2", G="#E69F00", 'U'="#D55E00")
  map <- c("#00811B"="#009E73", "#2000C7"="#0072B2", "#FFB32C"="#E69F00", "#D00001"="#D55E00")
  if(is(motifs, "list")){
    stopifnot(all(sapply(motifs, function(.ele) inherits(.ele, c("pcm", "pfm", "psam")))))
    lapply(motifs, function(.ele){
      if(all(names(.ele$color) %in% names(DNAcolors))){
        .ele$color <- DNAcolors
      }else{
        if(all(names(.ele$color) %in% names(RNAcolors))){
          .ele$color <- RNAcolors
        }else{
          .col <- map[.ele$color]
          .col[is.na(.col)] <- .ele$color[is.na(.col)]
          names(.col) <- names(.ele$color)
          .ele$color <- .col
        }
      }
      .ele
    })
  }else{
    stopifnot(inherits(motifs, c("pcm", "pfm", "psam")))
    if(all(names(motifs$color) %in% names(DNAcolors))){
      motifs$color <- DNAcolors
    }else{
      if(all(names(motifs$color) %in% names(RNAcolors))){
        motifs$color <- RNAcolors
      }else{
        .col <- map[.ele$color]
        .col[is.na(.col)] <- .ele$color[is.na(.col)]
        names(.col) <- names(.ele$color)
        .ele$color <- .col
      }
    }
    motifs
  }
}