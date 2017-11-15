safeColor <- function(motifs){
  #DNAcolors <- c(A="#009E73", C="#0072B2", G="#E69F00", 'T'="#D55E00")
  #RNAcolors <- c(A="#009E73", C="#0072B2", G="#E69F00", 'U'="#D55E00")
  #map <- c("#00811B"="#009E73", "#2000C7"="#0072B2", "#FFB32C"="#E69F00", "#D00001"="#D55E00")
  DNAcolors <- c(A="#40E0D0", C="#2000C7", G="#FFB32C", 'T'="#FF00FF")
  RNAcolors <- c(A="#40E0D0", C="#2000C7", G="#FFB32C", 'U'="#FF00FF")
  map <- c("#00811B"="#40E0D0", "#2000C7"="#2000C7", "#FFB32C"="#FFB32C", "#D00001"="#FF00FF")
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


# colorSet <- c("Dm"="#009E73", "b1h"="#009E73", "sw"="#008080", "bml"="darkgreen", 
#               "Mm"="#D55E00", "MmDREAM"="#56B4E9", "Ms"="#E69F00", 
#               "Hs"="#CC79A7")

colorSet <- c("Dm"="#40E0D0", "b1h"="#40E0D0", "sw"="#008080", "bml"="darkgreen", 
              "Mm"="brown", "MmDREAM"="blue", "Ms"="#F69156", 
              "Hs"="#D900D9")