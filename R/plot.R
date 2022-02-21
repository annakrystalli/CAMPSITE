plot_CAMPSITE <- function(traitmat, linmat, step_size, ylims = NULL) {
  #plots a simulation, with incipient lineages in red, good in black
  steps <- max(linmat[, 4])/step_size #calculate number of steps in simulation
  max_trait <- NULL
  min_trait <- NULL
  for (i in 1:nrow(linmat)) {
    max_trait <- c(max_trait,
                   max(traitmat[[i]][-(1:4)], 
                       na.rm = T)) #extract maximal trait values for each lineage
    min_trait <- c(min_trait,
                   min(traitmat[[i]][-(1:4)], 
                       na.rm = T)) #extract minimal trait values for each lineage
  }
  if (is.null(ylims)) {
    plot(1,
         type = "n",
         xlim = c(1, steps),
         ylim = c(min(min_trait), 
                  max(max_trait)),
         ylab = "trait value",
         xlab = "Time") #generate plot area
  }
  else {
    plot(1, type = "n",
         xlim = c(1, steps),
         ylim = ylims, 
         ylab = "trait value",
         xlab = "Time") #generate plot area with pre-determined limits for trait space
  }
  completion <- linmat[, 6] #extract vector of speciation or extinction times for each lineage
  #handle and plot each lineage based on its end point:
  for (i in 1:length(completion)) {
    if (is.na(completion[i])) { #extinct incipient lineages
      lines(x = (5:length(traitmat[[i]])), #branch length
            y = traitmat[[i]][5:length(traitmat[[i]])], #trait values
            col = "red", #incipient lineages in red
            cex = 0.5)
    }
    else if (linmat[i, 5] == -2 &&
             !is.na(linmat[i, 7]) &&
             completion[i] != linmat[i, 7]) { #extinct good lineages
      lines(x = 5:((linmat[i, 7]/step_size) + 5), #branch length of incipient stage
            y = traitmat[[i]][5:((linmat[i, 7]/step_size) + 5)], #trait values of incipient stage
            col = "red", #incipient stage in red
            cex = 0.5)
      lines(x = ((linmat[i, 7]/step_size) + 5):length(traitmat[[i]]), #branch length of good stage
            y = traitmat[[i]][((linmat[i, 7]/step_size) + 5):length(traitmat[[i]])], #trait value of good stage
            col = "black", #good stage in black
            cex = 0.5)
    }
    else { #extant good and incipient lineages
      lines(x = 5:((completion[i]/step_size) + 5), #branch length of incipient stage
            y = traitmat[[i]][5:((completion[i]/step_size) + 5)], #trait values of incipient stage
            col = "red", #incipient stage in red
            cex = 0.5)
      if (length(traitmat[[i]]) > (round(completion[i]/step_size) +  5)) { #if completed speciation before end of process
        lines(x = ((completion[i]/step_size) + 5):length(traitmat[[i]]), #branch length of good stage
              y = traitmat[[i]][((completion[i]/step_size) + 5):length(traitmat[[i]])], #trait values of good stage
              col = "black", #good stage in black
              cex = 0.5)
      }
    }
  }
}
