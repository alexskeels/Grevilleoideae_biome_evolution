
setUpETT2 <- function (states, CET) 
{
  state.no <- length(states)
  ETT <- data.frame(time = NA, node = 0, event = NA, extant_sp = 0)
  for (j in 1:state.no) {
    ETT <- cbind(ETT, states[j])
  }
  names(ETT)[5:length(ETT)] <- as.character(states)
  ETT[1, 5:length(ETT)] <- rep(0, times = state.no)
  root.node <- CET[CET$node.type == "root", "node"]
  root.value <- CET[CET$node.type == "root", "states"]
  ETT[1, 1] <- CET[CET$node.type == "root", "time_bp"]
  ETT[1, "node"] <- root.node
  ETT
}


getTableIndexFromEventTime2 <- function (CET, ET, ETT) {
  if (length(which(CET$time_bp %in% ET)) == 1) {
    index <- which(CET$time_bp %in% ET)
    node <- CET$node[index]
  }
  else {
    index <- which(CET$time_bp %in% ET & !CET$node %in% 
                     ETT$node)[1]
    node <- CET$node[index]
  }
  return(c(index, node))
}

nextEventTime2 <- function (CET,  ET) {
  event.timing <- c(CET$time_bp)
  order.event <- event.timing[order(event.timing)]
  nextEventTime <- order.event[which(order.event == ET) - 1]
  if (order.event[which(order.event == ET)] == head(order.event)[1]) {
    nextEventTime = 0
  }
  return(nextEventTime)
}



whichStateChangesClado2 <- function (x, ETT, ET, index) {
  node <- x$node[index]
  daughters <- x$node[which(x$parent==node)]
  
  if (is.na(x$states[index])) {
    return(ETT)
  }
  
  # if the state of both descendants is the same as parent - increase state by 1
  if (x$states[index] == x$states[which(x$node %in% daughters)][1] & 
      x$states[index] == x$states[which(x$node %in% daughters)][2]) {
    state.value <- as.numeric(ETT[which(ETT[, 1] == ET), which(names(ETT) == as.character(x$states[index]))])
    ETT[which(ETT[, 1] == ET), which(names(ETT) == as.character(x$states[index]))] <- state.value +  1
    return(ETT)
  }
  
  # if the state of both descendants is different same as parent - decrease OG state by 1 and add 1 to each daughter state
  if (x$states[index] != x$states[which(x$node %in% daughters)][1] & 
      x$states[index] != x$states[which(x$node %in% daughters)][2]) {
    
    # minus 1 from OG
    state.value <- as.numeric(ETT[which(ETT[, 1] == ET), which(names(ETT) ==  as.character(x$states[index]))])
    ETT[which(ETT[, 1] == ET), which(names(ETT) == as.character(x$states[index]))] <- state.value - 1
    
    # add 1 to each new species
    state.value <-as.numeric( ETT[which(ETT[, 1] == ET), which(names(ETT) == as.character(x$states[which(x$node %in% daughters)][1]))])
    ETT[which(ETT[, 1] == ET), which(names(ETT) == as.character(x$states[which(x$node %in% daughters)][1]))] <- state.value + 1
    state.value <- as.numeric(ETT[which(ETT[, 1] == ET), which(names(ETT) == as.character(x$states[which(x$node %in% daughters)][2]))])
    ETT[which(ETT[, 1] == ET), which(names(ETT) == as.character(x$states[which(x$node %in% daughters)][2]))] <- state.value + 1
    return(ETT)
  }
  
  # if state of one daughter differs, add to that daughter
  if (x$states[index] != x$states[which(x$node %in% daughters)][2] & 
      x$states[index] == x$states[which(x$node %in% daughters)][1]) {
    state.value <- as.numeric(ETT[which(ETT[, 1] == ET), which(names(ETT) == as.character(x$states[which(x$node %in% daughters)][2]))])
    ETT[which(ETT[, 1] == ET), which(names(ETT) == as.character(x$states[which(x$node %in% daughters)][2]))] <- state.value +  1
    return(ETT)
  }
  
  # if state of the other daughter differs
  if (x$states[index] == x$states[which(x$node %in% daughters)][2] & 
      x$states[index] != x$states[which(x$node %in% daughters)][1]) {
    state.value <- as.numeric(ETT[which(ETT[, 1] == ET), which(names(ETT) ==  as.character(x$states[which(x$node %in% daughters)][1]))])
    ETT[which(ETT[, 1] == ET), which(names(ETT) == as.character(x$states[which(x$node %in% daughters)][1]))] <- state.value +  1
    return(ETT)
  }
}

getEventTiming2 <- function (tree, CET) {
  CET_k <- CET
  states <- unique(c(CET_k$states))
  ETT <- setUpETT2(states, CET_k)
  event.time <- ETT[1, 1]
  index <- which(CET_k$node.type == "root")
  for (j in 1:nrow(CET_k)) {
    index.node <- getTableIndexFromEventTime2(CET_k, ET = event.time,  ETT)
    index <- index.node[1]
    
    if (CET_k$node.type[index] == "tip") {
      event.time <- nextEventTime2(CET_k,  event.time)
      next
    } else {
      ETT <- whichStateChangesClado2(CET_k, ETT, event.time, 
                                     index)
    }
    ETT[length(ETT[, 1]), "node"] <- index.node[2]
    ETT[length(ETT[, 1]), "extant_sp"] <- sum(as.numeric(ETT[length(ETT[, 1]), 5:length(ETT)]))
    event.time <- nextEventTime2(CET_k, event.time)
    ETT <- rbind(ETT, ETT[length(ETT[, 1]), ], make.row.names = F)
    ETT[length(ETT[, 1]), "time"] <- event.time
  }
  ETT
}

correctNodeHeights2<- function (CET, sd = 0.01) {
  CET$time_bp <- CET$time_bp + rnorm(length(CET$time_bp),  0, sd)
  while (any(duplicated(CET$time_bp))) {
    CET$time_bp <- CET$time_bp + rnorm(length(CET$time_bp),  0, sd)
  }
  return(CET)
}

getLTSTDataTable2 <- function (ETT) {
  require(stringr)
  ranges <- colnames(ETT)[5:ncol(ETT)]
  unique.biogeo <- unique(unlist(stringr::str_extract_all(ranges, 
                                                          boundary("character"))))
  ETT2 <- as.data.frame(matrix(ncol = (length(unique.biogeo) + 
                                         1), nrow = length(ETT[, 1])))
  colnames(ETT2) <- c("time", unique.biogeo)
  ETT_tmp <- ETT
  for (j in 1:length(unique.biogeo)) {
    
    region_richness <- ETT_tmp[5:length(ETT_tmp)][grepl(unique.biogeo[j], colnames(ETT_tmp)[5:length(ETT_tmp)])]
    region_richness <- sapply(region_richness,  as.numeric)
    ETT2[, which(colnames(ETT2) %in% unique.biogeo[j])] <- rowSums(region_richness)
  }
  ETT2$time <- ETT$time
  ETT2
}
