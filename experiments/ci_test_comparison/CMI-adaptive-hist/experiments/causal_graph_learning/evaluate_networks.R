##### Helper functions that help to evaluate the correctness of a network

getEdges = function(res, nnodes,V){
    g = res@graph@edgeL
    pcList = list()
    from = c()
    to = c()
    for(i in 1:nnodes){
        ne = g[[V[i]]]$edges
        if(length(ne) > 0){
            for(c in V[ne]){
                from = c(from, V[i])
                to = c(to, c)
            }
        }
    }
    arcs = cbind(from, to, dir=rep("->", length(to)))
    arcs = condenseArrows(arcs)
    return(arcs)
}
isInConflictWithEdge = function(test, ref){
    if(test[3] == "--" | ref[3] == "--"){
        # one of them is undirected
        if((test[1] == ref[1] & test[2] == ref[2]) | (test[1] == ref[2] & test[2] == ref[1])){
            return(TRUE)
        }else{
            return(FALSE)
        }
    }else{
        # both directed
        if(test[1] == ref[2] & test[2] == ref[1]){
            return(TRUE)
        }else{
            return(FALSE)
        }
    }
}
Dim = function(arcs){
    if(is.null(arcs)){
        return(c(0,0))
    }else{
        d = dim(arcs)
        if(is.null(d)){
            d2 = dim(data.frame(t(arcs)))
            if(is.null(d2)){
                return(c(0,0))
            }else{
                return(d2)
            }
        }else{
            return(d)
        }
    }
}
condenseArrows = function(arcs){
    t = 1
    if(is.null(arcs)){
        return(arcs)
    }else{
        while(t <= Dim(arcs)[1] & Dim(arcs)[1] > 1){
            currConflicts = 0
            currArc = arcs[t,]
            keep = rep(TRUE, dim(arcs)[1])
            if(currArc[3] == "--"){
                if(t+1 <= Dim(arcs)[1]){
                    for(i in (t+1):(Dim(arcs)[1])){
                        if(isInConflictWithEdge(arcs[i,], currArc)){
                            keep[i] = FALSE
                        }
                    }
                }
            }else{
                if(t+1 <= Dim(arcs)[1]){
                    for(i in (t+1):(Dim(arcs)[1])){
                        if(isInConflictWithEdge(arcs[i,], currArc)){
                            if(currConflicts == 0){
                                arcs[t,3] = "--"
                                currArc = arcs[t,]
                            }else{
                                keep[i] = FALSE
                            }
                            currConflicts = currConflicts + 1
                        }
                    }
                }
            }
            arcs = arcs[keep,]
            if(sum(keep == FALSE) > 0 | currConflicts > 0){
                t = t
            }else{
                t = t + 1
            }
        }}
    return(arcs)
}
edgeExistsUndir = function(edge, arcs){
    for(i in 1:dim(arcs)[1]){
        arc = arcs[i,]
        if(sum(edge[1:2] %in% arc) == 2){
            return(TRUE)
        }
    }
    return(FALSE)
}
compareToTrueNet = function(arcs, trueArcs, score=edgeExistsUndir){
    ### get true network
    tp = 0
    if(Dim(arcs)[1] > 0){
        if(Dim(arcs)[1] == 1){
            if(score(edge=arcs, arcs=trueArcs)){
                tp = tp + 1
            }
        }else{
            for(i in 1:dim(arcs)[1]){
                if(score(edge=arcs[i,], arcs=trueArcs)){
                    tp = tp + 1
                }
            }
        }
    }
    if(Dim(arcs)[1] == 1){
        return(list(found=1, tp=tp, correct=dim(trueArcs)[1]))
    }else{
        return(list(found=dim(arcs)[1], tp=tp, correct=dim(trueArcs)[1]))
    }
}

F1 = function(prec, rec){
  if(prec + rec == 0){
    return(0)
  }else{
    return((2 * prec * rec)/(prec + rec))
  }
}
