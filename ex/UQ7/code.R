library(ReliabilityTheory)

gr = graph.formula(s -- 1 -- 2 -- 3 -- t, 1 -- 4:5 -- 6 -- t,
		s -- 7 -- 8 -- t, s -- 9 -- 10 -- 11 -- t, 7 -- 10 -- 8)

V(gr)$compType[match( c("1","6","11") , V(gr)$name  )] <- "G1" 
V(gr)$compType[match( c("2","3","9") , V(gr)$name  )] <- "G2"
V(gr)$compType[match( c("4","5","10") , V(gr)$name  )] <- "G3"
V(gr)$compType[match( c("7","8") , V(gr)$name  )] <- "G4" 
#V(gr)$compType 

sig <- computeSystemSurvivalSignature(gr)
