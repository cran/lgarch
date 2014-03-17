.LGARCHRECURSION1 <- function (iStart, iEnd, phi1, theta1, yzeroadj, innov, 
    lny2adj, uadj) 
.C("LGARCHRECURSION1", iStart = as.integer(iStart), 
    iEnd = as.integer(iEnd), phi1 = as.double(phi1), theta1 = as.double(theta1), 
    yzeroadj = as.double(yzeroadj), innov = as.double(innov), 
    lny2adj = as.double(lny2adj), uadj = as.double(uadj), PACKAGE = "lgarch")
