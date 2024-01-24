#' ---
#' title: "Get venn diagrams"
#' subtitle: "Longitudinal MVMR in POPS"
#' author: "Janne Pott"
#' date: "Last compiled on `r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'   html_document:
#'     toc: true
#'     number_sections: true
#'     toc_float: true
#'     code_folding: show
#' ---
#'
#' # Introduction ####
#' ***
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../SourceFile_HPC.R")
.libPaths()

tag = format(Sys.time(), "%Y-%m-%d")
tag = gsub("202.-","23-",tag)
tag = gsub("-","",tag)

#' # Load data ####
#' ***
load("../results/08_1_Associations_SNPs_exposure_linMix_231205.RData")
load("../results/08_2_Associations_SNPs_exposure_gamlss_231205.RData")
load("../results/08_5_Associations_SNPs_exposure_gamlss_IA_231211.RData")

myExposures = unique(myAsscos_X_linearMixed$phenotype)
myExposures = myExposures[c(7,12,13)]

data1 = myAsscos_X_linearMixed[phenotype %in% myExposures]
data2 = myAsscos_X_gamlss[phenotype %in% myExposures]
data3 = myAsscos_X_gamlssIA[phenotype %in% myExposures]

#' # Define functions ####
#' ***
#' Taken from Holger Kirsten!

venndiagram <- function(x=x, y=y, z=z, w=w, unique=T, title="Venn Diagram", labels=c("x", "y", "z", "w"), lines=1, lcol=1, tcol=1, diacol=1, plot=T, type="3", printsub=TRUE, ...) {
  ## Remove duplicates and NA fields in x, y, z and w
  if(unique==T) {
    x <- unique(x); x <- as.vector(stats::na.omit(x))
    y <- unique(y); y <- as.vector(stats::na.omit(y))
    if(!missing("z")) {
      z <- unique(z); z <- as.vector(stats::na.omit(z))
    }
    if(!missing("w")) {
      w <- unique(w); w <- as.vector(stats::na.omit(w))
    }
  }
  
  ## Check valid type selection
  if(!type %in% c("2", "2map", "3", "3map", "4", "4map", "4el", "4elmap")) {
    return("Error: the 'type' argument can only be set to one of these values: 2, 2map, 3, 3map, 4, 4map, 4el, 4elmap.")
  }
  
  ## Plot a 2-way venn diagram
  if(type=="2") {
    # Define ovelap queries
    q1 <- x[x %in% y]
    q2 <- x[!x %in% y]
    q3 <- y[!y %in% x]
    
    ## Store query vectors in list
    qlist <- list(q1=q1, q2=q2, q3=q3)
    
    ## Perfom query counts
    count <- unlist(lapply(qlist, length))
    countDF <- data.frame(query=names(count) , count=as.vector(count))
    olDF <- data.frame(x=c(5.0, 3.1, 7.0), y=c(6.1, 6.1, 6.1), count=countDF$count)
    if(printsub==TRUE) {mysub <- paste(paste("N unique: xy =", length(unique(c(x,y)))), paste("; x =", length(unique(x))), paste("; y =", length(unique(y))), sep="")} else {mysub <- ""}
    if(plot==T) {
      ## Plot the 2-way venn diagram
      graphics::symbols(x=c(4, 6), y = c(6, 6), circles=c(2, 2), xlim=c(0, 10), ylim=c(0, 10), inches=F, main=title, sub=mysub, xlab="", ylab="",  xaxt="n", yaxt="n", bty="n", fg=lines, ...);
      graphics::text(olDF$x, olDF$y, olDF$count, col=tcol, ...); graphics::text(c(2.0, 8.0), c(8.8, 8.8), labels[1:2], col=lcol, ...)
    }
    
    ## Return query list
    return(qlist)
  }
  
  ## Plot 2-way mapping venn diagram
  if(type=="2map") {
    olDFdebug <- data.frame(x=c(5.0, 3.1, 7.0), y=c(6.1, 6.1, 6.1), count=paste("q", 1:3, sep=""), ...)
    graphics::symbols(x=c(4, 6), y = c(6, 6), circles=c(2, 2), xlim=c(0, 10), ylim=c(0, 10), inches=F, main="Mapping Venn Diagram", xlab="", ylab="",  xaxt="n", yaxt="n", bty="n", fg=lines, ...);
    graphics::text(olDFdebug$x, olDFdebug$y, olDFdebug$count, col=tcol, ...); graphics::text(c(2.0, 8.0), c(8.8, 8.8), paste(labels[1:2], "=", c("x","y")), col=lcol, ...)
  }
  
  ## Plot a 3-way venn diagram
  if(type=="3") {
    ## Define ovelap queries
    q1 <- x[x %in% y & x %in% z]
    q2 <- x[x %in% z]; q2 <- q2[!q2 %in% y]
    q3 <- y[y %in% z]; q3 <- q3[!q3 %in% x]
    q4 <- x[x %in% y]; q4 <- q4[!q4 %in% z]
    q5 <- x[!x %in% y]; q5 <- q5[!q5 %in% z]
    q6 <- y[!y %in% z]; q6 <- q6[!q6 %in% x]
    q7 <- z[!z %in% x]; q7 <- q7[!q7 %in% y]
    
    ## Store query vectors in list
    qlist <- list(q1=q1, q2=q2, q3=q3, q4=q4, q5=q5, q6=q6, q7=q7)
    
    ## Perfom query counts
    count <- unlist(lapply(qlist, length))
    countDF <- data.frame(query=names(count) , count=as.vector(count))
    olDF <- data.frame(x=c(5.0, 3.8, 6.3, 5.0, 3.0, 7.0, 5.0), y=c(5.6, 4.6, 4.6, 6.9, 6.5, 6.5, 3.0), count=countDF$count)
    if(printsub==TRUE) {mysub <- paste(paste("N unique: xyz =", length(unique(c(x,y,z)))), paste("; x =", length(unique(x))), paste("; y =", length(unique(y))), paste("; z =", length(unique(z))), sep="")} else { mysub <- "" }
    if(plot==T) {
      ## Plot the 3-way venn diagram
      graphics::symbols(x=c(4, 6, 5), y = c(6, 6, 4), circles=c(2, 2, 2), xlim=c(0, 10), ylim=c(0, 10), inches=F, main=title, sub=mysub, xlab="", ylab="",  xaxt="n", yaxt="n", bty="n", fg=lines, ...);
      graphics::text(olDF$x, olDF$y, olDF$count, col=tcol, ...); graphics::text(c(2.0, 8.0, 6.0), c(8.8, 8.8, 1.1), labels[1:3], col=lcol, ...)
    }
    
    ## Return query list
    return(qlist)
  }
  
  ## Plot 3-way mapping venn diagram
  if(type=="3map") {
    olDFdebug <- data.frame(x=c(5.0, 3.8, 6.3, 5.0, 3.0, 7.0, 5.0), y=c(5.6, 4.6, 4.6, 6.9, 6.5, 6.5, 3.0), count=paste("q", 1:7, sep=""), ...)
    graphics::symbols(x=c(4, 6, 5), y = c(6, 6, 4), circles=c(2, 2, 2), xlim=c(0, 10), ylim=c(0, 10), inches=F, main="Mapping Venn Diagram", xlab="", ylab="",  xaxt="n", yaxt="n", bty="n", fg=lines, ...);
    graphics::text(olDFdebug$x, olDFdebug$y, olDFdebug$count, col=tcol, ...); graphics::text(c(2.0, 8.0, 6.0), c(8.8, 8.8, 1.1), paste(labels[1:3], "=", c("x","y","z")), col=lcol, ...)
  }
  
  ## Overlap queries for 4-way venn diagram
  if(type=="4" | type=="4el" | type=="4elmap") {
    ## Define ovelap queries
    xy <- x[x %in% y]; xz <-x[x %in% z]; xw <- x[x %in% w]; yz <- y[y %in% z]; yw <- y[y %in% w]; zw <- z[z %in% w]
    q1 <- xy[xy %in% zw]
    q2 <- xw[xw %in% z]; q2 <- q2[!q2 %in% y]
    q3 <- yz[yz %in% w]; q3 <- q3[!q3 %in% x]
    q4 <- yz[yz %in% x]; q4 <- q4[!q4 %in% w]
    q5 <- xw[xw %in% y]; q5 <- q5[!q5 %in% z]
    q6 <- xy[!xy %in% z]; q6 <- q6[!q6 %in% w]
    q7 <- zw[!zw %in% x]; q7 <- q7[!q7 %in% y]
    q8 <- xz[!xz %in% y]; q8 <- q8[!q8 %in% w]
    q9 <- yw[!yw %in% x]; q9 <- q9[!q9 %in% z]
    q10 <- x[!x %in% c(y,z,w)]
    q11 <- y[!y %in% c(x,z,w)]
    q12 <- z[!z %in% c(x,y,w)]
    q13 <- w[!w %in% c(x,y,z)]
    q14 <- xw[!xw %in% y]; q14 <- q14[!q14 %in% z]
    q15 <- yz[!yz %in% x]; q15 <- q15[!q15 %in% w]
    
    ## Store query vectors in list
    qlist <- list(q1=q1, q2=q2, q3=q3, q4=q4, q5=q5, q6=q6, q7=q7, q8=q8, q9=q9, q10=q10, q11=q11, q12=q12, q13=q13, q14=q14, q15=q15)
    
    ## Perfom query counts
    count <- unlist(lapply(qlist, length))
    countDF <- data.frame(query=names(count) , count=as.vector(count))
    olDF <- data.frame(x=c(4.8, 3.9, 5.7, 3.9, 5.7, 4.8, 4.8, 3.0, 6.5, 3.0, 6.5, 3.0, 6.5, 4.8, 4.8), y=c(5.2, 4.2, 4.2, 6.3, 6.3, 7.2, 3.2, 5.2, 5.2, 7.2, 7.2, 3.2, 3.2, 1.0, 0.4), count=countDF$count)
    
    if(printsub==TRUE) {mysub <- paste(paste("N unique: xyzw =", length(unique(c(x,y,z,w)))), paste("; x =", length(unique(x))), paste("; y =", length(unique(y))), paste("; z =", length(unique(z))), paste("; w =", length(unique(w))), sep="") } else { mysub <- "" }
    
    ## Plot 4-way venn diagram as circles
    if(plot==T & type=="4") {
      graphics::symbols(x=c(4, 5.5, 4, 5.5), y = c(6, 6, 4.5, 4.5), circles=c(2, 2, 2, 2), xlim=c(0, 10), ylim=c(0, 10), inches=F, main=title, sub=mysub, xlab="", ylab="",  xaxt="n", yaxt="n", bty="n", fg=lines, ...);
      graphics::text(olDF$x[1:13], olDF$y[1:13], olDF$count[1:13], col=tcol, ...) # rows 14-15 of olDF are printed in last step
      graphics::text(c(2.0, 7.5, 2.0, 7.5), c(8.3, 8.3, 2.0, 2.0), labels, col=lcol, ...)
      graphics::text(c(3.8, 3.8), c(1.0, 0.4), c(paste("Only in ", labels[1], " & ", labels[4], ": ", olDF$count[14], sep=""), paste("Only in ", labels[2], " & ", labels[3], ": ", olDF$count[15], sep="")), col=diacol, ...)
    }
    
    ## Plot 4-way venn diagram as ellipses
    if(plot==T & (type=="4el" | type=="4elmap")) {
      olDF <- data.frame(x=c(5.0, 4.2, 6.4, 3.6, 5.8, 2.9, 7.1, 3.1, 6.9, 1.5, 3.5, 6.5, 8.5, 5.0, 5.0), y=c(2.8, 1.4, 4.0, 4.0, 1.4, 5.9, 5.9, 2.2, 2.2, 4.8, 7.2, 7.2, 4.8, 0.7, 6.0), count=countDF$count)
      ## Plot ellipse
      plotellipse <- function (center=c(1,1), radius=c(1,2), rotate=1, segments=360, xlab="", ylab="", ...) {
        angles <- (0:segments) * 2 * pi/segments
        rotate <- rotate*pi/180
        ellipse <- cbind(radius[1] * cos(angles), radius[2] * sin(angles))
        ellipse <- cbind( ellipse[,1]*cos(rotate) + ellipse[,2]*sin(rotate), ellipse[,2]*cos(rotate) - ellipse[,1]*sin(rotate) )
        ellipse <- cbind(center[1]+ellipse[,1], center[2]+ellipse[,2])
        graphics::plot(ellipse, type = "l", xlim = c(0, 10), ylim = c(0, 10), xlab = "", ylab = "", ...)
      }
      ## Plot ellipse as 4-way venn diagram
      ellipseVenn <- function(lines=lines, olDF, title=title, labels=labels, sub=mysub, main, lcol=lcol, tcex=1.3, ...) {
        graphics::split.screen(c(1,1))
        plotellipse(center=c(3.5,3.6), radius=c(2,4), rotate=-35, segments=360, xlab="", ylab="", col=lines[1], axes=FALSE, main=title, sub=mysub, ...)
        graphics::screen(1, new=FALSE)
        plotellipse(center=c(4.7,4.4), radius=c(2,4), rotate=-35, segments=360, xlab="", ylab="", col=lines[2], axes=FALSE, ...)
        graphics::screen(1, new=FALSE)
        plotellipse(center=c(5.3,4.4), radius=c(2,4), rotate=35, segments=360, xlab="", ylab="", col=lines[3], axes=FALSE, ...)
        graphics::screen(1, new=FALSE)
        plotellipse(center=c(6.5,3.6), radius=c(2,4), rotate=35, segments=360, xlab="", ylab="", col=lines[4], axes=FALSE, ...)
        graphics::text(olDF[1:15,1], olDF[1:15,2], olDF[1:15,3], col=tcol, ...)
        graphics::text(c(0.4, 2.8, 7.5, 9.4), c(7.3, 8.3, 8.3, 7.3), labels, col=lcol, ...)
        graphics::close.screen(all=TRUE)
      }
      ## Plot 4-way ellipse venn diagram
      if(type=="4el") {
        ellipseVenn(olDF=olDF, lcol=lcol, lines=lines, labels=labels, title=title, ...)
      }
      
      ## Plot 4-way ellipse mapping venn diagram
      if(type=="4elmap") {
        olDFdebug <- data.frame(x=c(5.0, 4.2, 6.4, 3.6, 5.8, 2.9, 7.1, 3.1, 6.9, 1.5, 3.5, 6.5, 8.5, 5.0, 5.0), y=c(2.8, 1.4, 4.0, 4.0, 1.4, 5.9, 5.9, 2.2, 2.2, 4.8, 7.2, 7.2, 4.8, 0.7, 6.0), count=paste("q", 1:15, sep=""), ...)
        ellipseVenn(olDF=olDFdebug, lcol=lcol, lines=lines, labels=paste(labels, "=", c("x","y","z","w")), title="Mapping Venn Diagram", ...)
      }
    }
    
    ## Return query list
    return(qlist)
  }
  
  ## Plot 4-way circle mapping venn diagram
  if(type=="4map") {
    olDFdebug <- data.frame(x=c(4.8, 3.9, 5.7, 3.9, 5.7, 4.8, 4.8, 3.0, 6.5, 3.0, 6.5, 3.0, 6.5, 4.8, 4.8), y=c(5.2, 4.2, 4.2, 6.3, 6.3, 7.2, 3.2, 5.2, 5.2, 7.2, 7.2, 3.2, 3.2, 1.0, 0.4), count=paste("q", 1:15, sep=""), ...)
    graphics::symbols(x=c(4, 5.5, 4, 5.5), y = c(6, 6, 4.5, 4.5), circles=c(2, 2, 2, 2), xlim=c(0, 10), ylim=c(0, 10), inches=F, main="Mapping Venn Diagram", xlab="", ylab="",  xaxt="n", yaxt="n", bty="n", fg=lines, ...);
    graphics::text(olDFdebug$x[1:13], olDFdebug$y[1:13], olDFdebug$count[1:13], col=tcol, ...); graphics::text(c(2.0, 7.5, 2.0, 7.5), c(8.3, 8.3, 2.0, 2.0), paste(labels, "=", c("x","y","z","w")), col=lcol, ...)
    graphics::text(c(3.8, 3.8), c(0.97, 0.36), c(paste("Only in ", labels[1], " & ", labels[4], ": ", olDFdebug$count[14], sep=""), paste("Only in ", labels[2], " & ", labels[3], ": ", olDFdebug$count[15], sep="")), col=tcol, ...)
  }
  
}


venn2 = function(x1,y1, mytitle="2-Way Venn Diagram", mylabels = NA, plotte =T, venntype = '2')
{
  # 28/2/13 plotte par
  # 150119 vector check
  if(all(is.vector(x1) | is.factor(x1),is.vector(y1)|is.factor(y1))==F) stop("All input data must be vectors...")
  if(is.na(mylabels[1])) mylabels = c(deparse(substitute(x1)), deparse(substitute(y1)))
  qlist <- venndiagram(x=x1, y=y1, unique=T, title=mytitle, labels= mylabels, plot=plotte, lines=c(2,3), lcol=c(2,3), tcol=c(1,1,1), lwd=3, cex=1.3, printsub=T, type=venntype)
  qlist
}

venn3 = function(x1,y1,z1, mytitle="3-Way Venn Diagram", mylabels = NA,  plotte =T, venntype = '3')
{
  # 28/2/13 plotte par
  # 150119 vector check
  if(all(is.vector(x1)|is.factor(x1),is.vector(y1)|is.factor(y1),is.vector(z1)|is.factor(z1))==F) stop("All input data must be vectors...")
  
  if(is.na(mylabels[1])) mylabels = c(deparse(substitute(x1)), deparse(substitute(y1)), deparse(substitute(z1)))
  qlist <- venndiagram(x=x1, y=y1, z=z1, unique=T, title=mytitle, labels= mylabels, plot=plotte, lines=c(2,3,4), lcol=c(2,3,4), tcol=c(1,1,1,1,1,1,1), lwd=3, cex=1.3, printsub=T, type=venntype)
  qlist
}

#' # Loop 1 ####
#' ***
#' Venn diagramm per exposure (EFW absolute, Z-score, Centiles): check overlap of associated SNPs between exposure types (mean, var, slope)

dumTab = foreach(i=1:3)%do%{
  #i=1
  set_mean_lmm = data1[pval_mean < 0.05 & phenotype == myExposures[i],SNP]
  set_mean_gam = data2[pval_mean < 0.05 & phenotype == myExposures[i],SNP]
  set_mean_gamIA = data3[pval_mean < 0.05 & phenotype == myExposures[i],SNP]
  
  set_var_gam = data2[pval_var < 0.05 & phenotype == myExposures[i],SNP]
  set_var_gamIA = data3[pval_var < 0.05 & phenotype == myExposures[i],SNP]
  
  set_slope_lmm = data1[pval_var < 0.05 & phenotype == myExposures[i],SNP]
  set_slope_gamIA = data3[pval_slope < 0.05 & phenotype == myExposures[i],SNP]
  
  d1=venn3(x1 = set_mean_gamIA,y1 = set_mean_lmm,z1 = set_mean_gam, 
           mytitle = paste0("Venn ",myExposures[i]," mean p<0.05"),
           mylabels = c("gamlssIA","lmm","gamlss"))
  d2=venn2(x1 = set_var_gamIA,y1 = set_var_gam,
           mytitle = paste0("Venn ",myExposures[i]," var p<0.05"),  
           mylabels = c("gamlssIA","gamlss"))
  d3=venn2(x1 = set_slope_gamIA,y1 = set_slope_lmm,
           mytitle = paste0("Venn ",myExposures[i]," slope p<0.05"),  
           mylabels = c("gamlssIA","lmm"))
  
}


#' # Loop 2 ####
#' ***
#' Venn diagramm per exposure (EFW absolute, Z-score, Centiles): check overlap of associated SNPs between GX methods (lmm, gamlss, gamlssIA)

dumTab = foreach(i=1:3)%do%{
  #i=1
  set_mean_lmm = data1[pval_mean < 0.05 & phenotype == myExposures[i],SNP]
  set_slope_lmm = data1[pval_var < 0.05 & phenotype == myExposures[i],SNP]

  set_mean_gam = data2[pval_mean < 0.05 & phenotype == myExposures[i],SNP]
  set_var_gam = data2[pval_var < 0.05 & phenotype == myExposures[i],SNP]
  
  set_mean_gamIA = data3[pval_mean < 0.05 & phenotype == myExposures[i],SNP]
  set_slope_gamIA = data3[pval_slope < 0.05 & phenotype == myExposures[i],SNP]
  set_var_gamIA = data3[pval_var < 0.05 & phenotype == myExposures[i],SNP]
  
  d1=venn3(x1 = set_mean_gamIA,y1 = set_var_gamIA,z1 = set_slope_gamIA, 
           mytitle = paste0("Venn ",myExposures[i]," mean p<0.05"),
           mylabels = c("mean","var","slope"))
  d2=venn2(x1 = set_mean_lmm,y1 = set_slope_lmm,
           mytitle = paste0("Venn ",myExposures[i]," var p<0.05"),  
           mylabels = c("mean","slope"))
  d3=venn2(x1 = set_mean_gam,y1 = set_var_gam,
           mytitle = paste0("Venn ",myExposures[i]," slope p<0.05"),  
           mylabels = c("mean","var"))
  
}


#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

