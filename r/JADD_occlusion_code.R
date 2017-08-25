# !!! Make sure wd() is ./r/ !!! #

## ---- options ----
options(digits=3)
options(scipen=999)
options(contrasts=c('contr.sum', 'contr.poly'))

## ---- loadPackages ----
# Note: many are redundant and unused, this is taken directly from some general code from my thesis
Check.Packages <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE, repos="http://cran.rstudio.com/")
  sapply(pkg, require, character.only = TRUE)
}
packages <- c("ggplot2", "afex", "Hmisc", "pander", "plyr", "memisc", "apsrtable", "nFactors", 
              "psych", "data.table", "aspace", "moments", "diptest", "QuantPsyc", "relaimpo", "xtable",
              "pequod", "lmSupport", "lsmeans", "MPDiR","modelfree", "missForest", "lavaan", "mokken",
              "MissMech", "mediation", "reshape2", "tidyr", "gridExtra", "dplyr", "plotrix", "ggExtra",
              "devtools", "GGally", "VIM", "stringr", "cowplot", "MVN", "royston", "RcmdrMisc", "VennDiagram")
Check.Packages(packages)

## ---- loadFunctions ----
# Some of these function may not be used in the code below, as they are "helper" functions used throughout my thesis
# Replace missing crosstabulated factor levels with NA
findDiff<-function(x, y){
  x.p<-do.call("paste", x)
  y.p<-do.call("paste", y)
  x[!x.p %in% y.p,]$x<-NA
  return(x)
}

# Round all numeric variables in a data frame
roundDF <- function(x, digits) {
  # x: data frame 
  # digits: number of digits to round
  numericColumns <- sapply(x, class) == 'numeric'
  x[numericColumns] <-  round(x[numericColumns], digits)
  x
}

# Highlight significant p values
highlightSig<-function(x){
  ifelse(x[,'p']<0.05 | x[,'p'] == "<0.001", paste0("\\textbf{",x[,'p'],"}"), x[,'p']) # highlight p values < 0.05
}

# Function for standardised coefficients
lm.beta.lmer <- function(mod) {
  b <- fixef(mod)[-1]
  sd.x <- sd(getME(mod,"X")[,-1])
  sd.y <- sd(getME(mod,"y"))
  b*sd.x/sd.y
}

# Match multiple graph width, even when only one has axis labels
matchGridWidth<-function(x,y) {
  x<-ggplotGrob(x)
  y<-ggplotGrob(y)
  maxWidth <- grid::unit.pmax(x$widths[2:5], y$widths[2:5])
  x$widths[2:5]<-as.list(maxWidth)
  y$widths[2:5]<-as.list(maxWidth)
  p <- arrangeGrob(x, y, ncol = 2)
}

# Create pretty ("formatted") and functional ("raw") tables for mixed model objects
mixedModelTables<-function(x) {
  x.table<-summary(x$full.model)$"coefficients"
  x.table<-cbind(x.table, append(as.numeric(x$anova.table$p.value), "", 0))
  x.table<-cbind(x.table, append(as.numeric(lm.beta.lmer(x$full.model)), "", 0))
  x.table<-apply(x.table,2,as.numeric)
  colnames(x.table)<-c("B", "B Std. Error", "t value", "p", "beta")
  x.table<-x.table[,c(5,1:4)]
  rownames(x.table)<-c("Intercept",x$anova.table$Effect)
  x.table[,'p']<-format.pval(x.table[,'p'], digits = 3, eps = 0.001)
  x.table.rownames<-rownames(x.table)
  x.table<-as.data.frame(x.table, stringsAsFactors = FALSE)
  rownames(x.table)<-x.table.rownames
  x.table[,1:4]<-sapply(x.table[,1:4], as.numeric)
  x.table.raw<-x.table
  x.table[,'p']<-highlightSig(x.table)
  names(x.table)[names(x.table) == 'beta'] <- "$\\beta$"
  x.table.trunc <- x.table.raw[which(x.table.raw$p<0.05),]
  outputTables <- list("formatted" = x.table, "raw" = x.table.raw, "trunc" = x.table.trunc)
}

afexMixedTables<-function(x, adj.p = NULL, adj.p.n = NULL){
  if(!is.null(adj.p) & is.null(adj.p.n)){
    return("If adjusting p values, remember to include [adj.p.n]!")
  } else if(!is.null(adj.p) & !is.null(adj.p.n)){
    x.table.anova<-x$anova.table[,1:6]
    x.table.anova$p.value <- sapply(x.table.anova$p.value, function(p) p.adjust(p, method = "bonferroni", n = adj.p.n))
  } else if(is.null(adj.p)){
    x.table.anova<-x$anova.table[,1:6]
  }
  x.table.lme<-summary(x$full.model)$coefficients
  x.table.lme<-cbind.data.frame(Effect = rownames(summary(x$full.model)$coefficients), beta = append(lm.beta.lmer(x$full.model), NA,0), x.table.lme)
  rownames(x.table.lme)<-NULL
  x.table.lme<-as.data.frame(x.table.lme, stringsAsFactors = FALSE)
  colnames(x.table.lme)<-c("Effect", "beta", "B", "B Std. Error", "t value")
  x.table.lme <- roundDF(x.table.lme, 3)
  
  x.table.anova$ndf<-x.table.anova$ddf <- NULL
  x.table.anova$df<-paste(x$anova.table$ndf, round(x$anova.table$ddf, 3), sep = ", ")
  colnames(x.table.anova)<-c("Effect", "F", "F scaling", "p", "df")
  x.table.anova[,'p']<-format.pval(x.table.anova[,'p'], digits = 3, eps = 0.001)
  x.table.anova<-x.table.anova[,c(1:2,5,3:4)]
  x.table.anova.raw<-x.table.anova
  x.table.anova[,'p']<-highlightSig(x.table.anova)
  x.table.anova <- roundDF(x.table.anova, 3)
  
  x.table.anova.trunc <- x.table.anova.raw[which(x.table.anova.raw$p<0.05),]
  outputTables <- list("lme" = x.table.lme, "formatted" = x.table.anova, "raw" = x.table.anova.raw, "trunc" = x.table.anova.trunc)
  
  return(outputTables)
}

reportFTest <- function(table, effect){
  df <- table$raw[which(table$raw$Effect == effect), "df"]
  fValue <- table$raw[which(table$raw$Effect == effect), "F"]
  pValue <- table$raw[which(table$raw$Effect == effect), "p"]
  paste("F(", df, ") = ", round(fValue, 3), ", p = ", pValue, sep ="")
}

reportTTest <- function(posthoc, contrast, 
                        xName = NULL, xCase = NULL, 
                        yName = NULL, yCase = NULL){
  if((!is.null(xName) & is.null(xCase)) | (!is.null(yName) & is.null(yCase))){
    return("Remember to include the condition(s) you want to perform the t-test on!")
  } else if(is.null(xName) & is.null(yName)){
    posthoc.data.frame <- data.frame(summary(posthoc$contrasts))
  } else if(is.null(xName)){
    posthoc.data.frame <- data.frame(summary(posthoc$contrasts))[which(data.frame(summary(posthoc$contrasts))
                                                                       [yName] == yCase),]
  } else if(is.null(yName)){
    posthoc.data.frame <- data.frame(summary(posthoc$contrasts))[which(data.frame(summary(posthoc$contrasts))
                                                                       [xName] == xCase),]
  }
  df <- posthoc.data.frame[which(posthoc.data.frame$contrast == contrast), "df"]
  tValue <- posthoc.data.frame[which(posthoc.data.frame$contrast == contrast), "t.ratio"]
  pValue <- posthoc.data.frame[which(posthoc.data.frame$contrast == contrast), "p.value"]
  paste("t(", round(df, 3), ") = ", round(tValue, 3), ", p = ", format.pval(pValue, digits = 3, eps = 0.001), sep ="")
}

# Clean table columns by removing duplicates
cleanTableColumns <- function(x){     
  oldx <- c(FALSE, x[-1]==x[-length(x)])  
  # is the value equal to the previous?    
  res <- x
  res[oldx] <- NA
  return(res)
}

# Make whitespace in posthoc tables LaTeX-compatible
whitespaceLaTeXCompatible <- function(x){
  trim.leading <- function (x)  sub("^\\s+", "", x)
  add.whitespace.trailing <- function (x) sub("\\s+$", "\\\\hphantom{x}", x)
  x<-trim.leading(x)
  x<-add.whitespace.trailing(x)
  return(x)
}

# Create nice correlation table with p-value identifiers
corstars <- function(x, corrMethod = "pearson"){
  require(Hmisc)
  x <- as.matrix(x)
  R <- rcorr(x, type = corrMethod)$r
  p <- rcorr(x, type = corrMethod)$P
  mystars <- ifelse(p < .001, "***", ifelse(p < .01, "** ", ifelse(p < .05, "* ", " ")))
  R <- format(round(cbind(rep(-1.11, ncol(x)), R), 2))[,-1]
  Rnew <- matrix(paste(R, mystars, sep=""), ncol=ncol(x))
  diag(Rnew) <- paste(diag(R), " ", sep="")
  rownames(Rnew) <- colnames(x)
  colnames(Rnew) <- paste(colnames(x), "", sep="")
  Rnew <- as.data.frame(Rnew)
  return(Rnew)
}

# Adjust correlation matrix p-values for multi comparisons
corstars.adjusted <- function(x, corrMethod = "pearson"){
  require(Hmisc)
  x <- as.matrix(x)
  R <- rcorr.adjust(x, type = corrMethod)$R$r
  p <- rcorr.adjust(x, type = corrMethod)$P
  p[p=="<.0001"]<-"0"
  p <- apply(p,1,as.numeric)
  rownames(p) <- colnames(p)
  mystars <- ifelse(p < .001, "***", ifelse(p < .01, "** ", ifelse(p < .05, "* ", " ")))
  R <- format(round(cbind(rep(-1.11, ncol(x)), R),2))[,-1]
  Rnew <- matrix(paste(R, mystars, sep=""), ncol=ncol(x))
  diag(Rnew) <- "{--}"
  rownames(Rnew) <- colnames(x)
  colnames(Rnew) <- paste(colnames(x), "", sep="")
  Rnew <- as.data.frame(Rnew)
  Rnew[upper.tri(Rnew)] <- NA
  return(Rnew)
}

# Center (selected subset) of columns. Use cbind() to add these to current dataset
centerColumns <- function(x){
  xCenter <- x - rep(colMeans(x, na.rm=T), rep.int(nrow(x), ncol(x)))
  colnames(xCenter) <- paste("c", colnames(x), sep = ".")
  return(xCenter)
}

## ---- occlusion2-data-import ----
occlusion2.data<-NA
occlusion2.data$psychometric.parameters <- read.csv("../data/JADD_occlusion_PsychometricData.csv", stringsAsFactors = FALSE)
occlusion2.data$reactionData <- read.csv("../data/JADD_occlusion_ReactionData.csv", stringsAsFactors = FALSE)
occlusion2.data$descriptives <- read.csv("../data/JADD_occlusion_Descriptives.csv", stringsAsFactors = FALSE)

# Aggregated data can be found in '..\data\aggregate\JADD_occlusion_data_aggregate.csv'

## ---- occlusion2-data-munge ----
# Remove TD participants who have larger-than-threshold AQ score
occlusion2.td.highAQ <- occlusion2.data$descriptives[which(occlusion2.data$descriptives$group != "asd" & occlusion2.data$descriptives$AQ>=32), "ID"]
occlusion2.data$psychometric.parameters.trunc<-occlusion2.data$psychometric.parameters[which(!occlusion2.data$psychometric.parameters$ID %in% occlusion2.td.highAQ),]
occlusion2.data$reactionData.trunc<-occlusion2.data$reactionData[which(!occlusion2.data$reactionData$ID %in% occlusion2.td.highAQ),]
occlusion2.data$descriptives.initial<-occlusion2.data$descriptives[which(!occlusion2.data$descriptives$ID %in% occlusion2.td.highAQ),]

# Who passed criteria for ADOS and AQ? ASD C = 2, S = 4, T = 7; Autism C = 3, S = 6, T = 10
occlusion2.data$descriptives$ADOSclassification <- NA
occlusion2.data$descriptives[which(occlusion2.data$descriptives$ADOScomm >= 2 & occlusion2.data$descriptives$ADOSsocial >= 4 & 
                                     (occlusion2.data$descriptives$ADOScomm + occlusion2.data$descriptives$ADOSsocial) >=7),
                             "ADOSclassification"] <- "asd"
occlusion2.data$descriptives[which(occlusion2.data$descriptives$ADOScomm >= 3 & occlusion2.data$descriptives$ADOSsocial >= 6 & 
                                     (occlusion2.data$descriptives$ADOScomm + occlusion2.data$descriptives$ADOSsocial) >=10),
                             "ADOSclassification"] <- "autism"

# Remove participants who could not do the task
occlusion2.unable <- data.frame(ID=c("A05", "A07", "A19", "A22", "A26"),
           group=NA, 
           stringsAsFactors=FALSE)
occlusion2.unable[grep("^A", occlusion2.unable$ID), ]$group<-"asd"
occlusion2.data$psychometric.parameters.trunc<-occlusion2.data$psychometric.parameters.trunc[which(!occlusion2.data$psychometric.parameters.trunc$ID %in% occlusion2.unable$ID),]
occlusion2.data$reactionData.trunc<-occlusion2.data$reactionData.trunc[which(!occlusion2.data$reactionData.trunc$ID %in% occlusion2.unable$ID),]
occlusion2.data$descriptives.trunc<-occlusion2.data$descriptives.initial[which(!occlusion2.data$descriptives.initial$ID %in% occlusion2.unable$ID),]

# Make it easier to analyse reaction time
Renaming<-function(df){
  if(length(df$intensity) == 5){
    df$intensity<-c(-2.5,-1,0,1,2.5)
  } else if(length(df$intensity) == 3){
    df$intensity<-c(-1,0,1)
  }
  
  df
}

occlusion2.data$reactionData.levels.trunc<-ddply(occlusion2.data$reactionData.trunc,.(disparitySign, condition, ID), Renaming)
occlusion2.data$reactionData.levels.trunc$intensity<-as.factor(occlusion2.data$reactionData.levels.trunc$intensity)

occlusion2.data$psychometric.parameters.trunc.difference<-subset(occlusion2.data$psychometric.parameters.trunc, condition != "baseline", select = c("ID", "disparitySign", "condition", "threshold", "slope", "group"))
baseline.psychometric <-subset(occlusion2.data$psychometric.parameters.trunc,condition == "baseline", select = c("ID", "disparitySign", "threshold"))
colnames(baseline.psychometric)<-c("ID", "disparitySign", "baseline")
occlusion2.data$psychometric.parameters.trunc.difference<-merge(occlusion2.data$psychometric.parameters.trunc.difference, baseline.psychometric, by = c("ID", "disparitySign"))
baseline.psychometric <- NULL
occlusion2.data$psychometric.parameters.trunc.difference$difference<-occlusion2.data$psychometric.parameters.trunc.difference$baseline - occlusion2.data$psychometric.parameters.trunc.difference$threshold

occlusion2.data$reactionData.levels.trunc.difference<-subset(occlusion2.data$reactionData.levels.trunc, condition != "baseline", select = c("ID", "intensity", "disparitySign", "condition", "rt", "sd", "rtReciprocal", "sdReciprocal", "mu", "sigma", "tau", "group"))

covariates <- subset(occlusion2.data$psychometric.parameters.trunc.difference,,select = c("ID", "disparitySign", "baseline"))[!duplicated(subset(occlusion2.data$psychometric.parameters.trunc.difference,,select = c("ID", "disparitySign", "baseline"))),]
covariates <- dcast(covariates, ...~disparitySign, value.var = "baseline")
colnames(covariates) <- c("ID", "crossedBaseline", "uncrossedBaseline")

occlusion2.data$psychometric.parameters.trunc.difference<-merge(occlusion2.data$psychometric.parameters.trunc.difference, covariates, by = c("ID"))
occlusion2.data$reactionData.levels.trunc.difference<-merge(occlusion2.data$reactionData.levels.trunc.difference, covariates, by = c("ID"))

covariates <- NULL

occlusion2.data$psychometric.parameters.trunc.difference$condition <- as.factor(as.character(occlusion2.data$psychometric.parameters.trunc.difference$condition))
occlusion2.data$reactionData.levels.trunc.difference$condition <- as.factor(as.character(occlusion2.data$reactionData.levels.trunc.difference$condition))

occlusion2.data$psychometric.parameters.outliersRemoved<-ddply(occlusion2.data$psychometric.parameters.trunc.difference, .(disparitySign, condition), function(d){
  limitsThreshold = median(log10(d$threshold)) + 2.5*c(-1, 1)*mad(log10(d$threshold))
  limitsSlope = median(d$slope) + 2.5*c(-1, 1)*mad(d$slope)
  d$threshold[which(((log10(d$threshold) - limitsThreshold[1])*(limitsThreshold[2] - log10(d$threshold))) <= 0)]<-NA  
  d$slope[which(((d$slope - limitsSlope[1])*(limitsSlope[2] - d$slope)) <= 0)]<-NA
  return(d)
})

occlusion2.data$reactionData.levels.outliersRemoved<-ddply(occlusion2.data$reactionData.levels.trunc.difference, .(disparitySign, condition,intensity), function(d){
  limitsRT = median(d$rt) + 2.5*c(-1, 1)*mad(d$rt)
  limitsSD = median(d$sd) + 2.5*c(-1, 1)*mad(d$sd)
  limitsRTReciprocal = median(d$rtReciprocal) + 2.5*c(-1, 1)*mad(d$rtReciprocal)
  limitsSDReciprocal = median(d$sdReciprocal) + 2.5*c(-1, 1)*mad(d$sdReciprocal)
  limitsMu = median(d$mu) + 2.5*c(-1, 1)*mad(d$mu)
  limitsSigma = median(d$sigma) + 2.5*c(-1, 1)*mad(d$sigma)
  limitsTau = median(d$tau) + 2.5*c(-1, 1)*mad(d$tau)
  d$rt[which(((d$rt - limitsRT[1])*(limitsRT[2] - d$rt)) <= 0)]<-NA  
  d$sd[which(((d$sd - limitsSD[1])*(limitsSD[2] - d$sd)) <= 0)]<-NA  
  d$rtReciprocal[which(((d$rtReciprocal - limitsRTReciprocal[1])*(limitsRTReciprocal[2] - d$rtReciprocal)) <= 0)]<-NA  
  d$sdReciprocal[which(((d$sdReciprocal - limitsSDReciprocal[1])*(limitsSDReciprocal[2] - d$sdReciprocal)) <= 0)]<-NA  
  d$mu[which(((d$mu - limitsMu[1])*(limitsMu[2] - d$mu)) <= 0)]<-NA  
  d$sigma[which(((d$sigma - limitsSigma[1])*(limitsSigma[2] - d$sigma)) <= 0)]<-NA  
  d$tau[which(((d$tau - limitsTau[1])*(limitsTau[2] - d$tau)) <= 0)]<-NA
  return(d)
})

occlusion2.data$psychometric.parameters.outliersRemoved$c.log.crossedBaseline<-scale(log10(occlusion2.data$psychometric.parameters.outliersRemoved$crossedBaseline), scale = F, center = T)
occlusion2.data$psychometric.parameters.outliersRemoved$c.log.uncrossedBaseline<-scale(log10(occlusion2.data$psychometric.parameters.outliersRemoved$uncrossedBaseline), scale = F, center = T)
occlusion2.data$reactionData.levels.outliersRemoved$c.log.crossedBaseline<-scale(log10(occlusion2.data$reactionData.levels.outliersRemoved$crossedBaseline), scale = F, center = T)
occlusion2.data$reactionData.levels.outliersRemoved$c.log.uncrossedBaseline<-scale(log10(occlusion2.data$reactionData.levels.outliersRemoved$uncrossedBaseline), scale = F, center = T)

occlusion2.data$psychometric.parameters.outliersRemoved$group<- factor(occlusion2.data$psychometric.parameters.outliersRemoved$group, levels  = c("td", "asd"))
occlusion2.data$psychometric.parameters.outliersRemoved$disparitySign<- factor(occlusion2.data$psychometric.parameters.outliersRemoved$disparitySign, levels  = c("crossed", "uncrossed"))
occlusion2.data$psychometric.parameters.outliersRemoved$condition<- factor(occlusion2.data$psychometric.parameters.outliersRemoved$condition, levels  = c("targetlowerlimit","occlusion", "conflict"))
occlusion2.data$reactionData.levels.outliersRemoved$group<- factor(occlusion2.data$reactionData.levels.outliersRemoved$group, levels  = c("td", "asd"))
occlusion2.data$reactionData.levels.outliersRemoved$disparitySign<- factor(occlusion2.data$reactionData.levels.outliersRemoved$disparitySign, levels  = c("crossed", "uncrossed"))
occlusion2.data$reactionData.levels.outliersRemoved$condition<- factor(occlusion2.data$reactionData.levels.outliersRemoved$condition, levels  = c("targetlowerlimit","occlusion", "conflict"))

## ---- occlusion2-qqplot ----
occlusion2.data$qqplot <- melt(dcast(occlusion2.data$psychometric.parameters.outliersRemoved[,c("ID", "disparitySign", 
                                                                                                "condition", "threshold", "baseline")], 
                                     ID + disparitySign + baseline ~ ..., value.var = "threshold"), 
                               id=c("ID", "disparitySign"),variable.name = "condition", value.name = "threshold")

occlusion2.data$qqplot$log.threshold <- log10(occlusion2.data$qqplot$threshold)
names(occlusion2.data$qqplot) <- c("ID", "disparitySign", "condition" ,"untransformed", "log-transformed")

levels(occlusion2.data$qqplot$condition) <- c("baseline", "beside", "occlusion", "conflict")

occlusion2.data$qqplot.comparison <- melt(occlusion2.data$qqplot, id=c("ID", "disparitySign", "condition"), variable.name = "transformed", value.name = "threshold")
occlusion2.data$qqplot.comparison.stats <- ddply(occlusion2.data$qqplot.comparison, .(condition, transformed), summarize,
                                                 skew = format(round(skewness(threshold, na.rm = TRUE), 3), nsmall = 3),
                                                 kurt = format(round(kurtosis (threshold, na.rm = TRUE), 3), nsmall = 3),
                                                 jb = format(round(jarque.test(threshold[!is.na(threshold)])$statistic, 3), nsmall = 3),
                                                 p.value = format.pval(jarque.test(threshold[!is.na(threshold)])$p.value, digits = 3, eps = 0.001))

occlusion2.data$qqplot.comparison.resid <- ddply(.data = occlusion2.data$qqplot.comparison, .variables = .(condition, transformed), function(dat){
  q <- qqnorm(dat$threshold, plot = FALSE)
  dat$xq <- q$x
  dat
}
)

occlusion2.qqplot <- ggplot(data = occlusion2.data$qqplot.comparison.resid, aes(x = xq, y = threshold)) +
  geom_smooth(method = "lm", se = FALSE, size = 1) +
  geom_point() +
  #   geom_text(data = occlusion2.data$qqplot.comparison.stats, aes(label=paste("JB = ", jb), size = 12),
  #             x=Inf, y = -Inf, , hjust=1.02, vjust=-5.4) +
  #   geom_text(data = occlusion2.data$qqplot.comparison.stats, aes(label=paste("p = ", p.value), size = 12),
  #             x=Inf, y = -Inf, , hjust=1.02, vjust=-4) +
  geom_text(data = occlusion2.data$qqplot.comparison.stats, aes(label=paste("S = ", skew), size = 14),
            x=Inf, y = -Inf, , hjust=1.1, vjust=-2.6) +
  geom_text(data = occlusion2.data$qqplot.comparison.stats, aes(label=paste("C = ", kurt), size = 14),
            x=Inf, y = -Inf, , hjust=1.1, vjust=-1.2) +
  xlab("theoretical quantiles") +
  ylab("sample quantiles") +
  facet_grid(transformed ~ condition, scales = "free") +
  theme_bw(base_size = 14)+
  theme(legend.position="none",
        panel.background = element_blank(), 
        axis.title.y = element_text(vjust=1.5),
        axis.title.x = element_text(vjust=0),
        strip.background = element_rect(fill = "black", color = "black"),
        strip.text = element_text(face = "bold", colour = "white"))

## ---- fig-occlusion2-qqplot ----
occlusion2.qqplot

## ---- occlusion2-descriptives ----
occlusion2.matching<-sapply(c("initial", "final"),function(x) NULL)
occlusion2.matching$initial<-sapply(c("age", "verbal", "perf", "full","aq"),function(x) NULL)
occlusion2.matching$initial$age<-t.test(age~group, occlusion2.data$descriptives.initial)
occlusion2.matching$initial$verbal<-t.test(verbalRaw~group, occlusion2.data$descriptives.initial)
occlusion2.matching$initial$perf<-t.test(perfRaw~group, occlusion2.data$descriptives.initial)
occlusion2.matching$initial$full<-t.test(fullRaw~group, occlusion2.data$descriptives.initial)
occlusion2.matching$initial$aq<-t.test(AQ~group, occlusion2.data$descriptives.initial, na.rm=T)

occlusion2.matching$final<-sapply(c("age", "verbal", "perf", "full","aq"),function(x) NULL)
occlusion2.matching$final$age<-t.test(age~group, occlusion2.data$descriptives.trunc)
occlusion2.matching$final$verbal<-t.test(verbalRaw~group, occlusion2.data$descriptives.trunc)
occlusion2.matching$final$perf<-t.test(perfRaw~group, occlusion2.data$descriptives.trunc)
occlusion2.matching$final$full<-t.test(fullRaw~group, occlusion2.data$descriptives.trunc)
occlusion2.matching$final$aq<-t.test(AQ~group, occlusion2.data$descriptives.trunc, na.rm=T)

## ---- occlusion2-baseline-stereo-ttest ----
occlusion2.baseline.comparisons<-sapply(c("tno", "crossed", "uncrossed", "asd", "td", "correlations"),function(x) NULL)
occlusion2.data$baseline <- merge(occlusion2.data$descriptives.trunc[,c(1,12:13)], unique(occlusion2.data$psychometric.parameters.outliersRemoved[,c(1,6,9:10)]), by = c("ID", "group"))
occlusion2.data$baseline$TNOScore <- occlusion2.data$baseline$TNOScore/3600
occlusion2.data$baseline <- melt(occlusion2.data$baseline, id = c("ID", "group"), variable.name = "measure", value.name = "threshold")

occlusion2.baseline.comparisons$tno<-t.test(log10(TNOScore/3600)~group, data = occlusion2.data$descriptives.trunc, alternative = "greater")
occlusion2.baseline.comparisons$crossed<-t.test(log10(baseline)~group, data = unique(subset(occlusion2.data$psychometric.parameters.outliersRemoved, disparitySign == "crossed" , select = c('ID', 'disparitySign', 'baseline', 'group'))), alternative = "less")
occlusion2.baseline.comparisons$uncrossed<-t.test(log10(baseline)~group, data = unique(subset(occlusion2.data$psychometric.parameters.outliersRemoved, disparitySign == "uncrossed" , select = c('ID', 'disparitySign', 'baseline', 'group'))), alternative = "less")
occlusion2.baseline.comparisons$asd<-t.test(log10(baseline)~disparitySign, data = unique(subset(occlusion2.data$psychometric.parameters.outliersRemoved, group == "asd" , select = c('ID', 'disparitySign', 'baseline', 'group'))), paired = TRUE, var.equal = TRUE)
occlusion2.baseline.comparisons$td<-t.test(log10(baseline)~disparitySign, data = unique(subset(occlusion2.data$psychometric.parameters.outliersRemoved, group == "td" , select = c('ID', 'disparitySign', 'baseline', 'group'))), paired = TRUE, var.equal = TRUE)
occlusion2.baseline.comparisons$correlations <- corstars(dcast(occlusion2.data$baseline, ID+group~measure)[,c("TNOScore", "crossedBaseline", "uncrossedBaseline")])
occlusion2.baseline.comparisons$tno$descriptives<-data.table(occlusion2.data$descriptives.trunc)[,list(mean=mean(log10(TNOScore/3600), na.rm=T), sd=sd(log10(TNOScore/3600), na.rm=T)), by=group]
occlusion2.baseline.comparisons$uncrossed$descriptives<-data.table(unique(subset(occlusion2.data$psychometric.parameters.outliersRemoved, disparitySign == "uncrossed" , select = c('ID', 'disparitySign', 'baseline', 'group'))))[,list(mean=mean(log10(baseline), na.rm=T), sd=sd(log10(baseline), na.rm=T)), by=group]

## ---- occlusion2-mixed-models ----
# This model fails to converge, but max(abs(relgrad)) is < 0.001 (0.0000562) so it's okay -- see https://github.com/lme4/lme4/issues/120
occlusion2.mixed.threshold<-mixed(log10(threshold)~group*condition*disparitySign*c.log.crossedBaseline + group*condition*disparitySign*c.log.uncrossedBaseline+(1+disparitySign+condition|ID),occlusion2.data$psychometric.parameters.outliersRemoved, control = lmerControl(optCtrl = list(maxfun = 100000)))
occlusion2.mixed.slope<-mixed(slope~group*condition*disparitySign*c.log.crossedBaseline + group*condition*disparitySign*c.log.uncrossedBaseline+(1+disparitySign+condition|ID),occlusion2.data$psychometric.parameters.outliersRemoved, control = lmerControl(optCtrl = list(maxfun = 100000)))

occlusion2.mixed.rt<-mixed(rtReciprocal~group*condition*disparitySign*c.log.crossedBaseline + group*condition*disparitySign*c.log.uncrossedBaseline+(1+(disparitySign*condition)|ID),occlusion2.data$reactionData.levels.outliersRemoved, control = lmerControl(optCtrl = list(maxfun = 100000)))
occlusion2.mixed.sd<-mixed(sdReciprocal~group*condition*disparitySign*c.log.crossedBaseline + group*condition*disparitySign*c.log.uncrossedBaseline+(1+(disparitySign*condition)|ID),occlusion2.data$reactionData.levels.outliersRemoved, control = lmerControl(optCtrl = list(maxfun = 100000)))

## ---- occlusion2-mixed-tables ----
tableRowRename <- function(table){
  table$Effect<-gsub("c.log.crossedBaseline","crossedBaseline", table$Effect)
  table$Effect<-gsub("c.log.uncrossedBaseline","uncrossedBaseline", table$Effect)
  return(table)
}

occlusion2.mixed.threshold.tables <- afexMixedTables(occlusion2.mixed.threshold, adj.p = T, adj.p.n = 4)
occlusion2.mixed.threshold.tables$lme.formatted <- tableRowRename(occlusion2.mixed.threshold.tables$lme)
colnames(occlusion2.mixed.threshold.tables$lme.formatted) <- c("Effect", "$\\beta$", "B", "B Std. Error", "t value")
occlusion2.mixed.threshold.tables$lme.latex <- xtable(occlusion2.mixed.threshold.tables$lme.formatted, 
                                                 label="table:table-occlusion2-lmm-threshold", 
                                                 caption = c("Raw output from linear mixed-model analysis of 
                                          participant characteristics (including crossed and uncrossed baseline
                                          stereoacuity and presence of autism diagnosis), disparity sign, and occlusion configuration,
                                                  and the effect of these predictors upon relative disparity threshold.",
                                                  "Raw output from linear mixed-model analysis of relative disparity threshold"), 
                                                 digits = 3, align="llccrr")
occlusion2.mixed.threshold.tables$formatted <- tableRowRename(occlusion2.mixed.threshold.tables$formatted)
occlusion2.mixed.threshold.tables$latex <- xtable(occlusion2.mixed.threshold.tables$formatted, 
                                                  label="table:table-occlusion2-mixed-threshold", 
                                                  caption = c("F-test of fixed terms in a linear mixed model 
                                                  analysis of participant characteristics (including baseline stereoacuity 
                                                  and presence of \\gls{asd} diagnosis), occlusion cues, and sign of disparity,
                                                  and the effect of these predictors upon relative disparity threshold.",
                                                  "F-test of linear mixed-model analysis of relative disparity threshold"), 
                                                  digits = 3, align="llccrr")

occlusion2.mixed.slope.tables <- afexMixedTables(occlusion2.mixed.slope, adj.p = T, adj.p.n = 4)
occlusion2.mixed.slope.tables$lme.formatted <- tableRowRename(occlusion2.mixed.slope.tables$lme)
colnames(occlusion2.mixed.slope.tables$lme.formatted) <- c("Effect", "$\\beta$", "B", "B Std. Error", "t value")
occlusion2.mixed.slope.tables$lme.latex <- xtable(occlusion2.mixed.slope.tables$lme.formatted, 
                                                      label="table:table-occlusion2-lmm-slope", 
                                                      caption = c("Raw output from linear mixed-model analysis of 
                                          participant characteristics (including crossed and uncrossed baseline
                                          stereoacuity and presence of autism diagnosis), disparity sign, and occlusion configuration,
                                                  and the effect of these predictors upon the slope of the psychometric function.",
                                                  "Raw output from linear mixed-model analysis of psychometric function slope"), 
                                                      digits = 3, align="llccrr")
occlusion2.mixed.slope.tables$formatted <- tableRowRename(occlusion2.mixed.slope.tables$formatted)
occlusion2.mixed.slope.tables$latex <- xtable(occlusion2.mixed.slope.tables$formatted, 
                                              label="table:table-occlusion2-mixed-slope",
                                              caption = c("F-test of fixed terms in a linear mixed model 
                                              for predicting the slope of the psychometric function at
                                              threshold.",
                                              "F-test of linear mixed-model analysis of psychometric function slope"), 
                                              digits = 3, align="llccrr")

occlusion2.mixed.rt.tables <- afexMixedTables(occlusion2.mixed.rt, adj.p = T, adj.p.n = 4)
occlusion2.mixed.rt.tables$lme.formatted <- tableRowRename(occlusion2.mixed.rt.tables$lme)
colnames(occlusion2.mixed.rt.tables$lme.formatted) <- c("Effect", "$\\beta$", "B", "B Std. Error", "t value")
occlusion2.mixed.rt.tables$lme.latex <- xtable(occlusion2.mixed.rt.tables$lme.formatted, 
                                                      label="table:table-occlusion2-lmm-rt", 
                                                      caption = c("Raw output from linear mixed-model analysis of 
                                                      participant characteristics (including crossed and uncrossed baseline
                                                      stereoacuity and presence of autism diagnosis), disparity sign, and occlusion configuration,
                                                      and the effect of these predictors upon median reaction speed to
                                                      relative disparity.",
                                                      "Raw output from linear mixed-model analysis of reaction speed"), 
                                                      digits = 3, align="llccrr")
occlusion2.mixed.rt.tables$formatted <- tableRowRename(occlusion2.mixed.rt.tables$formatted)
occlusion2.mixed.rt.tables$latex <- xtable(occlusion2.mixed.rt.tables$formatted, 
                                           label="table:table-occlusion2-mixed-rt", 
                                           caption = c("F-test of fixed terms in a linear mixed model 
                                           for predicting speed of response to relative disparity.",
                                           "F-test of linear mixed-model analysis of reaction speed"), 
                                           digits = 3, align="llccrr")

occlusion2.mixed.sd.tables <- afexMixedTables(occlusion2.mixed.sd, adj.p = T, adj.p.n = 4)
occlusion2.mixed.sd.tables$lme.formatted <- tableRowRename(occlusion2.mixed.sd.tables$lme)
colnames(occlusion2.mixed.sd.tables$lme.formatted) <- c("Effect", "$\\beta$", "B", "B Std. Error", "t value")
occlusion2.mixed.sd.tables$lme.latex <- xtable(occlusion2.mixed.sd.tables$lme.formatted, 
                                               label="table:table-occlusion2-lmm-sd", 
                                               caption = c("Raw output from linear mixed-model analysis of 
                                                      pasdicipant characteristics (including crossed and uncrossed baseline
                                                      stereoacuity and presence of autism diagnosis), disparity sign, and occlusion configuration,
                                                      and the effect of these predictors upon reaction speed variability to
                                                      relative disparity.",
                                                      "Raw output from linear mixed-model analysis of reaction speed variability"), 
                                               digits = 3, align="llccrr")
occlusion2.mixed.sd.tables$formatted <- tableRowRename(occlusion2.mixed.sd.tables$formatted)
occlusion2.mixed.sd.tables$latex <- xtable(occlusion2.mixed.sd.tables$formatted, 
                                           label="table:table-occlusion2-mixed-sd", 
                                           caption = c("F-test of fixed terms in a linear mixed model 
                                           for predicting the variability of speed of response to relative disparity.",
                                           "F-test of linear mixed-model analysis of reaction speed variability"), 
                                           digits = 3, align="llccrr")

## ---- table-occlusion2-mixed-threshold ----
print(occlusion2.mixed.threshold.tables$latex, caption.placement = "top", only.contents = T, 
             include.colnames = FALSE, include.rownames = FALSE, 
             booktabs = TRUE, sanitize.text.function = I, tabular.environment="longtable",
             hline.after=c(-1,nrow(occlusion2.mixed.threshold.tables$formatted)),
             add.to.row=list(pos=list(0),command="\\midrule \\endhead "))

## ---- table-occlusion2-mixed-slope ----
print(occlusion2.mixed.slope.tables$latex, caption.placement = "top", only.contents = T, 
             include.colnames = FALSE, include.rownames = FALSE, 
             booktabs = TRUE, sanitize.text.function = I, tabular.environment="longtable",
             hline.after=c(-1,nrow(occlusion2.mixed.slope.tables$formatted)),
             add.to.row=list(pos=list(0),command="\\midrule \\endhead "))

## ---- table-occlusion2-mixed-rt ----
print(occlusion2.mixed.rt.tables$latex, caption.placement = "top", only.contents = T, 
             include.colnames = FALSE, include.rownames = FALSE, 
             booktabs = TRUE, sanitize.text.function = I, tabular.environment="longtable",
             hline.after=c(-1,nrow(occlusion2.mixed.rt.tables$formatted)),
             add.to.row=list(pos=list(0),command="\\midrule \\endhead "))

## ---- table-occlusion2-mixed-sd ----
print(occlusion2.mixed.sd.tables$latex, caption.placement = "top", only.contents = T, 
             include.colnames = FALSE, include.rownames = FALSE, 
             booktabs = TRUE, sanitize.text.function = I, tabular.environment="longtable",
             hline.after=c(-1,nrow(occlusion2.mixed.sd.tables$formatted)),
             add.to.row=list(pos=list(0),command="\\midrule \\endhead "))

## ---- table-occlusion2-lmm-threshold ----
print(occlusion2.mixed.threshold.tables$lme.latex, caption.placement = "top", only.contents = T, 
             include.colnames = FALSE, include.rownames = FALSE, 
             booktabs = TRUE, sanitize.text.function = I, tabular.environment="longtable",
             hline.after=c(-1,nrow(occlusion2.mixed.threshold.tables$lme)),
             add.to.row=list(pos=list(0),command="\\midrule \\endhead "))

## ---- table-occlusion2-lmm-slope ----
print(occlusion2.mixed.slope.tables$lme.latex, caption.placement = "top", only.contents = T, 
             include.colnames = FALSE, include.rownames = FALSE, 
             booktabs = TRUE, sanitize.text.function = I, tabular.environment="longtable",
             hline.after=c(-1,nrow(occlusion2.mixed.slope.tables$lme)),
             add.to.row=list(pos=list(0),command="\\midrule \\endhead "))

## ---- table-occlusion2-lmm-rt ----
print(occlusion2.mixed.rt.tables$lme.latex, caption.placement = "top", only.contents = T, 
             include.colnames = FALSE, include.rownames = FALSE, 
             booktabs = TRUE, sanitize.text.function = I, tabular.environment="longtable",
             hline.after=c(-1,nrow(occlusion2.mixed.rt.tables$lme)),
             add.to.row=list(pos=list(0),command="\\midrule \\endhead "))

## ---- table-occlusion2-lmm-sd ----
print(occlusion2.mixed.sd.tables$lme.latex, caption.placement = "top", only.contents = T, 
             include.colnames = FALSE, include.rownames = FALSE, 
             booktabs = TRUE, sanitize.text.function = I, tabular.environment="longtable",
             hline.after=c(-1,nrow(occlusion2.mixed.sd.tables$lme)),
             add.to.row=list(pos=list(0),command="\\midrule \\endhead "))

## ---- occlusion2-posthoc ----
occlusion2.posthoc.threshold.condition <- lsmeans(occlusion2.mixed.threshold$full.model, pairwise~condition)
occlusion2.posthoc.threshold.uncrossedBaseline <- lsmeans(occlusion2.mixed.threshold$full.model, pairwise~c.log.uncrossedBaseline, 
        at = list("c.log.uncrossedBaseline" = c(quantile(occlusion2.data$psychometric.parameters.outliersRemoved$c.log.uncrossedBaseline)[[2]],
                                                quantile(occlusion2.data$psychometric.parameters.outliersRemoved$c.log.uncrossedBaseline)[[4]])))

occlusion2.posthoc.slope.condition <- lsmeans(occlusion2.mixed.slope$full.model, pairwise~condition)
occlusion2.posthoc.slope.crossedBaseline <- lsmeans(occlusion2.mixed.slope$full.model, pairwise~c.log.crossedBaseline, 
                                                      at = list("c.log.crossedBaseline" = c(quantile(occlusion2.data$psychometric.parameters.outliersRemoved$c.log.crossedBaseline)[[2]],
                                                                                              quantile(occlusion2.data$psychometric.parameters.outliersRemoved$c.log.crossedBaseline)[[4]])))
occlusion2.posthoc.slope.uncrossedBaseline <- lsmeans(occlusion2.mixed.slope$full.model, pairwise~c.log.uncrossedBaseline, 
                                                          at = list("c.log.uncrossedBaseline" = c(quantile(occlusion2.data$psychometric.parameters.outliersRemoved$c.log.uncrossedBaseline)[[2]],
                                                                                                  quantile(occlusion2.data$psychometric.parameters.outliersRemoved$c.log.uncrossedBaseline)[[4]])))
occlusion2.posthoc.slope.condition.disparitySign <- lsmeans(occlusion2.mixed.slope$full.model, pairwise~disparitySign|condition)

occlusion2.posthoc.rt.group <- lsmeans(occlusion2.mixed.rt$full.model, pairwise~group)
occlusion2.posthoc.rt.condition <- lsmeans(occlusion2.mixed.rt$full.model, pairwise~condition)
occlusion2.posthoc.rt.uncrossedBaseline <- lsmeans(occlusion2.mixed.rt$full.model, pairwise~c.log.uncrossedBaseline, 
                                                      at = list("c.log.uncrossedBaseline" = c(quantile(occlusion2.data$psychometric.parameters.outliersRemoved$c.log.uncrossedBaseline)[[2]],
                                                                                              quantile(occlusion2.data$psychometric.parameters.outliersRemoved$c.log.uncrossedBaseline)[[4]])))
occlusion2.posthoc.rt.group.disparitySign <- lsmeans(occlusion2.mixed.rt$full.model, pairwise~group|disparitySign)

occlusion2.posthoc.sd.condition <- lsmeans(occlusion2.mixed.sd$full.model, pairwise~condition)
occlusion2.posthoc.sd.uncrossedBaseline <- lsmeans(occlusion2.mixed.sd$full.model, pairwise~c.log.uncrossedBaseline, 
                                                   at = list("c.log.uncrossedBaseline" = c(quantile(occlusion2.data$psychometric.parameters.outliersRemoved$c.log.uncrossedBaseline)[[2]],
                                                                                           quantile(occlusion2.data$psychometric.parameters.outliersRemoved$c.log.uncrossedBaseline)[[4]])))
occlusion2.posthoc.sd.group.condition <- lsmeans(occlusion2.mixed.sd$full.model, pairwise~condition|group)
occlusion2.posthoc.sd.disparitySign.uncrossedBaseline <- lstrends(occlusion2.mixed.sd$full.model, pairwise~disparitySign, var = "c.log.uncrossedBaseline")
## ---- occlusion2-graphs ----
# function for mean labels
mean.n <- function(x){
  return(c(y = median(x)*0.97, label = round(mean(x),2))) 
  # experiment with the multiplier to find the perfect position
}

occlusion2.plot.baseline.comparisons.group <- ggplot(data = occlusion2.data$baseline, 
                                                    aes(x = measure, y = log10(threshold), color = group, fill = group))+
  geom_point(position=position_jitterdodge(dodge.width=.9, jitter.width = 0.4, jitter.height = 0.2), alpha = .4, 
             size = 3) +
  stat_summary(aes(fill = NULL), fun.data = "mean_cl_boot", geom = "crossbar", size = 1, fatten = 2, position = "dodge")+
  
  geom_text(data = ddply(occlusion2.data$baseline, .(group, measure), function(x) smean.cl.boot(log10(x$threshold))),
            aes(label = round(Mean,2), x = measure, y = Upper, color = group), 
            position = position_dodge(width = .9), vjust = -.8, fontface = "bold", show_guide = F)+
  scale_x_discrete(name = "measure", breaks=c("TNOScore", "crossedBaseline", "uncrossedBaseline"), 
                   labels=c("TNO", "crossed", "uncrossed"))+
  ylab(expression(disparity ~ threshold ~ (log[10] ~ degrees)))+
  scale_fill_brewer(palette="Set2")+
  scale_color_brewer(palette="Set2")+
  theme_bw(base_size = 14)+
  theme(legend.justification=c(1,1), 
        legend.position=c(1,1),
        panel.background = element_blank(), 
        panel.grid.minor = element_blank(), 
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.title.y = element_text(vjust=1.5),
axis.title.x = element_text(vjust=0)) 

occlusion2.plot.threshold.condition.group <- ggplot(NULL, aes(x = condition, y = NULL))+
  geom_jitter(data = occlusion2.data$psychometric.parameters.outliersRemoved, aes(x = condition, y = log10(threshold)),
              position = position_jitter(w = .2, h = .2),
              alpha = .4, size = 3)+
  geom_crossbar(data = data.frame(summary(lsmeans(occlusion2.mixed.threshold$full.model, ~ condition|group))),
                aes(y = lsmean, ymax = upper.CL, ymin = lower.CL),
                size = 1, position = "dodge", fatten = 2)+
  facet_grid(.~group, scales = "free", space = "free")+
  scale_x_discrete(name = "occlusion status", breaks=c("occlusion", "conflict", "targetlowerlimit"), labels=c("congruent", "conflict", "beside"))+
  ylab(expression(disparity ~ threshold ~ (log[10] ~ degrees)))+
  geom_text(data = data.frame(condition = "conflict", y = Inf,
                              lab = sprintf("% 3s", c("TD", "ASD")),
                              group = as.character(levels(occlusion2.data$psychometric.parameters.outliersRemoved$group))),
            aes(y = y, label = lab), hjust = .8,
            vjust = 1.2, size = 16, color = "gray80",
            fontface = "bold") + theme_bw(base_size = 14)+
  geom_text(data = data.frame(summary(lsmeans(occlusion2.mixed.threshold$full.model, ~ condition|group))), aes(label = round(lsmean,2), x = condition, y = upper.CL), position = position_dodge(width = .9), vjust = -.8, fontface = "bold", show_guide = F)+
  theme(legend.position = "none",
        panel.background = element_blank(), 
        panel.grid.minor = element_blank(), 
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.title.y = element_text(vjust=1.5),
        axis.title.x = element_text(vjust=0))

occlusion2.plot.threshold.uncrossedBaseline <- ggplot(cbind(occlusion2.data$psychometric.parameters.outliersRemoved, data.frame(facet = c("fake"))),
                                                      aes(x = c.log.uncrossedBaseline, y = log10(threshold))) + geom_point() + 
  geom_abline(intercept = occlusion2.mixed.threshold.tables$lme[which(occlusion2.mixed.threshold.tables$lme$Effect == "(Intercept)"), "B"],
              slope = occlusion2.mixed.threshold.tables$lme[which(occlusion2.mixed.threshold.tables$lme$Effect == "c.log.uncrossedBaseline"), "B"],
              size = 1, color = "red") +
  ylab(expression(disparity ~ threshold ~ (log[10] ~ degrees))) +
  xlab(expression (atop("uncrossed baseline", paste("(log" [10], " degrees)")))) +
  facet_grid(.~facet) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        #         axis.text.y=element_blank(),
        #         axis.ticks.y=element_blank(),
        #         axis.title.y=element_blank(),
        panel.background = element_blank(), 
        panel.grid.minor = element_blank(), 
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.title.y = element_text(vjust=1.5),
        axis.title.x = element_text(vjust=0))

occlusion2.plot.slope.condition.disparitySign.group <- ggplot(NULL, aes(x = condition, y = NULL))+
  geom_point(data = occlusion2.data$psychometric.parameters.outliersRemoved, 
             aes(y = slope, fill = disparitySign, color = disparitySign),
             position = position_jitterdodge(jitter.width = .2, jitter.height = .2, dodge.width = .9),
             alpha = .25, size = 3)+
  geom_crossbar(data = data.frame(summary(lsmeans(occlusion2.mixed.slope$full.model, ~ condition|disparitySign*group))),
                aes(y = lsmean, ymax = upper.CL, ymin = lower.CL, color = disparitySign),
                size = 1, position = "dodge", fatten = 2)+
  facet_grid(.~group, scales = "free", space = "free")+
  scale_x_discrete(name = "occlusion status", breaks=c("occlusion", "conflict", "targetlowerlimit"), labels=c("congruent", "conflict", "beside"))+
  scale_fill_discrete(name="disparity sign")+
  scale_color_discrete(name="disparity sign")+
  ylab("slope of psychometric function")+
  geom_text(data = data.frame(condition = "conflict", y = Inf,
                              lab = sprintf("% 3s", c("TD", "ASD")),
                              group = as.character(levels(occlusion2.data$psychometric.parameters.outliersRemoved$group))),
            aes(y = y, label = lab), hjust = .6,
            vjust = 1.2, size = 18, color = "gray80",
            fontface = "bold") + theme_bw(base_size = 14)+
  geom_text(data = data.frame(summary(lsmeans(occlusion2.mixed.slope$full.model, ~ condition|disparitySign*group))), aes(label = round(lsmean,2), x = condition, y = upper.CL, color = disparitySign), position = position_dodge(width = .9), vjust = -.8, fontface = "bold", show_guide = F)+
  theme(legend.justification=c(0,1), 
        legend.position=c(0,1),
        panel.background = element_blank(), 
        panel.grid.minor = element_blank(), 
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.title.y = element_text(vjust=1.5),
        axis.title.x = element_text(vjust=0))

## Not run
# occlusion2.plot.rt.condition.disparitySign.group <- ggplot(NULL, aes(x = condition, y = NULL))+
#   geom_point(data = occlusion2.data$reactionData.levels.outliersRemoved, 
#              aes(y = rtReciprocal, fill = disparitySign, color = disparitySign),
#              position = position_jitterdodge(jitter.width = .2, jitter.height = .2, dodge.width = .9),
#              alpha = .25, size = 3)+
#   geom_crossbar(data = data.frame(summary(lsmeans(occlusion2.mixed.rt$full.model, ~ condition|disparitySign*group))),
#                 aes(y = lsmean, ymax = upper.CL, ymin = lower.CL, color = disparitySign),
#                 size = 1, position = "dodge", fatten = 2)+
#   facet_grid(.~group, scales = "free", space = "free")+
#   scale_x_discrete(name = "occlusion status", breaks=c("occlusion", "conflict", "targetlowerlimit"), labels=c("congruent", "conflict", "beside"))+
#   scale_fill_discrete(name="disparity sign")+
#   scale_color_discrete(name="disparity sign")+
#   ylab("median reaction speed")+
#   geom_text(data = data.frame(condition = "conflict", y = Inf,
#                               lab = sprintf("% 3s", c("TD", "ASD")),
#                               group = as.character(levels(occlusion2.data$psychometric.parameters.outliersRemoved$group))),
#             aes(y = y, label = lab), hjust = .6,
#             vjust = 1.2, size = 18, color = "gray80",
#             fontface = "bold") + theme_bw(base_size = 14)+
#   geom_text(data = data.frame(summary(lsmeans(occlusion2.mixed.rt$full.model, ~ condition|disparitySign*group))), aes(label = round(lsmean,2), x = condition, y = upper.CL, color = disparitySign), position = position_dodge(width = .9), vjust = -.8, fontface = "bold", show_guide = F)+
#   theme(legend.justification=c(0,1), 
#         legend.position=c(0,1),
#         panel.background = element_blank(), 
#         panel.grid.minor = element_blank(), 
#         strip.background = element_blank(),
#         strip.text.x = element_blank(),
#         axis.title.y = element_text(vjust=1.5),
#         axis.title.x = element_text(vjust=0))

occlusion2.table.rt.group.disparitySign<-data.frame(summary(lsmeans(occlusion2.mixed.rt$full.model, ~disparitySign|group)))[,c("disparitySign", "group", "lsmean", "SE")]

occlusion2.table.rt.group.disparitySign<-rbind(occlusion2.table.rt.group.disparitySign, 
                                               cbind(data.frame(summary(lsmeans(
                                                 occlusion2.mixed.rt$full.model, ~disparitySign)))
                                                 [,c("disparitySign", "lsmean", "SE")], 
                                                 group = c("mean")),cbind(data.frame(summary(lsmeans(
                                                   occlusion2.mixed.rt$full.model, ~group)))
                                                   [,c("group", "lsmean", "SE")], 
                                                   disparitySign = c("overall")), 
                                               c("overall", "mean", mean(occlusion2.data$reactionData.levels.outliersRemoved$rtReciprocal, na.rm = T),
                                                 std.error(occlusion2.data$reactionData.levels.outliersRemoved$sdReciprocal, na.rm = T)))

occlusion2.table.rt.group.disparitySign[,c(-1,-2)] <- apply(occlusion2.table.rt.group.disparitySign[,c(-1,-2)], 2, function(x) as.numeric(x));
occlusion2.table.rt.group.disparitySign[,c(-1,-2)] <- round(occlusion2.table.rt.group.disparitySign[,c(-1, -2)], digits = 2)
occlusion2.table.rt.group.disparitySign$SE <- format.pval(occlusion2.table.rt.group.disparitySign$SE, digits = 2, eps = 0.01)
occlusion2.table.rt.group.disparitySign$stats <- apply( occlusion2.table.rt.group.disparitySign[ , c("lsmean", "SE") ] , 1 , paste , collapse = " $\\pm$ " )

occlusion2.table.rt.group.disparitySign <- acast(occlusion2.table.rt.group.disparitySign[,c("disparitySign", "group", "stats")], disparitySign ~ group)

## Not run
# occlusion2.plot.sd.group.condition <- ggplot(NULL, aes(x = condition, y = NULL))+
#   geom_jitter(data = subset(occlusion2.data$reactionData.levels.outliersRemoved, intensity !=0), 
#               aes(x = condition, y = sdReciprocal),
#               position = position_jitter(w = .2, h = .2),
#               alpha = .4, size = 3)+
#   geom_crossbar(data = data.frame(summary(lsmeans(occlusion2.mixed.sd$full.model, ~ condition|group))),
#                 aes(y = lsmean, ymax = upper.CL, ymin = lower.CL),
#                 size = 1, position = "dodge", fatten = 2)+
#   facet_grid(.~group, scales = "free", space = "free")+
#   scale_x_discrete(name = "occlusion status", breaks=c("occlusion", "conflict", "targetlowerlimit"), labels=c("congruent", "conflict", "beside"))+
#   ylab("standard deviation of response speed")+
#   geom_text(data = data.frame(condition = "conflict", y = Inf,
#                               lab = sprintf("% 3s", c("TD", "ASD")),
#                               group = as.character(levels(occlusion2.data$psychometric.parameters.outliersRemoved$group))),
#             aes(y = y, label = lab), hjust = .6,
#             vjust = 1.2, size = 18, color = "gray80",
#             fontface = "bold") + theme_bw(base_size = 14)+
#   theme(legend.position = "none",
#         panel.background = element_blank(), 
#         panel.grid.minor = element_blank(), 
#         strip.background = element_blank(),
#         strip.text.x = element_blank(),
#         axis.title.y = element_text(vjust=1.5),
#         axis.title.x = element_text(vjust=0))

occlusion2.posthoc.sd.group.condition.table <- data.frame(cld(occlusion2.posthoc.sd.group.condition, by = NULL, sort = FALSE, Letters = c("abcdef")))
occlusion2.posthoc.sd.group.condition.table <- occlusion2.posthoc.sd.group.condition.table[,c(2,1,3:8)]
occlusion2.posthoc.sd.group.condition.table$condition<- as.character(unlist(occlusion2.posthoc.sd.group.condition.table$condition))
occlusion2.posthoc.sd.group.condition.table$condition[occlusion2.posthoc.sd.group.condition.table$condition == "targetlowerlimit"]<-"beside"
colnames(occlusion2.posthoc.sd.group.condition.table)<-c("Group", "Condition", "\\gls{lsm}", "\\gls{se}", "df", "lower CL", "upper CL", "Posthoc")
occlusion2.posthoc.sd.group.condition.table$Group<-cleanTableColumns(occlusion2.posthoc.sd.group.condition.table$Group)
occlusion2.posthoc.sd.group.condition.table$Posthoc<-whitespaceLaTeXCompatible(occlusion2.posthoc.sd.group.condition.table$Posthoc)

## ---- fig-occlusion2-baseline-comparisons ----
occlusion2.plot.baseline.comparisons.group

## ---- fig-occlusion2-threshold ----
plot_grid(occlusion2.plot.threshold.condition.group,occlusion2.plot.threshold.uncrossedBaseline, 
          labels = "AUTO", align = 'h', rel_widths = c(1, 0.5))

## ---- fig-occlusion2-slope-condition-disparitySign ----
occlusion2.plot.slope.condition.disparitySign.group

## ---- fig-occlusion2-rt-group-condition-disparitySign ----
# occlusion2.plot.rt.condition.disparitySign.group

## ---- occlusion2-rt-group-disparitySign-descriptives ----
print.xtable(occlusion2.table.rt.group.disparitySign, only.contents = T, booktabs = T, sanitize.text.function = I)

## ---- fig-occlusion2-sd-group-condition ----
# occlusion2.plot.sd.group.condition

## ---- table-occlusion2-sd-group-condition ----
print(xtable(occlusion2.posthoc.sd.group.condition.table, label="table:table-occlusion2-sd-group-condition", caption = c("Pair-wise comparisons of least-squared means for two-way interaction between diagnostic group and occlusion configuration (the latter denoted by \"condition\" in the table), using Tukey's honest significant difference test with $\\alpha = .05$. \\emph{Note}. Rows containing the same letter are not significantly different to each other \\protect\\cite{Piepho:2004ku}. Acronyms: \\acrfull{td}, \\acrfull{asd}, \\acrfull{lsm}, \\acrfull{se}, \\acrfull{cl}.","Pair-wise comparisons for two way interaction between diagnostic group and occlusion configuration for \\acrshort{td} and \\acrshort{asd}"), digits = 3, align="lllcccccc"), include.rownames = FALSE, booktabs = TRUE, caption.placement = "top", sanitize.text.function = I)