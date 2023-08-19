library(readxl)
library(ggplot2)
library(ggbeeswarm)
library(tidyverse, quietly = TRUE)
library(arm)  # st. error
library(mdscore)  # for calculating LRT
library(ROCit)
library(flextable)
library(dplyr)
library(lmtest)
library(pscl)

df<-as.data.frame(read_xlsx("D:/Basel/master/7 Data/parous rate/parous_graph.xlsx"))
df<-as.data.frame(read_xlsx("D:/Basel/master/2 R/parous/geographic/climate.xlsx"))

mda<-mda %>% mutate(insecticide_control = relevel(factor(insecticide_control), "f"),
                    #country = relevel(factor(country), "Thailand"),
                    season = relevel(factor(season), "rainy"),
                    complex = relevel(factor(complex), "Maculatus"),
                    location_combine = relevel(factor(location_combine), "location_combine_indoor"),
                    host_sampling = relevel(factor(host_sampling), "host_sampling_biting_w"),
                    #         zone=relevel(factor(zone), "tro_savannah"),
                    climate=relevel(factor(climate),"tro"))

#reference group
df<-df %>%
  mutate(zone = case_when(zone == 1 ~ "tro_rainforest",
                          zone == 2 ~ "tro_monsoon",
                          zone == 3 ~ "tro_savannah",
                          zone == 11 ~ "tem_dry win_hot sum",
                          zone == 12~ "tem_dry win_warm sum",
                          zone == 14 ~ "tem_no dry",
                          zone == 21~ "cold_dry win_hot sum",
                          zone == 25~ "cold_no dry_hot sum"))


df$country = relevel(factor(df$country), ref = "Thailand")
df$insecticide_control = relevel(factor(df$insecticide_control), ref = "f")
df$season = relevel(factor(df$season), ref = "rainy")
df$complex = relevel(factor(df$complex), ref = "Maculatus")
df$location_combine = relevel(factor(df$location_combine), ref = "indoor")
df$host_sampling = relevel(factor(df$host_sampling), ref = "biting_w")
df$zone=relevel(factor(df$zone), ref = "tro_savannah")
df$climate=relevel(factor(df$climate), ref = "tro")


#A quick (and dirty) way to check for strong correlations between variables
Corr <- cor(sapply(df, as.numeric),
            use = "pairwise.complete.obs", method = "spearman")
corrplot::corrplot(Corr, method = "square", type = "upper",
                   tl.col = "black")


#model
model3<-glm(parity_rate~insecticide_control+season+location_combine+species_complex*host_sampling+zone,family="binomial",data=df,weights=parity_total)
summary(model3)

model4<-glm(parity_rate~insecticide_control+season+location_combine+complex+host_sampling+climate,family="binomial",data=df,weights=parity_total)
summary(model4)

model5<-glm(parity_rate~insecticide_control+season+location_combine*species_complex+host_sampling+zone,family="binomial",data=df,weights=parity_total)
summary(model5)

model6<-glm(parity_rate~insecticide_control+season+location_combine+host_sampling+complex*climate,family="binomial",data=df,weights=parity_total)
summary(model6)

lrtest(model3,model4)

library(epiDisplay)
lrtest(model4,model41)

tibble(
  "risk factor"=names(coef(model4)),
  "Coefficient" = round(coef(model4), 2),
  "SE" = round(se.coef(model4), 2),
  "OR" = round(exp(coef(model4)), 2),
  "95%CI Low limit OR" = round(exp(confint.default(model4, level = 0.95)[,1]), 2),
  "95%CI Upper limit OR" = round(exp(confint.default(model4, level =
                                                       0.95)[,2]), 2),
  "wald_z" = coef(model4)/se.coef(model4)
)%>%
  rowwise() %>%
  mutate("wald_p_value" = round(2*min(pnorm(1-wald_z), pnorm(wald_z)),2)) %>%    
  flextable()


#missing data
library(naniar)
library(GGally)
library(VIM)
library(FactoMineR)
library(missMDA)
library(mice)
library(finalfit)
library(naniar)
library(FactoMineR)
library(reshape2)

explanatory = c("complex","climate","insecticide_control","season","location_combine","host_sampling")
explanatory = c("insecticide_control","season","location_combine","host_sampling")
dependent = "parity_rate"

df %>% filter(df$parity_total>0)%>%
  missing_pairs(dependent, explanatory, position = "fill")
df %>%
  ff_glimpse(dependent, explanatory)

df %>% 
  summary_factorlist(dependent, explanatory, 
                     na_include=TRUE, p=TRUE)

df %>% 
  finalfit.glm(dependent,explanatory,weights = "parity_total") ->t
t

df%>%
  ff_interaction(complex,climate) %>% 
  finalfit.glm(dependent, c(explanatory, "complex_climate"),weights = "parity_total") -> h
h


#mice
keep <- c("country","host_sampling","latitude","longitude","insecticide_control","season","complex","location_combine","zone","climate","parity_total","parity_rate")
dfmiss = df[,(names(df) %in% keep)]%>% filter(df$parity_total>0)

#missMDA
#keep <- c("country","host_sampling","latitude","longitude","insecticide_control","season","complex","location_combine","zone","climate","parity_total","parity_rate")
keep <- c("host_sampling","insecticide_control","season","complex","location_combine","climate","parity_total","parity_rate")

dfmiss = df[,(names(df) %in% keep)]%>% filter(df$parity_total>0)



dim(na.omit(dfmiss))


gg_miss_var(dfmiss)
hi<-summary(aggr(dfmiss, sortVar=TRUE))$combinations
vis_miss(dfmiss, sort_miss = TRUE)
mcar_test(dfmiss)


summary(dfmiss)
marginplot(dfmiss[c(4,9)])

#getting to know data
head(dfmiss)
str(dfmiss)
md.pattern(dfmiss)
res<-summary(aggr(dfmiss, sortVar=TRUE))

#number and proportion of complete cases
Ncc <- cbind(
  "#" = table(complete.cases(dfmiss)),
  "%" = round(100 * table(complete.cases(dfmiss))/nrow(dfmiss), 2)
)
rownames(Ncc) <- c("incompl.", "complete")
Ncc

# number and proportion of missing values per variable
cbind("# NA" = sort(colSums(is.na(dfmiss))),
      "% NA" = round(sort(colMeans(is.na(dfmiss))) * 100, 2))

# categorical and continuous
nc <- max(5, ceiling(sqrt(ncol(dfmiss))))
nr <- ceiling(ncol(dfmiss) / nc)
par(mfrow = c(nr, nc), mgp = c(2, 0.6, 0), mar = c(2, 3, 3, 0.5))
for (i in 1:ncol(dfmiss)) {
  if (is.numeric(dfmiss[, i])) {
    hist(dfmiss[, i], nclass = 50, xlab = "",
         main = paste0(names(dfmiss[i]), " (",
                       round(mean(is.na(dfmiss[, i])) * 100, 2), "% NA)")
    )
  } else {
    barplot(table(dfmiss[, i]), ylab = "Frequency",
            main = paste0(names(dfmiss[i]), " (",
                          round(mean(is.na(dfmiss[, i])) * 100, 2), "% NA)"))
  }
}

#imputation Set-up run
imp0 <- mice(dfmiss, maxit = 0)

#predictor matrix
pred <- imp0$predictorMatrix
pred
pred[c("country"),] <- 1
pred[,c("country")] <- 1
pred[c("country"),c("country")] <- 0


meth <- imp0$meth
meth

#Visit sequence
visSeq <- imp0$visitSequence


#run the imputation
imp1 <- mice(dfmiss, method = meth, predictorMatrix = pred,
            maxit = 10, m = 50,
            seed = 2023,printFlag = FALSE)

#post peocessing
post <- imp1$post
#post["year_start"] <- "imp[[j]][,i] <- squeeze(imp[[j]][,i], c(1984, 1990))"

imp2 <- update(imp1, post = post, maxit = 10, seed = 2023)


#check
class(imp2)
imp2$method
imp2$predictorMatrix
imp2$visitSequence
imp2$loggedEvents

#plot for convergence

plot(imp1, layout = c(3, 6))

ggplot(melt(imp2$chainMean[c("insecticide_control","location_combine","host_sampling","year_start"),,]), aes(x = Var2, y = value, color = Var3)) +
  geom_line() +
  facet_wrap("Var1", scales = 'free') +
  theme(legend.position = 'none') +
  xlab("iteration") +
  ylab("imputed value")

#impute values

longDF <- complete(imp1, "long")
summary(longDF)

lapply(imp2$imp2, function(x) 
  summary(unlist(x))
)


propplot <- function(x, formula, facet = "wrap", ...) {
  library(ggplot2)
  
  cd <- data.frame(mice::complete(x, "long", include = TRUE))
  cd$.imp <- factor(cd$.imp)
  
  r <- as.data.frame(is.na(x$data))
  
  impcat <- x$meth != "" & sapply(x$data, is.factor)
  vnames <- names(impcat)[impcat]
  
  if (missing(formula)) {
    formula <- as.formula(paste(paste(vnames, collapse = "+",
                                      sep = ""), "~1", sep = ""))
  }
  
  tmsx <- terms(formula[-3], data = x$data)
  xnames <- attr(tmsx, "term.labels")
  xnames <- xnames[xnames %in% vnames]
  
  if (paste(formula[3]) != "1") {
    wvars <- gsub("[[:space:]]*\\|[[:print:]]*", "", paste(formula)[3])
    # wvars <- all.vars(as.formula(paste("~", wvars)))
    wvars <- attr(terms(as.formula(paste("~", wvars))), "term.labels")
    if (grepl("\\|", formula[3])) {
      svars <- gsub("[[:print:]]*\\|[[:space:]]*", "", paste(formula)[3])
      svars <- all.vars(as.formula(paste("~", svars)))
    } else {
      svars <- ".imp"
    }
  } else {
    wvars <- NULL
    svars <- ".imp"
  }
  
  for (i in seq_along(xnames)) {
    xvar <- xnames[i]
    select <- cd$.imp != 0 & !r[, xvar]
    cd[select, xvar] <- NA
  }
  
  
  for (i in which(!wvars %in% names(cd))) {
    cd[, wvars[i]] <- with(cd, eval(parse(text = wvars[i])))
  }
  
  meltDF <- reshape2::melt(cd[, c(wvars, svars, xnames)], id.vars = c(wvars, svars))
  meltDF <- meltDF[!is.na(meltDF$value), ]
  
  
  wvars <- if (!is.null(wvars)) paste0("`", wvars, "`")
  
  a <- plyr::ddply(meltDF, c(wvars, svars, "variable", "value"), plyr::summarize,
                   count = length(value))
  b <- plyr::ddply(meltDF, c(wvars, svars, "variable"), plyr::summarize,
                   tot = length(value))
  mdf <- merge(a,b)
  mdf$prop <- mdf$count / mdf$tot
  
  plotDF <- merge(unique(meltDF), mdf)
  plotDF$value <- factor(plotDF$value,
                         levels = unique(unlist(lapply(x$data[, xnames], levels))),
                         ordered = T)
  
  p <- ggplot(plotDF, aes(x = value, fill = get(svars), y = prop)) +
    geom_bar(position = "dodge", stat = "identity") +
    theme(legend.position = "bottom", ...) +
    ylab("proportion") +
    scale_fill_manual(name = "",
                      values = c("black",
                                 colorRampPalette(
                                   RColorBrewer::brewer.pal(9, "Blues"))(x$m + 3)[1:x$m + 3])) +
    guides(fill = guide_legend(nrow = 1))
  
  if (facet == "wrap")
    if (length(xnames) > 1) {
      print(p + facet_wrap(c("variable", wvars), scales = "free"))
    } else {
      if (is.null(wvars)) {
        print(p)
      } else {
        print(p + facet_wrap(wvars, scales = "free"))
      }
    }
  
  if (facet == "grid")
    if (!is.null(wvars)) {
      print(p + facet_grid(paste(paste(wvars, collapse = "+"), "~ variable"),
                           scales = "free"))
    }
}
propplot(imp1)
propplot(imp1, host_sampling~location_combine)
propplot(imp1, location_combine~host_sampling)
propplot(imp1, location_combine~country+host_sampling)
propplot(imp1, insecticide_control~country)
#fitting model on mids objects
model1 <- with(imp1, glm(parity_rate~insecticide_control+season+location_combine+species_complex+host_sampling+zone,family="binomial",weights=parity_total))
class(model1)
names(model1)
lapply(model1, class)
model1$analyses
summary(model1$analyses[[5]])

#pooling
pooled1 <- pool(model1)
summary(pooled1)

###missMDA
keep <- c("country","host_sampling","latitude","longitude","year_start","insecticide_control","season","complex","location_combine","climate","zone","parity_total","parity_rate")
keep <- c("host_sampling","insecticide_control","season","complex","location_combine","climate","parity_total","parity_rate")
dfmiss = df[,(names(df) %in% keep)]%>% filter(df$parity_total>0)



nbdim <- estim_ncpFAMD(dfmiss) # estimate the number of dimensions to impute 
res.famd <- imputeFAMD(dfmiss, ncp = nbdim$ncp)
mda<-as.data.frame(res.famd$completeObs)


mda<-mda %>% mutate(insecticide_control = relevel(factor(insecticide_control), "f"),
                  #country = relevel(factor(country), "Thailand"),
                  season = relevel(factor(season), "rainy"),
                  complex = relevel(factor(complex), "Maculatus"),
                  location_combine = relevel(factor(location_combine), "location_combine_indoor"),
                  host_sampling = relevel(factor(host_sampling), "host_sampling_biting_w"),
        #         zone=relevel(factor(zone), "tro_savannah"),
                  climate=relevel(factor(climate),"tro"))

model4<-glm(parity_rate~insecticide_control+season+location_combine+complex+host_sampling+climate,family="binomial",data=mda,weights=parity_total)
summary(model4)

#multiple imputation
nbdim <- estim_ncpFAMD(dfmiss)
res.mi<-MIFAMD(dfmiss,ncp = 3 ,nboot=50)
nbdim$ncp


impfamd1<-prelim(res.mi,dfmiss)
fitfamd1 <- with(data=impfamd1,exp=glm(parity_rate~insecticide_control+season+location_combine+complex+host_sampling+climate,family="binomial",weights=parity_total))
res.poolfamd1<-pool(fitfamd1)
summary(res.poolfamd1)
## End(Not run)


