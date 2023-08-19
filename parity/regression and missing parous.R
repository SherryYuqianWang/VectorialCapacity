library(readxl)
library(dplyr)
library(tidyverse, quietly = TRUE)
library(flextable) #create table
library(arm)  # calculate st. error

library(mice)#Multiple imputation with chain equation
library(jomo)#Joint modelling multiple imputation


library(missMDA)#MDA imputation

library(finalfit)#generate table
library(VIM)#plot missing data

library(ggplot2)

library(mdscore)  # for calculating LRT
library(ROCit)


library(lmtest)

library(pscl)

library(kableExtra)
library(lme4)



library(naniar)
library(GGally)

library(FactoMineR)

library(mice)

library(naniar)
library(FactoMineR)
library(reshape2)

df<-as.data.frame(read_xlsx("malaria/Statistics/parity/parous_climate.xlsx"))

df1<-df %>% filter(df$parity_total>9) %>% mutate(insecticide = relevel(factor(insecticide), "f"),
                  #land = relevel(factor(land), "urban"),
                  land_use = relevel(factor(land_use), "herb"),
                  #IGBP = relevel(factor(IGBP), "wet land"),
                    country = relevel(factor(country), "Thailand"),
                    season = relevel(factor(season), "dry"),
                    species_complex = relevel(factor(species_complex), "Maculatus group"),
                    location = relevel(factor(location), "indoor"),
                    method = relevel(factor(method), "biting_w"),
                    #zone=relevel(factor(zone), "tro_savannah"),
                    climate=relevel(factor(climate),"temperate"),
                  site=as.factor(site))

keep <- c("insecticide","land_use","method","season","species_complex","location","climate","parity_total","parity_rate")
dfmiss = df1[,(names(df1) %in% keep)]%>% filter(df1$parity_total>9) #only include data point with no less than 10 samples


#know the number and proportion of (in)complete case
cctab <- table(complete.cases(dfmiss))
cbind(
  "#" = setNames(cctab, c('incomplete', 'complete')),
  "%" = round(100 * cctab/nrow(dfmiss), 2)
)
cbind("# NA" = sort(colSums(is.na(dfmiss))),
      "% NA" = round(sort(colMeans(is.na(dfmiss))) * 100, 2))



#MICE (passive imputation considering interaction)
keep <- c("country","insecticide","land_use","method","season","species_complex","location","climate","parity_total","parity_rate")
dfmiss = df1[,(names(df1) %in% keep)]%>% filter(df1$parity_total>9)

#make dummy variables
dfmiss$dirtro <- ifelse(dfmiss$climate=="tropical" & dfmiss$species_complex=="An. dirus",1,0)
dfmiss$mintro <- ifelse(dfmiss$climate=="tropical" & dfmiss$species_complex=="An. minimus",1,0)
dfmiss$sintro <- ifelse(dfmiss$climate=="tropical" & dfmiss$species_complex=="An. sinensis",1,0)

#imputation Set-up run
ini <- mice(dfmiss, max = 0, print = FALSE)

meth <- ini$method
meth["dirus"] <- "~I(species_complex*climate.cold)"
meth["tropical"] <- "~I(species_complex*climate.tropical)"


#predictor matrix
pred <-ini$predictorMatrix
pred
#pred["insecticide","country"] <- 0

meth <- ini$method
meth

#Visit sequence
visSeq <- ini$visitSequence


#run the imputation
imp1 <- mice(dfmiss, method = meth, predictorMatrix = pred,
             maxit = 10, m = 50,
             seed = 2023, printFlag = FALSE)

#post peocessing
#post <- imp1$post
#post["year_start"] <- "imp[[j]][,i] <- squeeze(imp[[j]][,i], c(1984, 1990))"

#imp2 <- update(imp1, post = post, maxit = 10, seed = 2023)


#plot for convergence

plot(imp1, layout = c(3, 6))

#ggplot(melt(imp2$chainMean[c("insecticide","location","method","year_start"),,]), aes(x = Var2, y = value, color = Var3)) +
#  geom_line() +
#  facet_wrap("Var1", scales = 'free') +
#  theme(legend.position = 'none') +
#  xlab("iteration") +
#  ylab("imputed value")

#impute values

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
propplot(imp1, method~location)
propplot(imp1, location~method)
propplot(imp1, location~country+method)
propplot(imp1, insecticide~country)
#fitting model on mids objects
micemodel1 <- with(imp1, glm(parity_rate~insecticide+season+location+method+land_use+species_complex*climate,family="binomial",weights=parity_total))

#pooling
poolmice1 <- pool(micemodel1)
summary(poolmice1)

mice1 = poolmice1 %>%                                  
  fit2df(estimate_name = "OR (MICE)", exp = TRUE)

#another way: use dummy variable
micemodel2 <- with(imp1, glm(parity_rate~insecticide+season+location+method+species_complex+climate+dirtro+mintro+sintro,family="binomial",weights=parity_total))
poolmice2 <- pool(micemodel2)

mice2 = poolmice2 %>%                                  
  fit2df(estimate_name = "OR (MICE)", exp = TRUE)



#Joint modeling (jomo)
keep <- c("insecticide","method","season","species_complex","location","climate","parity_total","parity_rate")
dfmiss = df1[,(names(df1) %in% keep)]%>% filter(df1$parity_total>9)

mylevel<-c(1,1,1,1,1,1,1,1,1,1,1,1)
formula<-as.formula(parity_rate~insecticide+season+location+method+species_complex*climate)
imp<-jomo.glm(formula, dfmiss,  nburn=10, nbetween=10, family = "binomial")


level=mylevel,



###missMDA
#keep <- c("country","method","latitude","longitude","year_start","insecticide","season","complex","location","climate","zone","parity_total","parity_rate")
#eep <- c("method","insecticide","season","complex","location","climate","parity_total","parity_rate")

keep <- c("country","land_use","insecticide","method","season","species_complex","location","climate","parity_total","parity_rate")
dfmiss = df1[,(names(df1) %in% keep)]%>% filter(df1$parity_total>9)

#make dummy variables
dfmiss$dirtro <- ifelse(dfmiss$climate=="tropical" & dfmiss$species_complex=="An. dirus",1,0)
dfmiss$mintro <- ifelse(dfmiss$climate=="tropical" & dfmiss$species_complex=="An. minimus",1,0)
dfmiss$sintro <- ifelse(dfmiss$climate=="tropical" & dfmiss$species_complex=="An. sinensis",1,0)


#single imputation
nbdim <- estim_ncpFAMD(dfmiss) # estimate the number of dimensions to impute 
res.famd <- imputeFAMD(dfmiss, ncp = nbdim$ncp)
mda<-as.data.frame(res.famd$completeObs)

mda<-mda %>% mutate(insecticide = relevel(factor(insecticide), "f"),
                    #country = relevel(factor(country), "Thailand"),
                    season = relevel(factor(season), "rainy"),
                    species_complex = relevel(factor(species_complex), "Maculatus group"),
                    location = relevel(factor(location), "location_indoor"),
                    method = relevel(factor(method), "method_biting_w"),
                    #         zone=relevel(factor(zone), "tro_savannah"),
                    climate=relevel(factor(climate),"tropical"))

model4<-glm(parity_rate~insecticide+season+location+species_complex+method+climate,family="binomial",data=mda,weights=parity_total)
summary(model4)

#multiple imputation
nbdim <- estim_ncpFAMD(dfmiss)
res.mi<-MIFAMD(dfmiss,ncp = nbdim$ncp, nboot=50)
res.mi<-MIFAMD(dfmiss,ncp = 1, nboot=50)

dfmiss <- droplevels(dfmiss)

mda1<-prelim(res.mi,dfmiss)
propplot(mda1, method~location)
propplot(mda1, insecticide~climate)
mdamodel1 <- with(data=mda1,exp=glm(parity_rate~insecticide+season+location+method+land_use+species_complex*climate,family="binomial",weights=parity_total))
poolmda1<-pool(mdamodel1)
summary(poolmda1)

missMDA1 = poolmda1 %>%                                  
  fit2df(estimate_name = "OR (missMDA)", exp = TRUE)

explanatory = c("insecticide","land_use","season","species_complex*climate","location","method")
#explanatory = c("insecticide","land_use","season","location","method","climate")

dependent = "parity_rate"


## summary
dfmiss%>%
  summary_factorlist(dependent, explanatory, fit_id=TRUE,weights="total",digits = c(2, 2, 3, 1, 1)) ->summary1 

## univariable
dfmiss %>% 
    glmuni(dependent, explanatory, weights = "parity_total") %>%
      fit2df(estimate_suffix=" (univariable)")-> uni1

#write.csv(uni1,"parous_uni.csv")
 
## Add multivariable
dfmiss %>% 
    glmmulti(dependent,explanatory,family="binomial",weights = "parity_total") %>%
      fit2df(estimate_suffix=" (multivariable)")->multi1

##
fit <- glm(
  ff_formula(dependent, explanatory), 
  data=dfmiss, family = binomial, weights="parity_total"
)

summary1 %>% 
  ff_merge(uni1) %>% 
  #ff_merge(multi1) %>%
  #ff_merge(mice1)%>%
  #ff_merge(missMDA1)%>%
  knitr::kable(row.names=FALSE)%>%
  remove_column(c(4,5,6))%>%t

write.csv(t,"parous_regression.csv")
class(t)

%>% 
  dependent_label(dfmiss, dependent) -> t
t





#regerssion
model4<-glm(parity_rate~insecticide+season+location+land_use+species_complex*climate+method,family="binomial",data=dfmiss,weights=parity_total)
summary(model4)

mix<-glmer(parity_rate~insecticide+season+land_use+location+species_complex*climate+method+(1|site),family="binomial",weights = parity_total,data=dfmiss,control=glmerControl(optimizer="bobyqa"))
summary(mix)


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

tibble(
  "risk factor"=names(fixef(mix)),
  "Coefficient" = round(fixef(mix), 2),
  #"SE" = round(se.coef(model4), 2),
  "OR" = round(exp(fixef(mix)), 2),
  #"95%CI Low limit OR" = round(exp(confint.default(model4, level = 0.95)[,1]), 2),
  #"95%CI Upper limit OR" = round(exp(confint.default(model4, level =
                           #                            0.95)[,2]), 2),
  #"wald_z" = coef(model4)/se.coef(model4)
)%>%
  rowwise() %>%
 # mutate("wald_p_value" = round(2*min(pnorm(1-wald_z), pnorm(wald_z)),2)) %>%    
  flextable()




## summary
dfmiss %>%
  summary_factorlist(dependent, explanatory, fit_id=TRUE,weights="total",digits = c(2, 2, 3, 1, 1)) ->summary1 

## univariable
dfmiss %>% 
  glmuni(dependent, explanatory, weights = "parity_total") %>%
  fit2df(estimate_suffix=" (univariable)")-> uni1

## Add multivariable
dfmiss %>% 
  glmmulti(dependent,explanatory,family="binomial",weights = "parity_total") %>%
  fit2df(estimate_suffix=" (multivariable)")->multi1


summary1 %>% 
  ff_merge(uni1) %>% 
  ff_merge(multi1) %>%
  ff_merge(mice1)%>%
  ff_merge(missMDA1)%>%
  knitr::kable(row.names=FALSE, align = c("l", "l", "r", "r", "r","r", "r", "r"))