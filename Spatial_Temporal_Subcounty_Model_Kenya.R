#############################################################
inla.getOption("num.threads")

# Load libraries
library(raster)
library(rgdal)
library(sp)
library(spdep)
library(INLA)
library(mgcv)
library(ggplot2)
library(ggregplot)
library(MASS)
library(dplyr)
library(plotly)
library(hrbrthemes)
library(dygraphs)
library(xts)          # To make the conversion data-frame / xts format
library(tidyverse)
library(lubridate)
source("C:/Kenya_Spatial_Temporal_Model/R code/HighstatLibV11.R")
library(ggregplot)    # Plot fixed effects
## edit function from ggregplot
Efxplot.Val <- function (ModelList, Sig = TRUE, StarLoc = NULL, Alpha1 = 1, 
                         Alpha2 = 1, PointOutline = F, ModelNames = NULL, VarNames = NULL, 
                         VarOrder = NULL, Intercept = TRUE, Size = 1, tips = 0.2) 
{
  require(dplyr)
  require(ggplot2)
  require(INLA)
  require(MCMCglmm)
  Graphlist <- list()
  if (!class(ModelList) == "list") {
    ModelList <- list(ModelList)
  }
  for (i in 1:length(ModelList)) {
    model <- ModelList[[i]]
    if (class(model) == "inla") {
      Graph <- as.data.frame(summary(model)$fixed)
      colnames(Graph)[which(colnames(Graph) %in% c("0.025quant", 
                                                   "0.975quant"))] <- c("Lower", "Upper")
      colnames(Graph)[which(colnames(Graph) %in% c("0.05quant", 
                                                   "0.95quant"))] <- c("Lower", "Upper")
      colnames(Graph)[which(colnames(Graph) %in% c("mean"))] <- c("Estimate")
    }
    if (class(model) == "MCMCglmm") {
      Graph <- as.data.frame(summary(model)$solutions)
      colnames(Graph)[1:3] <- c("Estimate", "Lower", 
                                "Upper")
    }
    Graph$Model <- i
    Graph$Factor <- rownames(Graph)
    Graphlist[[i]] <- Graph
  }
  Graph <- bind_rows(Graphlist)
  Graph$Sig <- with(Graph, ifelse(Lower * Upper > 0, "*", 
                                  ""))
  Graph$Model <- as.factor(Graph$Model)
  if (!is.null(ModelNames)) {
    levels(Graph$Model) <- ModelNames
  }
  position <- ifelse(length(unique(Graph$Model)) == 1, "none", 
                     "right")
  if (is.null(VarOrder)) 
    VarOrder <- rev(unique(Graph$Factor))
  if (is.null(VarNames)) 
    VarNames <- VarOrder
  Graph$Factor <- factor(Graph$Factor, levels = VarOrder)
  levels(Graph$Factor) <- VarNames
  Graph %<>% as.data.frame %>% filter(!is.na(Factor))
  if (!Intercept) {
    VarNames <- VarNames[!str_detect(VarNames, "ntercept")]
    Graph <- Graph %>% filter(Factor %in% VarNames)
  }
  Graph$starloc <- NA
  min <- min(Graph$Lower, na.rm = T)
  max <- max(Graph$Upper, na.rm = T)
  if (Sig == TRUE) {
    Graph$starloc <- max + (max - min)/10
  }
  if (!is.null(StarLoc)) {
    Graph$starloc <- StarLoc
  }
  Graph$Alpha <- with(Graph, ifelse(Lower * Upper > 0, Alpha1, 
                                    Alpha2))
  Graph <- Graph %>% mutate(SigAlpha = factor(as.numeric(Lower * 
                                                           Upper > 0), levels = c(1, 0)))
  if (PointOutline) {
    PointOutlineAlpha <- Alpha1
  }
  else {
    PointOutlineAlpha <- 0
  }
  ggplot(Graph, aes(x = as.factor(Factor), y = Estimate, group = Model, 
                    colour = Model, alpha = SigAlpha)) + 
    labs(y = "Linear fixed effect (mean & 95% CI)") +
    theme(text = element_text(size = 18)) +
    geom_point(position = position_dodge(w = 0.5), 
               size = Size) + geom_errorbar(position = position_dodge(w = 0.5), 
                                            aes(ymin = Lower, ymax = Upper), size = 0.3, width = tips) + 
    geom_hline(aes(yintercept = 0), lty = 2) + labs(x = NULL) + 
    coord_flip() + theme(legend.position = position) + geom_text(aes(label = Sig, 
                                                                     y = starloc), position = position_dodge(w = 0.5), show.legend = F) + 
    scale_alpha_manual(values = c(Alpha1, Alpha2)) + guides(alpha = "none") + 
    geom_point(colour = "black", aes(group = Model), 
               position = position_dodge(w = 0.5), size = 4, alpha = PointOutlineAlpha) + 
    geom_errorbar(aes(ymin = Lower, ymax = Upper, group = Model), 
                  width = 0.1, position = position_dodge(w = 0.5), 
                  colour = "black", alpha = PointOutlineAlpha) + 
    geom_point(position = position_dodge(w = 0.5), size = 3, 
               alpha = PointOutlineAlpha)
}
#############################################################




#############################################################
# Load data
setwd("C:/Kenya_Spatial_Temporal_Model")
data <- read.csv(file = "Spatial_temporal_model_data.csv", header = TRUE)

# cases are the number of outbreaks per sub-county per month from 2006 to 2020
names(data)
head(data)
str(data)
data$Date <- as.factor(data$Date)
data$Date <- as.Date(data$Date, format = "%d/%m/%Y")
str(data)
# Add extra variables:
data$County  <- data$ID
#############################################################




#############################################################
# DATA EXPLORATION:


### Time series:
# Aggregate the sum of cases for each month  across all sub-counties 
# so that we have total cases per month across the country from the year 
# 2006 to 2020 (12 months * 15 years = 180 time points)
colnames(data)
data.ts <- data[,c("cases", "Date")] # Extract monthly cases and dates
head(data.ts)
data.ts <- aggregate(. ~ Date, data.ts, sum) # Aggregate sum of cases per month
head(data.ts)
# Usual area chart
p <- data.ts %>%
  ggplot( aes(x=Date, y=cases)) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  geom_area(fill="#00CCCC", alpha=0.5) + # col=#00CCCC
  #geom_line(color="#009999") +
  ylab("Reported Anthrax Cases Per Month") +
  xlab("Reporting Time") +
  theme_classic()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) #,face="bold"
p




### Aggregate observations: 
# From cases per months to number of cases per sub-county per quarter:
# Currently, cases are the number of outbreaks per sub-county per month 
# from Jan-2006 to Dec-2020. So, we aggregate sub-county cases by quarter.
# Resulting in 290 sub-counties * 4 quarters * 15 years = 17400 observations
# 
# Add quarter variable
quarter <- vector() 
for(i in 1:60){
  a <- rep(i,870)
  quarter <- c(quarter, a)
}
# Add as a new column
data$quarter <- quarter

## Aggregate
colnames(data)
data <- data[,c("ADM2_EN","ADM1_EN","ID","Year","cases", "Total.Population",
                "Male.Population","Female.Population","Population.Density",
                "Dairy.Exotic.Cattle","Beef.Exotic.Cattle","Indigenous.Cattle",
                "Agricultural.Land.Area","Farming.Households","Crop.Producing.Households",
                "Livestock.Prod.Households","County", "quarter")]
head(data)
data <- aggregate(. ~ ADM2_EN+ADM1_EN+ID+Year+Total.Population+
                      Male.Population+Female.Population+Population.Density+
                      Dairy.Exotic.Cattle+Beef.Exotic.Cattle+Indigenous.Cattle+
                      Agricultural.Land.Area+Farming.Households+Crop.Producing.Households+
                      Livestock.Prod.Households+County+quarter, data, sum) # Aggregate sum of cases per month
head(data)
dim(data)
## Adding response variables for Bernoulli part:
data$cases01 <- ifelse(data$cases==0, 0, 1) # new response value Bernoulli part





### Outliers:
# These cleveland dotplots show the values of the variables across the different case 
# locations. variables containing outliers should be log transformed prior to use. 
par(mfrow = c(1,1), mar = c(5,5,2,2), cex.lab = 1.5)
plot(x = data$cases,
     y = 1:nrow(data),
     xlab = "Anthrax outbreaks",
     ylab = "Order of the data")
colnames(data)
MyVar <- c("Total.Population","Male.Population","Female.Population",
           "Population.Density","Dairy.Exotic.Cattle","Beef.Exotic.Cattle","Indigenous.Cattle",
           "Agricultural.Land.Area","Farming.Households","Crop.Producing.Households",
           "Livestock.Prod.Households")   
library(lattice)
#Mydotplot(data[,MyVar])





### Collinearity:
# Used VIF values, Pearson correlations, and scatterplots
MyVar <- c("Total.Population","Male.Population","Female.Population",
           "Population.Density","Dairy.Exotic.Cattle","Beef.Exotic.Cattle","Indigenous.Cattle",
           "Agricultural.Land.Area","Farming.Households","Crop.Producing.Households",
           "Livestock.Prod.Households")  
ke.vif1 <- corvif(data[,MyVar])
# We have some serious collinearity!
# remove  "Male.Population","Female.Population","Beef.Exotic.Cattle",
# "Indigenous.Cattle""Farming.Households","Crop.Producing.Households", 
MyVar <- c("Total.Population","Population.Density","Dairy.Exotic.Cattle",
           "Agricultural.Land.Area","Livestock.Prod.Households")  
ke.vif2 <- corvif(data[,MyVar])
# Save
setwd("C:/Kenya_Spatial_Temporal_Model/Species data")
write.csv(ke.vif1, file = "VIF Results 1.csv")
write.csv(ke.vif2, file = "VIF Results 2.csv")


# Use Pearson correlations and scatterplots
# Corrplot
library(corrplot)
# All variables
MyVar1 <- data.frame(
  Total.Population = data$Total.Population, 
  Male.Population = data$Male.Population, 
  Female.Population = data$Female.Population,
  Population.Density = data$Population.Density, 
  Dairy.Exotic.Cattle = data$Dairy.Exotic.Cattle, 
  Beef.Exotic.Cattle = data$Beef.Exotic.Cattle, 
  Indigenous.Cattle = data$Indigenous.Cattle, 
  Agricultural.Land.Area = data$Agricultural.Land.Area, 
  Farming.Households = data$Farming.Households, 
  Crop.Producing.Households = data$Crop.Producing.Households, 
  Livestock.Prod.Households = data$Livestock.Prod.Households)
M = cor(MyVar1)
#corrplot(M, method = 'number', sig.level = 0.01) # colorful number
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(M, method="color", col=col(200),  
         type="upper", order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         # Combine with significance
         sig.level = 0.01, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE 
)

# Scatter plot + Pearson coefficients
library(GGally)
#ggpairs(data[,MyVar])





### Relationships:
MyVar <- c("Total.Population","Population.Density","Dairy.Exotic.Cattle",
           "Agricultural.Land.Area","Livestock.Prod.Households") 
# Incidence/counts
MyMultipanel.ggp2(Z = data, 
                  varx = MyVar, 
                  vary = "cases", 
                  ylab = "Number of outbreaks",
                  addSmoother = TRUE,
                  addRegressionLine = FALSE,
                  addHorizontalLine = FALSE)
# Occurrence: Presence/absence
MyMultipanel.ggp2(Z = data, 
                  varx = MyVar, 
                  vary = "cases01", 
                  ylab = "Number of outbreaks",
                  addSmoother = TRUE,
                  addRegressionLine = FALSE,
                  addHorizontalLine = FALSE)





### Zero inflation:
# Lots of 0s since many areas did not report outbreaks
par(mfrow=c(1,1))
plot(table(data$cases), type = "h", 
     xlab = "Number of outbreaks", 
     ylab = "Number of observations", 
     main = "Frequency of subcounties with given outbreak counts")
round(100 * table(data$cases)/nrow(data), digits = 2)
# Percentages of zeros
sum(data$cases == 0) / nrow(data) # 97.6% of the data has zero counts




# Summary data exploration:
# 1. Clustering: spatial and temporal.
# 2. Correlation between covariates & No clear relationships between response and covariates
# 3. Lots of zeros: zero-inflation
#############################################################




#############################################################
# START INLA ANALYSIS: 

# Bernoulli and ZIP GLM with Spatial-temporal auto-correlation



# Standardize all parametric covariates
MyStd <- function(x) {(x - mean(x)) / sd(x)}
data$Total.Population.std              <- MyStd(data$Total.Population)
data$Population.Density.std            <- MyStd(data$Population.Density)
data$Dairy.Exotic.Cattle.std           <- MyStd(data$Dairy.Exotic.Cattle)
data$Agricultural.Land.Area.std        <- MyStd(data$Agricultural.Land.Area)
data$Livestock.Prod.Households.std     <- MyStd(data$Livestock.Prod.Households)

# Training data 
train <- data # training data
dim(train) # 

# Add polygon for the BYM correlation
setwd("C:/Kenya_Spatial_Temporal_Model")
ken.shp <- readOGR("Kenya_subcounty.shp")
plot (ken.shp)
names(ken.shp@data)
# Calculating the number and identity of neighbors for each sub-county
ken.nb <- poly2nb(ken.shp)
ken.nb
Coords <- coordinates(ken.shp)
plot(ken.shp, border=grey(0.5))  
plot(ken.nb, coords = Coords,add = TRUE, pch = 16, lwd = 2)
nb2INLA("ken.graph", ken.nb)
ken.adj <- paste(getwd(), "ken.graph", sep ="/")
# Identity of neighbors
ken.Inla.nb <- inla.read.graph(filename = "ken.graph")
inla.debug.graph("ken.graph")

# Add priors for BYM correlation hyperparameters
HyperBYM <- list(
  prec.unstruct=list(prior = "pc.prec",param = c(0.01, 0.5)), #param = c(0.01, 0.5))
  prec.spatial=list(prior = "pc.prec",param = c(0.01, 0.5))) #param = c(0.01, 0.5))

#############################################################




#############################################################
# BERNOULLI MODEL - OCCURRENCE:



# Baselines (Intercept only)
I1.b <- inla(formula = cases01 ~  1, 
            family = "binomial", 
            data = train,
            control.compute = list(dic = TRUE))

# Plus covariates
I2.b <- inla(formula = cases01 ~ Total.Population.std + Population.Density.std +   
               Dairy.Exotic.Cattle.std + Agricultural.Land.Area.std + 
               Livestock.Prod.Households.std, 
             family = "binomial", 
             data = train,
             control.compute = list(dic = TRUE))

# Plus covariates + SRF
I3.b <- inla(formula = cases01 ~ Total.Population.std + Population.Density.std +   
               Dairy.Exotic.Cattle.std + Agricultural.Land.Area.std + 
               Livestock.Prod.Households.std + 
               f(County, 
                 model = "bym", 
                 graph = ken.adj,
                 hyper = HyperBYM,
                 adjust.for.con.comp = FALSE,
                 constr = TRUE,
                 scale.model = TRUE), 
            family = "binomial", 
            data = train,
            #verbose=TRUE,
            control.compute = list(dic = TRUE))

# Plus covariates + SRF + Temporal Trend
I4.b <- inla(formula = cases01 ~ Total.Population.std + Population.Density.std +   
               Dairy.Exotic.Cattle.std + Agricultural.Land.Area.std + 
               Livestock.Prod.Households.std + 
               f(County, 
                 model = "bym", 
                 graph = ken.adj,
                 hyper = HyperBYM,
                 adjust.for.con.comp = FALSE,
                 constr = TRUE,
                 scale.model = TRUE)  +
               f(quarter,model = "rw1",
                 scale.model = TRUE,
                 hyper = list(theta = list(prior="pc.prec", param=c(1,0.05)))), 
             family = "binomial", 
             data = train,
             control.compute = list(dic = TRUE))





### Model Selection and Validation: 


# Compare DIC
dic  <- c(I1.b$dic$dic,I2.b$dic$dic,I3.b$dic$dic,I4.b$dic$dic)   
Result     <- cbind(dic)
rownames(Result) <- c("Bernoulli GLM: Baseline",  
                      "Bernoulli GLM: Baseline + covariates",
                      "Bernoulli GLM: Baseline + covariates + SRF",
                      "Bernoulli GLM: Baseline + covariates + SRF + temporal trend")
Result
# The final model with the lowest DIC is **I4.b**
## Save
setwd("C:/Kenya_Spatial_Temporal_Model")
class(Result)
Result <- as.data.frame(Result)
write.csv(Result, file = "Quarterly Model Validation Results - Bernoulli.csv")



# Test robustness of final model by geographical cross-validation:

## Exclude observations from sub-counties reporting high incidence
# Model validation  for robustness
train1 <- train[!(train$ADM2_EN %in% c("Langata","Roysambu", "Westlands" )),] # minus Nairobi sub-counties
train2 <- train[!(train$ADM2_EN %in% c("Githunguri","Kabete","Kiambu","Kikuyu","Limuru","Ruiru","Thika Town" )),] # minus Kiambu sc
train3 <- train[!(train$ADM2_EN %in% c("Kiharu","Maragwa")),] # minus Muranga sc
train4 <- train[!(train$ADM2_EN %in% c("Mathira","Mukurweni")),] # minus Nyeri sc
train5 <- train[!train$ADM2_EN == "Maara",] # minus Tharaka-Nithi sc
train6 <- train[!(train$ADM2_EN %in% c("Nakuru Town East","Nakuru Town West", "Rongai" )),] # minus Nakuru sc
train7 <- train[!train$ADM2_EN == "Ainamoi",] # minus Kericho sc
train8 <- train[!(train$ADM2_EN %in% c("Kapseret","Soy", "Turbo")),] # minus Uasin Gishu sc
train9 <- train[!train$ADM2_EN == "Sotik",] # minus Bomet sc
train10 <-train[!train$ADM2_EN == "Kaloleni",] # minus Kilifi sc
train11 <- train[!train$ADM2_EN == "Kisauni",] # minus Mombasa


# Run (test robustness of best/Mb4 model)
Best.1 <- inla(formula = cases01 ~ Total.Population.std + Population.Density.std +   
                 Dairy.Exotic.Cattle.std + Agricultural.Land.Area.std + 
                 Livestock.Prod.Households.std + 
                 f(County, 
                   model = "bym", 
                   graph = ken.adj,
                   hyper = HyperBYM,
                   adjust.for.con.comp = FALSE,
                   constr = TRUE,
                   scale.model = TRUE)  +
                 f(quarter,model = "rw1",
                   scale.model = TRUE,
                   hyper = list(theta = list(prior="pc.prec", param=c(1,0.05)))), 
               family = "binomial",
               data = train1,
               control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE))

Best.2 <- inla(formula = cases01 ~ Total.Population.std + Population.Density.std +   
                 Dairy.Exotic.Cattle.std + Agricultural.Land.Area.std + 
                 Livestock.Prod.Households.std + 
                 f(County, 
                   model = "bym", 
                   graph = ken.adj,
                   hyper = HyperBYM,
                   adjust.for.con.comp = FALSE,
                   constr = TRUE,
                   scale.model = TRUE)  +
                 f(quarter,model = "rw1",
                   scale.model = TRUE,
                   hyper = list(theta = list(prior="pc.prec", param=c(1,0.05)))), 
               family = "binomial",
               data = train2,
               control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE))

Best.3 <- inla(formula = cases01 ~ Total.Population.std + Population.Density.std +   
                 Dairy.Exotic.Cattle.std + Agricultural.Land.Area.std + 
                 Livestock.Prod.Households.std + 
                 f(County, 
                   model = "bym", 
                   graph = ken.adj,
                   hyper = HyperBYM,
                   adjust.for.con.comp = FALSE,
                   constr = TRUE,
                   scale.model = TRUE)  +
                 f(quarter,model = "rw1",
                   scale.model = TRUE,
                   hyper = list(theta = list(prior="pc.prec", param=c(1,0.05)))), 
               family = "binomial", 
               data = train3,
               control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE))

Best.4 <- inla(formula = cases01 ~ Total.Population.std + Population.Density.std +   
                 Dairy.Exotic.Cattle.std + Agricultural.Land.Area.std + 
                 Livestock.Prod.Households.std + 
                 f(County, 
                   model = "bym", 
                   graph = ken.adj,
                   hyper = HyperBYM,
                   adjust.for.con.comp = FALSE,
                   constr = TRUE,
                   scale.model = TRUE)  +
                 f(quarter,model = "rw1",
                   scale.model = TRUE,
                   hyper = list(theta = list(prior="pc.prec", param=c(1,0.05)))), 
               family = "binomial", 
               data = train4,
               control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE))

Best.5 <- inla(formula = cases01 ~ Total.Population.std + Population.Density.std +   
                 Dairy.Exotic.Cattle.std + Agricultural.Land.Area.std + 
                 Livestock.Prod.Households.std + 
                 f(County, 
                   model = "bym", 
                   graph = ken.adj,
                   hyper = HyperBYM,
                   adjust.for.con.comp = FALSE,
                   constr = TRUE,
                   scale.model = TRUE)  +
                 f(quarter,model = "rw1",
                   scale.model = TRUE,
                   hyper = list(theta = list(prior="pc.prec", param=c(1,0.05)))), 
               family = "binomial", 
               data = train5,
               control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE))

Best.6 <- inla(formula = cases01 ~ Total.Population.std + Population.Density.std +   
                 Dairy.Exotic.Cattle.std + Agricultural.Land.Area.std + 
                 Livestock.Prod.Households.std + 
                 f(County, 
                   model = "bym", 
                   graph = ken.adj,
                   hyper = HyperBYM,
                   adjust.for.con.comp = FALSE,
                   constr = TRUE,
                   scale.model = TRUE)  +
                 f(quarter,model = "rw1",
                   scale.model = TRUE,
                   hyper = list(theta = list(prior="pc.prec", param=c(1,0.05)))), 
               family = "binomial", 
               data = train6,
               control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE))

Best.7 <- inla(formula = cases01 ~ Total.Population.std + Population.Density.std +   
                 Dairy.Exotic.Cattle.std + Agricultural.Land.Area.std + 
                 Livestock.Prod.Households.std + 
                 f(County, 
                   model = "bym", 
                   graph = ken.adj,
                   hyper = HyperBYM,
                   adjust.for.con.comp = FALSE,
                   constr = TRUE,
                   scale.model = TRUE)  +
                 f(quarter,model = "rw1",
                   scale.model = TRUE,
                   hyper = list(theta = list(prior="pc.prec", param=c(1,0.05)))), 
               family = "binomial", 
               data = train7,
               control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE))

Best.8 <- inla(formula = cases01 ~ Total.Population.std + Population.Density.std +   
                 Dairy.Exotic.Cattle.std + Agricultural.Land.Area.std + 
                 Livestock.Prod.Households.std + 
                 f(County, 
                   model = "bym", 
                   graph = ken.adj,
                   hyper = HyperBYM,
                   adjust.for.con.comp = FALSE,
                   constr = TRUE,
                   scale.model = TRUE)  +
                 f(quarter,model = "rw1",
                   scale.model = TRUE,
                   hyper = list(theta = list(prior="pc.prec", param=c(1,0.05)))), 
               family = "binomial", 
               data = train8,
               control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE))

Best.9 <- inla(formula = cases01 ~ Total.Population.std + Population.Density.std +   
                 Dairy.Exotic.Cattle.std + Agricultural.Land.Area.std + 
                 Livestock.Prod.Households.std + 
                 f(County, 
                   model = "bym", 
                   graph = ken.adj,
                   hyper = HyperBYM,
                   adjust.for.con.comp = FALSE,
                   constr = TRUE,
                   scale.model = TRUE)  +
                 f(quarter,model = "rw1",
                   scale.model = TRUE,
                   hyper = list(theta = list(prior="pc.prec", param=c(1,0.05)))), 
               family = "binomial", 
               data = train9,
               control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE))

Best.10 <- inla(formula = cases01 ~ Total.Population.std + Population.Density.std +   
                  Dairy.Exotic.Cattle.std + Agricultural.Land.Area.std + 
                  Livestock.Prod.Households.std + 
                  f(County, 
                    model = "bym", 
                    graph = ken.adj,
                    hyper = HyperBYM,
                    adjust.for.con.comp = FALSE,
                    constr = TRUE,
                    scale.model = TRUE)  +
                  f(quarter,model = "rw1",
                    scale.model = TRUE,
                    hyper = list(theta = list(prior="pc.prec", param=c(1,0.05)))), 
                family = "binomial", 
               data = train10,
               control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE))

Best.11 <- inla(formula = cases01 ~ Total.Population.std + Population.Density.std +   
                  Dairy.Exotic.Cattle.std + Agricultural.Land.Area.std + 
                  Livestock.Prod.Households.std + 
                  f(County, 
                    model = "bym", 
                    graph = ken.adj,
                    hyper = HyperBYM,
                    adjust.for.con.comp = FALSE,
                    constr = TRUE,
                    scale.model = TRUE)  +
                  f(quarter,model = "rw1",
                    scale.model = TRUE,
                    hyper = list(theta = list(prior="pc.prec", param=c(1,0.05)))), 
                family = "binomial", 
               data = train11,
               control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE))

# compare fixed effects/betas:
Efxplot.Val(list(I4.b,Best.1,Best.2,Best.3,Best.4,Best.5,Best.6,Best.7,Best.8,Best.9,Best.10,Best.11), 
        ModelNames = c("Best Occurrence Model","Nairobi", "Kiambu","Muranga","Nyeri","Tharaka-Nithi",
                       "Nakuru","Kericho","Uasin Gishu","Bomet","Kilifi","Mombasa"), Intercept = FALSE)
####################




### Model interpretation:


# plot fixed effects for best model:
Efxplot.Val(list(I4.b))
round(I4.b$summary.fixed[, c("mean", "0.025quant", "0.975quant")], 3)



# Spatial random effects:Bernoulli
# Plot the spatial random effect
ken.shp$SRFOcc <- I4.b$summary.random$County$mean[1:290]
library(RColorBrewer)
crq = brewer.pal(7, "YlGnBu")
range(ken.shp$SRFOcc)
rng = seq(-3, 4, 1)
spplot(ken.shp, "SRFOcc", col = "white", at = rng, 
       col.regions = crq,
       colorkey = list(space = "bottom", labels = paste(rng)),
       par.settings = list(axis.line = list(col = NA)),
       sub = "Spatial random effects - Occurrence")



# Random walk trend/temporal trend and 95% CI: (pg. 270; V1)
# The trend, mu, is a random effect that changes over time.
seasonsm <- I4.b$summary.random$quarter
plot(seasonsm[,1:2], type='l',
     xlab = "Quarter",
     ylab = "Random walk trend",
     ylim = c(-1, 1))
abline(h = 0, lty = 3)
lines (seasonsm[,c(1,4)], lty = 2)
lines (seasonsm[,c(1,6)], lty = 2)


# Posterior distribution of the sigma
# We have one variance in the model, but R-INLA uses precision parameters,
# tau_v = 1/sigma_v^2 (where sigma_e^2 is the variance for the pure noise term (v) in the random walk trend - mu_t = mu_t-1 + v) 

# All the numerical output is for precision, but we can convert it to posterior
# info on sigma_v. The marginal distribution of the precision parameter 
# is obtained by the marginals.hyperpar
hp <- I4.b$marginals.hyperpar
tau.v <- hp$`Precision for quarter`
# We want the standard deviation parameter sigma_v. First we calculate the posterior
# mean value of sigma_v
Tau2Sigma <- function(x){ sqrt(1/x)}
sigma.v <- inla.emarginal(fun=Tau2Sigma, marg = tau.v)
sigma.v # 0.4998675
# Marginal posterior distribution
md.sigma.v <- inla.tmarginal(fun=Tau2Sigma, marg = tau.v)
plot (x = md.sigma.v[,1], y = md.sigma.v[,2], type = "l",
      xlab = expression(sigma),
      ylab = expression(paste("P(", sigma ," | Data)")))





# fitted values - Occurrence/Bernoulli:
train2020 <- train[17111:17400,] # data for last quarter of 2020
train2020$fittedOcc <- I4.b$summary.fitted.values[17111:17400,"mean"]*5 # fitted values 
# Plot the spatial fitted values
ken.shp$fittedOcc <- train2020$fittedOcc
# plot
library(RColorBrewer)
crq = brewer.pal(9, "YlGnBu")
range(ken.shp$fittedOcc) # 0.001467997 0.897064860
rng = seq(0, 1, 0.1)
spplot(ken.shp, "fittedOcc", col = "white", at = rng, 
       col.regions = crq,
       colorkey = list(space = "bottom", labels = paste(rng)),
       par.settings = list(axis.line = list(col = NA)),
       sub = "Fitted Mean - Probability of Occurrence")
# Save
writeOGR(ken.shp, 
         dsn="C:/Kenya_Spatial_Temporal_Model/Fitted values", 
         layer="OccurrenceProbability2020_PlusSRF_bernoulli",
         driver="ESRI Shapefile")
#############################################################




#############################################################
# ZIP MODEL - OCCURRENCE:


# Baselines (Intercept only)
I1.zip <- inla(formula = cases ~  1, 
             family = "zeroinflatedpoisson1", 
             data = train,
             control.compute = list(dic = TRUE))

# Plus covariates
I2.zip <- inla(formula = cases ~ Total.Population.std + Population.Density.std +   
               Dairy.Exotic.Cattle.std + Agricultural.Land.Area.std + 
               Livestock.Prod.Households.std, 
             family = "zeroinflatedpoisson1", 
             data = train,
             control.compute = list(dic = TRUE))

# Plus covariates + SRF
I3.zip <- inla(formula = cases ~ Total.Population.std + Population.Density.std +   
               Dairy.Exotic.Cattle.std + Agricultural.Land.Area.std + 
               Livestock.Prod.Households.std + 
               f(County, 
                 model = "bym", 
                 graph = ken.adj,
                 hyper = HyperBYM,
                 adjust.for.con.comp = FALSE,
                 constr = TRUE,
                 scale.model = TRUE), 
             family = "zeroinflatedpoisson1", 
             data = train,
             #verbose=TRUE,
             control.compute = list(dic = TRUE))

# Plus covariates + SRF + Annual Trend
I4.zip <- inla(formula = cases ~ Total.Population.std + Population.Density.std +   
               Dairy.Exotic.Cattle.std + Agricultural.Land.Area.std + 
               Livestock.Prod.Households.std + 
               f(County, 
                 model = "bym", 
                 graph = ken.adj,
                 hyper = HyperBYM,
                 adjust.for.con.comp = FALSE,
                 constr = TRUE,
                 scale.model = TRUE)  +
               f(quarter,model = "rw1",
                 scale.model = TRUE,
                 hyper = list(theta = list(prior="pc.prec", param=c(1,0.05)))), 
             family = "zeroinflatedpoisson1", 
             data = train,
             control.compute = list(dic = TRUE))





### Model Selection and Validation: 


# Compare DIC
dic  <- c(I1.zip$dic$dic, I2.zip$dic$dic, I3.zip$dic$dic, I4.zip$dic$dic)   
Result     <- cbind(dic)
rownames(Result) <- c("ZIP GLM: Baseline",  
                      "ZIP GLM: Baseline + covariates",
                      "ZIP GLM: Baseline + covariates + SRF",
                      "ZIP GLM: Baseline + covariates + SRF + temporal trend")
Result
## The final model with the lowest DIC is **I4.zip**
# Save
setwd("C:/Kenya_Spatial_Temporal_Model")
class(Result)
Result <- as.data.frame(Result)
write.csv(Result, file = "Quarterly Model Validation Results - ZIP.csv")




### Test robustness by geographical cross-validation

train1 <- train[!(train$ADM2_EN %in% c("Langata","Roysambu", "Westlands" )),] # minus Nairobi sub-counties
train2 <- train[!(train$ADM2_EN %in% c("Githunguri","Kabete","Kiambu","Kikuyu","Limuru","Ruiru","Thika Town" )),] # minus Kiambu sc
train3 <- train[!(train$ADM2_EN %in% c("Kiharu","Maragwa")),] # minus Muranga sc
train4 <- train[!(train$ADM2_EN %in% c("Mathira","Mukurweni")),] # minus Nyeri sc
train5 <- train[!train$ADM2_EN == "Maara",] # minus Tharaka-Nithi sc
train6 <- train[!(train$ADM2_EN %in% c("Nakuru Town East","Nakuru Town West", "Rongai" )),] # minus Nakuru sc
train7 <- train[!train$ADM2_EN == "Ainamoi",] # minus Kericho sc
train8 <- train[!(train$ADM2_EN %in% c("Kapseret","Soy", "Turbo")),] # minus Uasin Gishu sc
train9 <- train[!train$ADM2_EN == "Sotik",] # minus Bomet sc
train10 <-train[!train$ADM2_EN == "Kaloleni",] # minus Kilifi sc
train11 <- train[!train$ADM2_EN == "Kisauni",] # minus Mombasa

# Run 
Best.1.zip <- inla(formula = cases ~ Total.Population.std + Population.Density.std +   
                     Dairy.Exotic.Cattle.std + Agricultural.Land.Area.std + 
                     Livestock.Prod.Households.std + 
                     f(County, 
                       model = "bym", 
                       graph = ken.adj,
                       hyper = HyperBYM,
                       adjust.for.con.comp = FALSE,
                       constr = TRUE,
                       scale.model = TRUE)  +
                     f(quarter,model = "rw1",
                       scale.model = TRUE,
                       hyper = list(theta = list(prior="pc.prec", param=c(1,0.05)))), 
                   family = "zeroinflatedpoisson1", 
                   data = train1,
                   control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE))

Best.2.zip <- inla(formula = cases ~ Total.Population.std + Population.Density.std +   
                     Dairy.Exotic.Cattle.std + Agricultural.Land.Area.std + 
                     Livestock.Prod.Households.std + 
                     f(County, 
                       model = "bym", 
                       graph = ken.adj,
                       hyper = HyperBYM,
                       adjust.for.con.comp = FALSE,
                       constr = TRUE,
                       scale.model = TRUE)  +
                     f(quarter,model = "rw1",
                       scale.model = TRUE,
                       hyper = list(theta = list(prior="pc.prec", param=c(1,0.05)))), 
                   family = "zeroinflatedpoisson1", 
                 data = train2,
                 control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE))

Best.3.zip <- inla(formula = cases ~ Total.Population.std + Population.Density.std +   
                     Dairy.Exotic.Cattle.std + Agricultural.Land.Area.std + 
                     Livestock.Prod.Households.std + 
                     f(County, 
                       model = "bym", 
                       graph = ken.adj,
                       hyper = HyperBYM,
                       adjust.for.con.comp = FALSE,
                       constr = TRUE,
                       scale.model = TRUE)  +
                     f(quarter,model = "rw1",
                       scale.model = TRUE,
                       hyper = list(theta = list(prior="pc.prec", param=c(1,0.05)))), 
                   family = "zeroinflatedpoisson1", 
                 data = train3,
                 control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE))

Best.4.zip <- inla(formula = cases ~ Total.Population.std + Population.Density.std +   
                     Dairy.Exotic.Cattle.std + Agricultural.Land.Area.std + 
                     Livestock.Prod.Households.std + 
                     f(County, 
                       model = "bym", 
                       graph = ken.adj,
                       hyper = HyperBYM,
                       adjust.for.con.comp = FALSE,
                       constr = TRUE,
                       scale.model = TRUE)  +
                     f(quarter,model = "rw1",
                       scale.model = TRUE,
                       hyper = list(theta = list(prior="pc.prec", param=c(1,0.05)))), 
                   family = "zeroinflatedpoisson1", 
                 data = train4,
                 control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE))

Best.5.zip <- inla(formula = cases ~ Total.Population.std + Population.Density.std +   
                     Dairy.Exotic.Cattle.std + Agricultural.Land.Area.std + 
                     Livestock.Prod.Households.std + 
                     f(County, 
                       model = "bym", 
                       graph = ken.adj,
                       hyper = HyperBYM,
                       adjust.for.con.comp = FALSE,
                       constr = TRUE,
                       scale.model = TRUE)  +
                     f(quarter,model = "rw1",
                       scale.model = TRUE,
                       hyper = list(theta = list(prior="pc.prec", param=c(1,0.05)))), 
                   family = "zeroinflatedpoisson1", 
                 data = train5,
                 control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE))

Best.6.zip <- inla(formula = cases ~ Total.Population.std + Population.Density.std +   
                     Dairy.Exotic.Cattle.std + Agricultural.Land.Area.std + 
                     Livestock.Prod.Households.std + 
                     f(County, 
                       model = "bym", 
                       graph = ken.adj,
                       hyper = HyperBYM,
                       adjust.for.con.comp = FALSE,
                       constr = TRUE,
                       scale.model = TRUE)  +
                     f(quarter,model = "rw1",
                       scale.model = TRUE,
                       hyper = list(theta = list(prior="pc.prec", param=c(1,0.05)))), 
                   family = "zeroinflatedpoisson1", 
                 data = train6,
                 control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE))

Best.7.zip <- inla(formula = cases ~ Total.Population.std + Population.Density.std +   
                     Dairy.Exotic.Cattle.std + Agricultural.Land.Area.std + 
                     Livestock.Prod.Households.std + 
                     f(County, 
                       model = "bym", 
                       graph = ken.adj,
                       hyper = HyperBYM,
                       adjust.for.con.comp = FALSE,
                       constr = TRUE,
                       scale.model = TRUE)  +
                     f(quarter,model = "rw1",
                       scale.model = TRUE,
                       hyper = list(theta = list(prior="pc.prec", param=c(1,0.05)))), 
                   family = "zeroinflatedpoisson1", 
                 data = train7,
                 control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE))

Best.8.zip <- inla(formula = cases ~ Total.Population.std + Population.Density.std +   
                     Dairy.Exotic.Cattle.std + Agricultural.Land.Area.std + 
                     Livestock.Prod.Households.std + 
                     f(County, 
                       model = "bym", 
                       graph = ken.adj,
                       hyper = HyperBYM,
                       adjust.for.con.comp = FALSE,
                       constr = TRUE,
                       scale.model = TRUE)  +
                     f(quarter,model = "rw1",
                       scale.model = TRUE,
                       hyper = list(theta = list(prior="pc.prec", param=c(1,0.05)))), 
                   family = "zeroinflatedpoisson1", 
                 data = train8,
                 control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE))

Best.9.zip <- inla(formula = cases ~ Total.Population.std + Population.Density.std +   
                     Dairy.Exotic.Cattle.std + Agricultural.Land.Area.std + 
                     Livestock.Prod.Households.std + 
                     f(County, 
                       model = "bym", 
                       graph = ken.adj,
                       hyper = HyperBYM,
                       adjust.for.con.comp = FALSE,
                       constr = TRUE,
                       scale.model = TRUE)  +
                     f(quarter,model = "rw1",
                       scale.model = TRUE,
                       hyper = list(theta = list(prior="pc.prec", param=c(1,0.05)))), 
                   family = "zeroinflatedpoisson1", 
                 data = train9,
                 control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE))

Best.10.zip <- inla(formula = cases ~ Total.Population.std + Population.Density.std +   
                      Dairy.Exotic.Cattle.std + Agricultural.Land.Area.std + 
                      Livestock.Prod.Households.std + 
                      f(County, 
                        model = "bym", 
                        graph = ken.adj,
                        hyper = HyperBYM,
                        adjust.for.con.comp = FALSE,
                        constr = TRUE,
                        scale.model = TRUE)  +
                      f(quarter,model = "rw1",
                        scale.model = TRUE,
                        hyper = list(theta = list(prior="pc.prec", param=c(1,0.05)))), 
                    family = "zeroinflatedpoisson1",
                  data = train10,
                  control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE))

Best.11.zip <- inla(formula = cases ~ Total.Population.std + Population.Density.std +   
                      Dairy.Exotic.Cattle.std + Agricultural.Land.Area.std + 
                      Livestock.Prod.Households.std + 
                      f(County, 
                        model = "bym", 
                        graph = ken.adj,
                        hyper = HyperBYM,
                        adjust.for.con.comp = FALSE,
                        constr = TRUE,
                        scale.model = TRUE)  +
                      f(quarter,model = "rw1",
                        scale.model = TRUE,
                        hyper = list(theta = list(prior="pc.prec", param=c(1,0.05)))), 
                    family = "zeroinflatedpoisson1", 
                  data = train11,
                  control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE))


# compare fixed effects/betas:
Efxplot.Val(list(I4.zip,Best.1.zip,Best.2.zip, Best.3.zip, Best.4.zip, Best.5.zip, Best.6.zip,
             Best.7.zip, Best.8.zip, Best.9.zip, Best.10.zip, Best.11.zip), 
        ModelNames = c("Best Incidence Model","Nairobi", "Kiambu","Muranga","Nyeri","Tharaka-Nithi",
                       "Nakuru","Kericho","Uasin Gishu","Bomet","Kilifi","Mombasa"),Intercept = FALSE)
#########################




### Model Interpretation: 

# Plot fixed effects:
Efxplot.Val(list(I4.zip))
round(I4.zip$summary.fixed[, c("mean", "0.025quant", "0.975quant")], 3)


# Spatial random effects: ZIP
# Plot the spatial random effect
ken.shp$SRFInc     <- I4.zip$summary.random$County$mean[1:290]
library(RColorBrewer)
crq = brewer.pal(7, "YlGnBu")
range(ken.shp$SRFInc)
rng = seq(-3, 4, 1)
spplot(ken.shp, "SRFInc", col = "white", at = rng, 
       col.regions = crq,
       colorkey = list(space = "bottom", labels = paste(rng)),
       par.settings = list(axis.line = list(col = NA)),
       sub = "Spatial random effects - Incidence")



# Random walk trend/temporal trend and 95% CI: (pg. 270; V1)
# The trend, mu, is a random effect that changes over time.
seasonsm <- I4.zip$summary.random$quarter
plot(seasonsm[,1:2], type='l',
     xlab = "Quarter",
     ylab = "Random walk trend",
     ylim = c(-1, 1))
abline(h = 0, lty = 3)
lines (seasonsm[,c(1,4)], lty = 2)
lines (seasonsm[,c(1,6)], lty = 2)

# Posterior distribution of the sigma:
# We have one variance in the model, but R-INLA uses precision parameters,
# tau_v = 1/sigma_v^2 (where sigma_e^2 is the variance for the pure noise term (v) in the random walk trend - mu_t = mu_t-1 + v) 
# All the numerical output is for precision, but we can convert it to posterior
# info on sigma_v. The marginal distribution of the precision parameter 
# is obtained by the marginals.hyperpar
hp <- I4.zip$marginals.hyperpar
tau.v <- hp$`Precision for quarter`
# We want the standard deviation parameter sigma_v. First we calculate the posterior
# mean value of sigma_v
Tau2Sigma <- function(x){ sqrt(1/x)}
sigma.v <- inla.emarginal(fun=Tau2Sigma, marg = tau.v)
sigma.v # 0.4998675
# Marginal posterior distribution
md.sigma.v <- inla.tmarginal(fun=Tau2Sigma, marg = tau.v)
plot (x = md.sigma.v[,1], y = md.sigma.v[,2], type = "l",
      xlab = expression(sigma),
      ylab = expression(paste("P(", sigma ," | Data)")))



# fitted values - ZIP:
# Formula for ZIP
# count_i = ZIP(mu,Pi)
# Exp(count) = (1-Pi)*mu
# Var(count_i) = (1-Pi)*(mu+Pi*mu^2)
# log(mu) = Intercept+covariates
# logit(Pi) = gamma
train2020 <- train[17111:17400,] # data for last quarter of 2020
# Extract mu:
train2020$fittedInc <- I4.zip$summary.fitted.values[17111:17400,"mean"]# fitted values for  2020
Pi <- I4.zip$summary.hyperpar[1,"mean"]
Pi
train2020$fittedInc <- (1-Pi)*train2020$fittedInc
# Plot the spatial fitted values (incidence per 100,000 people)
ken.shp$fittedInc <- ((train2020$fittedInc * 100000) /train2020$Total.Population)*100
# plot
library(RColorBrewer)
crq = brewer.pal(7, "YlGnBu")
range(ken.shp$fittedInc) # 
rng = seq(0, 14, 2)
spplot(ken.shp, "fittedInc", col = "white", at = rng, 
       col.regions = crq,
       colorkey = list(space = "bottom", labels = paste(rng)),
       par.settings = list(axis.line = list(col = NA)),
       sub = "Fitted Mean - Incidence")
# Save
writeOGR(ken.shp, 
         dsn="C:/Kenya_Spatial_Temporal_Model/Fitted values", 
         layer="Incidence2020_PlusSRF_ZIP",
         driver="ESRI Shapefile")

#############################################################




#############################################################
# Further model validation: ZIP MODEL:


# We first calculate the Pi = probability of false zeros
Pi <- I4.zip$summary.hyper[1,"mean"]
Pi
# Fitted model:
## count_i = ZIP(mu, Pi)
## E(count_i) = (1-Pi) * mu
## Var(count_i) = (1-Pi) * (mu + Pi * mu^2)
## log(mu) = intercept + covariate effects

mu <- I4.zip$summary.fitted.values[,"mean"]
ExpY <- (1-Pi) * mu
VarY <- (1-Pi) * (mu + Pi * mu^2)
E2 <- (train$cases - ExpY) / sqrt(VarY) # Pearson residuals

N <- nrow(train)
p <- length(I4.zip$names.fixed) + 1
Dispersion <- sum(E2^2) / (N-p)
Dispersion #0.45
  


### Fitted vs residuals
par(mfrow = c(1,1), mar = c(5,5,2,2))
plot(x = mu, 
     y = E2,
     xlab = "Fitted values",
     ylab = "Residuals",
     cex.lab = 1.5)
abline(h = 0, lty = 2)  

### Fitted vs observed:
par(mfrow = c(1,1), mar = c(5,5,2,2))
plot(x = mu, 
     y = train$cases,
     xlab = "Fitted values",
     ylab = "Observed values",
     cex.lab = 1.5)
abline(h = 0, lty = 2)     


# Plot residuals vs every covariate       
MyVar <- c("Total.Population","Population.Density","Dairy.Exotic.Cattle",
           "Agricultural.Land.Area","Livestock.Prod.Households")
train$E2 <- E2
MyMultipanel.ggp2(Z = train, 
                  varx = MyVar, 
                  vary = "E2", 
                  ylab = "Pearson residuals",
                  addSmoother = TRUE,
                  addRegressionLine = FALSE,
                  addHorizontalLine = TRUE)


# Not in the book:
# Difficult to assess. Let's use mgcv for this
T0 <- gam(E2 ~ 1, data = train)

Ttotal_pop                 <- gam(E2 ~ s(Total.Population), data = train)
Tpop_density               <- gam(E2 ~ s(Population.Density), data = train)
Texotic_cattleDairy        <- gam(E2 ~ s(Dairy.Exotic.Cattle), data = train)
Tagric_land_area           <- gam(E2 ~ s(Agricultural.Land.Area), data = train)
TlivestockP_households     <- gam(E2 ~ s(Livestock.Prod.Households), data = train)


AIC(T0, Ttotal_pop , Tpop_density, Texotic_cattleDairy, Tagric_land_area, TlivestockP_households)

# plot
par(mfrow = c(2,3))


plot(Ttotal_pop)
abline(h = 0, lty = 2)
#text(0.5, 1.55, "B", cex = 1.5)

plot(Tpop_density)
abline(h = 0, lty = 2)
#text(0.5, 1.55, "B", cex = 1.5)

plot(Texotic_cattleDairy)
abline(h = 0, lty = 2)
#text(0.5, 1.55, "B", cex = 1.5)

plot(Tagric_land_area)
abline(h = 0, lty = 2)
#text(0.5, 1.55, "B", cex = 1.5)

plot(TlivestockP_households)
abline(h = 0, lty = 2)
#text(0.5, 1.55, "B", cex = 1.5)


par(mfrow = c(1,1))

## Do the credible intervals of each covariate contain zero across the entire gradient
## of the covariates? Yes to all except indigenous cattle!! - consider adding 4 knots 




###########################
# Simulation study:
# Simulate data sets from the model to determine whether
# the model can cope with zero inflation (and check overdispersion
# in a more Bayesian way).




# ZIP model with spatial-temporal correlation
## run model
I4.zip.sim <- inla(formula = cases ~ Total.Population.std + Population.Density.std +   
                     Dairy.Exotic.Cattle.std + Agricultural.Land.Area.std + 
                     Livestock.Prod.Households.std + 
                     f(County, 
                       model = "bym", 
                       graph = ken.adj,
                       hyper = HyperBYM,
                       adjust.for.con.comp = FALSE,
                       constr = TRUE,
                       scale.model = TRUE)  +
                     f(quarter,model = "rw1",
                       scale.model = TRUE,
                       hyper = list(theta = list(prior="pc.prec", param=c(1,0.05)))), 
            family = "zeroinflatedpoisson1", 
            data = train,
            control.compute = list(config = TRUE))

NSim <- 1000
SimData.z <- inla.posterior.sample(n = NSim, result = I4.zip.sim)

# Extract Beta positions
MyParams.z <- c("(Intercept):1",  "Total.Population.std:1","Population.Density.std:1",
                "Dairy.Exotic.Cattle.std:1","Agricultural.Land.Area.std:1", "Livestock.Prod.Households.std:1")
MyID.z <- function(x){ which(rownames(SimData.z[[i]]$latent) == x) }
RowNum.Betas.z <- lapply(MyParams.z, MyID.z)
RowNum.Betas.z <- as.numeric(RowNum.Betas.z)
RowNum.Betas.z

#(rownames(SimData.z[[1]]$latent))[35000:35386]

N1 <- 290
MyParamsID.z <- paste("County:", 1:N1, sep = "")
RowNum.ID.z <- lapply(MyParamsID.z, MyID.z)
RowNum.ID.z <- unlist(RowNum.ID.z)
RowNum.ID.z


# Covariates
X <- model.matrix(~ 1 + Total.Population.std + Population.Density.std +   
                    Dairy.Exotic.Cattle.std + Agricultural.Land.Area.std + Livestock.Prod.Households.std ,
                    data = train)
X <- as.matrix(X)


library(plyr)
N  <- nrow(train)
Ysim <- matrix(nrow = N, ncol = NSim)
ZerosZIP <- vector(length = NSim) # 

# Run
library(VGAM)
for (i in 1: NSim){
  Betas <- SimData.z[[i]]$latent[RowNum.Betas.z]
  Ui <- SimData.z[[i]]$latent[RowNum.ID.z]
  Ui <- rep(Ui, 60)
  
  Pi <- SimData.z[[i]]$hyperpar[[1]]
  mu <- exp(X %*% Betas + Ui)
  Ysim[,i] <- rzipois(N,lambda = mu, pstr0 = Pi)
  ZerosZIP[i] <- sum(Ysim[,i]==0)
  
}


# Now we have 1000 simulated data sets from the model.
# What shall we do with these simulated data sets?
# We could calculate the number of zeros in each of the 1,000
# data sets.
zeros <- vector(length = NSim)
for(i in 1:NSim){
  zeros[i] <- sum(Ysim[,i] == 0)
}

summary(Ysim[,1])
#Let's plot this as a table

# Figure 24.11
par(mfrow = c(1,1), mar = c(5,5,2,2), cex.lab = 1.5)

Z <- table(zeros)
Range <- range(as.numeric(names(Z)))
SumZeros <- sum(train$cases == 0, na.rm=TRUE)
x1 <- min(Range[1], SumZeros)
x2 <- max(Range[2], SumZeros)

plot(table(zeros), 
     xlab = "How often do we have 0, 1, 2, 3, etc. number of zeros",
     ylab = "Frequency",
     xlim = c(0.998*x1,1.002*x2),
     main = "")
points(x = sum(train$cases == 0), 
       y = 0, 
       pch = 16, 
       cex = 3, 
       col = 2)
#The red dot is the number of zeros in the original data set.
#The data simulated from the NB model contains enough zeros.


#####################################################################




#####################################################################

# BOTH Results:

## edit function from ggregplot
Efxplot.Val <- function (ModelList, Sig = TRUE, StarLoc = NULL, Alpha1 = 1, 
                         Alpha2 = 1, PointOutline = F, ModelNames = NULL, VarNames = NULL, 
                         VarOrder = NULL, Intercept = TRUE, Size = 2, tips = 0.2) 
{
  require(dplyr)
  require(ggplot2)
  require(INLA)
  require(MCMCglmm)
  Graphlist <- list()
  if (!class(ModelList) == "list") {
    ModelList <- list(ModelList)
  }
  for (i in 1:length(ModelList)) {
    model <- ModelList[[i]]
    if (class(model) == "inla") {
      Graph <- as.data.frame(summary(model)$fixed)
      #Graph <- Graph[1:9,] ## edit 
      colnames(Graph)[which(colnames(Graph) %in% c("0.025quant", 
                                                   "0.975quant"))] <- c("Lower", "Upper")
      colnames(Graph)[which(colnames(Graph) %in% c("0.05quant", 
                                                   "0.95quant"))] <- c("Lower", "Upper")
      colnames(Graph)[which(colnames(Graph) %in% c("mean"))] <- c("Estimate")
    }
    if (class(model) == "MCMCglmm") {
      Graph <- as.data.frame(summary(model)$solutions)
      colnames(Graph)[1:3] <- c("Estimate", "Lower", 
                                "Upper")
    }
    Graph$Model <- i
    Graph$Factor <- rownames(Graph)
    Graphlist[[i]] <- Graph
  }
  Graph <- bind_rows(Graphlist)
  Graph$Sig <- with(Graph, ifelse(Lower * Upper > 0, "*", 
                                  ""))
  Graph$Model <- as.factor(Graph$Model)
  if (!is.null(ModelNames)) {
    levels(Graph$Model) <- ModelNames
  }
  position <- ifelse(length(unique(Graph$Model)) == 1, "none", 
                     "right")
  if (is.null(VarOrder)) 
    VarOrder <- rev(unique(Graph$Factor))
  if (is.null(VarNames)) 
    VarNames <- VarOrder
  Graph$Factor <- factor(Graph$Factor, levels = VarOrder)
  levels(Graph$Factor) <- VarNames
  Graph %<>% as.data.frame %>% filter(!is.na(Factor))
  if (!Intercept) {
    VarNames <- VarNames[!str_detect(VarNames, "ntercept")]
    Graph <- Graph %>% filter(Factor %in% VarNames)
  }
  Graph$starloc <- NA
  min <- min(Graph$Lower, na.rm = T)
  max <- max(Graph$Upper, na.rm = T)
  if (Sig == TRUE) {
    Graph$starloc <- max + (max - min)/10
  }
  if (!is.null(StarLoc)) {
    Graph$starloc <- StarLoc
  }
  Graph$Alpha <- with(Graph, ifelse(Lower * Upper > 0, Alpha1, 
                                    Alpha2))
  Graph <- Graph %>% mutate(SigAlpha = factor(as.numeric(Lower * 
                                                           Upper > 0), levels = c(1, 0)))
  if (PointOutline) {
    PointOutlineAlpha <- Alpha1
  }
  else {
    PointOutlineAlpha <- 0
  }
  ggplot(Graph, aes(x = as.factor(Factor), y = Estimate, group = Model, 
                    colour = Model, alpha = SigAlpha)) + 
    labs(y = "Linear fixed effect (mean & 95% CI)") +
    theme(text = element_text(size = 18)) +
    geom_point(position = position_dodge(w = 0.5), 
               size = Size) + geom_errorbar(position = position_dodge(w = 0.5), 
                                            aes(ymin = Lower, ymax = Upper), 
                                            size = 1.2, 
                                            width = tips) + 
    geom_hline(aes(yintercept = 0), lty = 2) + labs(x = NULL) + 
    coord_flip() + theme(legend.position = position) + geom_text(aes(label = Sig, 
                                                                     y = starloc), position = position_dodge(w = 0.5), show.legend = F) + 
    scale_alpha_manual(values = c(Alpha1, Alpha2)) + guides(alpha = "none") + 
    geom_point(colour = "black", aes(group = Model), 
               position = position_dodge(w = 0.5), size = 4, alpha = PointOutlineAlpha) + 
    geom_errorbar(aes(ymin = Lower, ymax = Upper, group = Model), 
                  width = 0.1, position = position_dodge(w = 0.5), 
                  colour = "black", alpha = PointOutlineAlpha) + 
    geom_point(position = position_dodge(w = 0.5), size = 3, 
               alpha = PointOutlineAlpha)
}

# Covariate effects
Efxplot.Val(list(I4.b,I4.zip), ModelNames = c("Occurrence Model","Incidence Model"), Intercept = FALSE)

###################################################################################


