#########################################
#Load support files and packages
library(lattice)
library(gstat)
library(raster)
library(rgdal)
library(sp)
library(spdep)
library(INLA)
library(ggmap)
library(rgeos)
library(fields)
library(mgcv)
library(ggplot2)
library(ggregplot)
library(MASS)
library('maps')
library('maptools')
library('mapdata')
data("worldHiresMapEnv")
source("C:/Kenya_Spatial_Temporal_Model/HighstatLibV10.R")
source("C:/Kenya_Spatial_Temporal_Model/HighstatLibV11.R")
source("C:/Kenya_Spatial_Temporal_Model/HighstatLibV13.R")

# Plot fixed effects
library(ggregplot)
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
      Graph <- Graph[1:2,] ## 
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
    theme(text = element_text(size = 18)) +
    geom_point(position = position_dodge(w = 0.5), 
               size = Size,
               col = "navyblue") + geom_errorbar(position = position_dodge(w = 0.5), 
                                                 aes(ymin = Lower, ymax = Upper), 
                                                 size = 1.2, 
                                                 width = tips, col = "navyblue") + 
    geom_hline(aes(yintercept = 0), lty = 2) + labs(x = NULL) + 
    coord_flip() + theme(legend.position = position) + geom_text(aes(label = Sig, 
                                                                     y = starloc), position = position_dodge(w = 0.5), show.legend = F, col="blue") + 
    scale_alpha_manual(values = c(Alpha1, Alpha2)) + guides(alpha = "none") + 
    geom_point(colour = "black", aes(group = Model), 
               position = position_dodge(w = 0.5), size = 4, alpha = PointOutlineAlpha) + 
    geom_errorbar(aes(ymin = Lower, ymax = Upper, group = Model), 
                  width = 0.1, position = position_dodge(w = 0.5), 
                  colour = "black", alpha = PointOutlineAlpha) + 
    geom_point(position = position_dodge(w = 0.5), size = 3, 
               alpha = PointOutlineAlpha)
}


#################################################
#Load the data
setwd("C:/Kenya_Spatial_Temporal_Model")
CR <- read.csv(file = "Spatial_model_data.csv", 
               header = TRUE)
names(CR)
str(CR)
colnames(CR) <- c("xcoord", "ycoord", "cases", "distance_water","elevation",
                  "calcium", "pH",   "soilwater", "rainfall",  "evi")
summary(CR) # 1 presence location has NAs; missing values for pH & calcium
dim(CR)
CR <- CR[!is.na(CR$pH),] # remove location with NA
summary(CR) # 
dim(CR)

# Change UTM meters to kilometers
CR$Xkm <- CR$xcoord / 1000
CR$Ykm <- CR$ycoord / 1000




#################################################
### DATA EXPLORATION

#* Outliers:
colnames (CR)
MyVar <- c("rainfall" , "evi" ,"distance_water" ,
           "elevation","calcium" , "pH" , "soilwater" )
Mydotplot(CR[,MyVar])


#* Collinearity:
#   - Use VIF values
#   - Use Pearson correlations
#   - Use scatterplots
# Variance inflation factors (VIF):
MyVar <- c("rainfall","evi","distance_water","elevation","calcium","pH","soilwater" )
corvif(CR[,MyVar])
Mypairs(CR[,MyVar])

# remove rainfall,  soil_pH 
MyVar <- c("evi" ,"distance_water","elevation","calcium", "soilwater")
corvif(CR[,MyVar])
Mypairs(CR[,MyVar])

# Scatter plot and Pearson correlation:
library(GGally)
ggpairs(CR[,MyVar])
# No collinearity >0.6!

# Corrplot
library(corrplot)
# corrplot 0.92 loaded
MyVar1 <- data.frame(
  Rainfall = CR$rainfall, 
  Evi = CR$evi, 
  Elevation = CR$elevation, 
  Distance.to.water = CR$distance_water, 
  Soil.calcium = CR$calcium,
  Soil.pH = CR$pH,
  Soil.water = CR$soilwater)

M = cor(MyVar1)
corrplot(M, method = 'number') # colorful number


#* Relationships:
MyVar <- c("evi" ,"distance_water","elevation","calcium", "soilwater")
MyMultipanel.ggp2(Z = CR, 
                  varx = MyVar, 
                  vary = "cases", 
                  ylab = "Anthrax Cases per Year",
                  addSmoother = TRUE,
                  addRegressionLine = FALSE,
                  addHorizontalLine = FALSE)


#* Zero-inflation:
cs<-table(CR$cases)
plot(cs/(nrow(CR)), ylab = "Proportion of counts", xlab="Counts")

#################################################
### INLA ANALYSIS

# Run a ZIP model 
#   - Without spatial autocorrelation and linear effects (ZIP-GLM)
#   - Without spatial autocorrelation and non-linear effects (ZIP-GAM)
#   - With spatial autocorrelation and non-linear effects (ZIP-GAM plus SRF)
# Compare DIC, WAIC, CPO



### Prepare the basis functions - smoothers:
#   - use an unpenalized cubic regression spline with 4 knots.


# Elevation
Lego.elevation <- smoothCon(mgcv::s(elevation, bs = "cr", k = 4, fx = TRUE),data = CR, absorb.cons = TRUE)[[1]]
Xelevation<- Lego.elevation$X 
Xelevation.df <- data.frame(Xelevation)
names(Xelevation.df) <- paste0("Xelevation", 1:ncol(Xelevation))
#lcs.elevation <- inla.make.lincombs(Xelevation.df)
lcs.elevation <- inla.make.lincombs(elevation1 = Xelevation.df[,"Xelevation1"],
                                    elevation2 = Xelevation.df[,"Xelevation2"],
                                    elevation3 = Xelevation.df[,"Xelevation3"])

# evi
Lego.evi <- smoothCon(mgcv::s(evi, bs = "cr", k = 4, fx = TRUE),data = CR, absorb.cons = TRUE)[[1]]
Xevi<- Lego.evi$X 
Xevi.df <- data.frame(Xevi)
names(Xevi.df) <- paste0("Xevi", 1:ncol(Xevi))
#lcs.evi <- inla.make.lincombs(Xevi.df)
lcs.evi <- inla.make.lincombs(evi1 = Xevi.df[,"Xevi1"],
                              evi2 = Xevi.df[,"Xevi2"],
                              evi3 = Xevi.df[,"Xevi3"])

# Give these linear combinations unique row names
names(lcs.elevation) <- paste(names(lcs.elevation), "elevation",sep = "")
names(lcs.evi) <- paste(names(lcs.evi), "evi",sep = "")
lcs <- c(lcs.elevation, lcs.evi)




### Standardize the linear covariates:
CR$distance_water.std  <- MyStd(CR$distance_water)
CR$elevation.std  <- MyStd(CR$elevation)
CR$evi.std  <- MyStd(CR$evi)
CR$soilwater.std  <- MyStd(CR$soilwater)


# Create a covariate matrix X that contains all continuous covariates
Xm <- data.frame(Intercept = rep(1,nrow(CR)),
                 distance_water.std = CR$distance_water.std,
                 elevation.std = CR$elevation.std,
                 evi.std = CR$evi.std,
                 soilwater.std = CR$soilwater,
                 elevation1 = Xelevation.df[,"Xelevation1"],
                 elevation2 = Xelevation.df[,"Xelevation2"],
                 elevation3 = Xelevation.df[,"Xelevation3"],
                 evi1 = Xevi.df[,"Xevi1"],
                 evi2 = Xevi.df[,"Xevi2"],
                 evi3 = Xevi.df[,"Xevi3"])
head(Xm)



### Create a mesh to calculate the spatial random effects:
# Estimate the distances between sites?
Loc <- cbind(CR$Xkm, CR$Ykm)
D   <- dist(Loc)
par(mfrow = c(1,1), mar = c(5,5,2,2), cex.lab = 1.5)
hist(D , freq = TRUE, main = "", xlab = "Distance between sites (km)", ylab = "Frequency")


# Make a mesh
#' What is the distance for which we would expect dependency? 
RangeGuess <- 100   # The distance for which we would expect dependency? = 100km
ConvHull <- inla.nonconvex.hull(Loc, convex = -0.05) # puts the inner boundary closer to the points.
MaxEdge <- RangeGuess / 5 # Recommended settings
mesh3 <- inla.mesh.2d(boundary = ConvHull, max.edge = c(1, 5) * MaxEdge, cutoff = MaxEdge / 5)
par(mfrow = c(1,1), mar=c(1, 1, 1, 1))
plot(mesh3, asp = 1) # Plot the mesh
points(Loc, col = "red", pch = 16, cex = 1)
mesh3$n # 


# Projector matrix
A3  <- inla.spde.make.A(mesh3, loc = Loc) # Define the weighting factors a_ik


# Define the SPDE
spde3 <- inla.spde2.pcmatern(mesh3, prior.range = c(100, 0.1), prior.sigma = c(0.1, 0.0001))


# Define w (the spatial field)
w3.index <- inla.spde.make.index(name = 'w', n.spde  = spde3$n.spde) # Define the spatial field


#' Make a stack for ZIP GAM with and without spatial correlation
StackGAM <- inla.stack(
  tag = "Fit",
  data = list(y = CR$cases),  
  A = list(1, 1, A3),                 
  effects = list(
    Intercept = rep(1, nrow(CR)),  #Intercept
    Xm        = Xm[,-1],           #All covariates in the model
    w         = w3.index))         #SRF




#Run the INLA models
# Execute the ZIP GLM without autocorrelation
I1 <- inla(formula(y ~ -1 + Intercept + distance_water.std + elevation.std + evi.std ),
           family = "zeroinflatedpoisson1", 
           control.inla = list(strategy = "gaussian"), #For faster calculations
           data = inla.stack.data(StackGAM),
           control.compute = list(dic = TRUE, waic = TRUE),
           control.predictor = list(A = inla.stack.A(StackGAM)))

# Execute the ZIP GAM without autocorrelation
I2 <- inla(formula(y ~ -1 + Intercept + distance_water.std + 
                     elevation1 + elevation2 + elevation3 +
                     evi1 + evi2 + evi3),
           family = "zeroinflatedpoisson1", 
           control.inla = list(strategy = "gaussian"), #For faster calculations
           data = inla.stack.data(StackGAM),
           control.compute = list(dic = TRUE, waic = TRUE),
           control.predictor = list(A = inla.stack.A(StackGAM)))

# Execute the ZIP GAM with spatial correlation
I3 <- inla(formula(y ~ -1 + Intercept + distance_water.std + 
                     elevation1 + elevation2 + elevation3 +
                     evi1 + evi2 + evi3 +
                     f(w, model = spde3)),
           lincomb = c(lcs.elevation, lcs.evi),  
           family = "zeroinflatedpoisson1", 
           control.inla = list(strategy = "gaussian"),  #For faster computation
           data = inla.stack.data(StackGAM),
           control.compute = list(dic = TRUE, waic = TRUE),
           control.predictor = list(A = inla.stack.A(StackGAM)))


#################################################

### Model Selection:

dic   <- c(I1$dic$dic, I2$dic$dic, I3$dic$dic)
waic  <- c(I1$waic$waic, I2$waic$waic, I3$waic$waic) 
Result <- cbind(dic, waic)
rownames(Result) <- c("ZIP GLM","ZIP GAM","ZIP GAM + SRF")
Result
## Save
setwd("C:/Kenya_Spatial_Temporal_Model")
class(Result)
Result <- as.data.frame(Result)
write.csv(Result, file = "Model Validation Results - Spatial ZIP.csv")       




#################################################

### Model Validation:



## Check for residual spatial dependency:

# Calculate the Pearson residuals
N   <- nrow(CR)
mu3 <- I3$summary.fitted.values[1:N,"mean"]
Pi  <- I3$summary.hyperpar["zero-probability parameter for zero-inflated poisson_1","mean"]
#' Recall that these are the mean and variance of a ZIP:
#'   E[anthrax_ti]   = (1 - Pi) * mu_ti
#'   var[anthrax_ti] = (1 - Pi) * (mu_ti + Pi * mu_ti^2)
ExpY <- (1-Pi) * mu3
VarY <- (1-Pi) * (mu3 + Pi * mu3^2)
# And that gives the following Pearson residuals
E3 <- (CR$cases - ExpY) / sqrt(VarY)
# A quick way to assess over-dispersion:
Npar <- nrow(I3$summary.fixed) + 1 + 1 #1 sigma and 1 Pi
sum(E3^2) / (N-Npar)

# Check for spatial patterns in the Pearson residuals with a variogram:
MyData <- data.frame(E3 = E3, Xkm = CR$Xkm, Ykm = CR$Ykm)
coordinates(MyData) <- c("Xkm", "Ykm")
V5 <- variogram(E3 ~ 1, MyData, cressie = TRUE, cutoff = 300)
ggplot() +
  geom_point(data = V5, aes(x = dist, y = gamma)) +
  geom_smooth(data = V5, aes(x = dist, y = gamma),col = "red") +
  xlab("Distance") + ylab("Semi-variogram") +
  theme(text = element_text(size = 15)) +
  theme(legend.position="none") +
  ylim(0, 1) 

# Plot Pearson residuals vs fitted values"
par(mfrow = c(1,1), mar = c(5,5,2,2), cex.lab = 1.5)
plot(x = ExpY, y = E3, xlab = "Fitted values", ylab = "Pearson residuals")

# Plot fitted values vs observed values
par(mfrow = c(1,1), mar = c(5,5,2,2), cex.lab = 1.5)
plot(x = ExpY, y = CR$cases, xlim = c(0, 85), ylim = c(0, 85),
     xlab = "Fitted values model I3",   ylab = "Observed data")





## Check for non-linear patterns/effects in the covariales:

# Plot residuals vs every covariate       
MyVar <- c("evi", "distance_water", "elevation")
CR$E3 <- E3
MyMultipanel.ggp2(Z = CR, 
                  varx = MyVar, 
                  vary = "E3", 
                  ylab = "Pearson residuals",
                  addSmoother = TRUE,
                  addRegressionLine = FALSE,
                  addHorizontalLine = TRUE)

# Difficult to assess. Use mgcv for this
T0 <- gam(E3 ~ 1, data = CR)
Tevi         <- gam(E3 ~ s(evi), data = CR)
Tdistw       <- gam(E3 ~ s(distance_water), data = CR)
Telev        <- gam(E3 ~ s(elevation), data = CR)
Tswater      <- gam(E3 ~ s(soilwater), data = CR)
AIC(T0, Tevi , Tdistw, Telev, Tswater)

# plot
par(mfrow = c(2,2))
plot(Tevi)
abline(h = 0, lty = 2)
#text(0.5, 1.55, "B", cex = 1.5)

plot(Tdistw)
abline(h = 0, lty = 2)
#text(0.5, 1.55, "B", cex = 1.5)

plot(Telev)
abline(h = 0, lty = 2)
#text(0.5, 1.55, "B", cex = 1.5)

#plot(Tcal)
#abline(h = 0, lty = 2)
#text(0.5, 1.55, "B", cex = 1.5)

plot(Tswater)
abline(h = 0, lty = 2)
text(0.5, 1.55, "B", cex = 1.5)

par(mfrow = c(1,1))

## Do the credible intervals of each covariate contain zero across the entire gradient
## of the covariates? Yes to all 





## Check for model's ability to handle zeros:

# Simulation study: Can the model handle zeros?
# Simulate data sets from the model to determine whether
# the model can cope with zero inflation (and check overdispersion
# in a more Bayesian way).

# ZIP model with spatial-temporal correlation
## run model
I3.sim <- inla(formula(y ~ -1 + Intercept + distance_water.std + soilwater.std +
                         elevation1 + elevation2 + elevation3 +
                         evi1 + evi2 + evi3 +
                         f(w, model = spde3)), 
               family = "zeroinflatedpoisson1", 
               control.inla = list(strategy = "gaussian"),  #For faster computation
               data = inla.stack.data(StackGAM),
               control.compute   = list(config = TRUE),
               control.predictor = list(A = inla.stack.A(StackGAM)))
NSim <- 1000
SimData <- inla.posterior.sample(n = NSim, result = I3.sim)

# Extract Beta positions
MyParams <- c("Intercept:1", "distance_water.std:1", "soilwater.std:1",
              "elevation1:1","elevation2:1","elevation3:1",
              "evi1:1", "evi2:1", "evi3:1")
MyID <- function(x){ which(rownames(SimData[[i]]$latent) == x) }
RowNum.Betas <- lapply(MyParams, MyID)
RowNum.Betas <- as.numeric(RowNum.Betas)
RowNum.Betas

N1 <- mesh3$n # - 1
MyParams <- paste("w:", 0:N1, sep = "")
RowNum.w <- lapply(MyParams, MyID)
RowNum.w <- unlist(RowNum.w)
RowNum.w

# Start a loop to extract betas, calculate
# the fitted values and simulate count data from 
# the model.
N  <- nrow(CR)
Ysim <- matrix(nrow = N, ncol = NSim)
mu.i <- matrix(nrow = N, ncol = NSim)
Xmat <- as.matrix(data.frame(Intercept = rep(1,N),
                             distance_water.std = CR$distance_water.std,
                             soilwater.std      = CR$soilwater.std,
                             elevation1 = Xelevation.df[,"Xelevation1"],
                             elevation2 = Xelevation.df[,"Xelevation2"],
                             elevation3 = Xelevation.df[,"Xelevation3"],
                             evi1 = Xevi.df[,"Xevi1"],
                             evi2 = Xevi.df[,"Xevi2"],
                             evi3 = Xevi.df[,"Xevi3"]))
A    <- as.matrix(A3)

library(VGAM)
for (i in 1: NSim){
  Betas <- SimData[[i]]$latent[RowNum.Betas]
  wi    <- SimData[[i]]$latent[RowNum.w] 
  Pi <- SimData[[i]]$hyperpar[[1]]
  
  eta   <- Xmat %*% Betas + A %*% wi
  mu.i[,i]    <- exp(eta)
  Ysim[,i]    <- rzipois(n = nrow(CR), lambda = mu.i[,i], pstr0 = Pi)
}

table(Ysim[,1])
table(Ysim[,2])
table(Ysim[,3])

# Now we have 1000 simulated data sets from the model.
# What shall we do with these simulated data sets?
# We could calculate the number of zeros in each of the 1,000
# data sets.
zeros <- vector(length = NSim)
for(i in 1:NSim){
  zeros[i] <- sum(Ysim[,i] == 0)
}
table(zeros)

#Let's plot this as a table
par(mar = c(5,5,2,2), cex.lab = 1.5)
plot(table(zeros), 
     #axes = FALSE,
     xlab = "How often do we have 0, 1, 2, 3, etc. number of zeros",
     ylab = "Number of zeros in 1000 simulated data sets",
     xlim = c(10,200),
     main = "Simulation results")
points(x = sum(CR$cases == 0), 
       y = 0, 
       pch = 16, 
       cex = 5, 
       col = 2)
#The red dot is the number of zeros in the original data set.
sum(sum(CR$cases == 0) > zeros) / 1000 # simulations with zero counts less than observed data
# [1] 0.805




###################################################

### Model assessment of performance:

# Cut-off/threshold value to evaluate model performance against withheld data 
# after predictions:

# minimum fitted value of the training presence data
mu <- I3$summary.fitted.values$mean[1:94] 
Pi  <- I3$summary.hyperpar["zero-probability parameter for zero-inflated poisson_1","mean"]
min.tp <- min((1-Pi) * mu)
min.tp    # 0.2052692





###################################################

### Model interpretation:

# Fixed effects best model:
Beta <- I3$summary.fixed[1:9, c("mean", "0.025quant", "0.975quant")] 
print(Beta, digits = 3)
# compare fixed effects/betas:
library(ggregplot)
Efxplot(list(I1, I2, I3), ModelNames = c("ZIP GLM","ZIP GAM", "ZIP GAM + SRF"))
Efxplot.Val(list(I3), ModelNames = c("ZIP GAM + SRF"))





# Smoothers for the best model:
# Fitted values: We first calculate the Pi = probability of false zeros
Pi <- I3$summary.hyper[1,"mean"]
Pi
# Fitted model:
## count_i = ZIP(mu, Pi)
## E(count_i) = (1-Pi) * mu
## Var(count_i) = (1-Pi) * (mu + Pi * mu^2)
## log(mu) = intercept + covariate effects
mu <- I3$summary.fitted.values[,"mean"]
ExpY <- (1-Pi) * mu

Ns <- nrow(CR) # Extract the relevant columns of I3$summary.lincomb.derived
# elevation smoother:
f.elevation    <- (1-Pi) * (exp(I3$summary.lincomb.derived[1:Ns + 0 * Ns, "mean"]))
SeLo.elevation <- (1-Pi) * (exp(I3$summary.lincomb.derived[1:Ns + 0 * Ns,"0.025quant"])) 
SeUp.elevation <- (1-Pi) * (exp(I3$summary.lincomb.derived[1:Ns + 0 * Ns,"0.975quant"]))
# evi smoother: rows 1+Ns to Ns+Ns
f.evi    <- (1-Pi) * (exp(I3$summary.lincomb.derived[1:Ns + 1 * Ns, "mean"])) 
SeLo.evi <- (1-Pi) * (exp(I3$summary.lincomb.derived[1:Ns + 1 * Ns,"0.025quant"])) 
SeUp.evi <- (1-Pi) * (exp(I3$summary.lincomb.derived[1:Ns + 1 * Ns,"0.975quant"]))

MyData <- data.frame(
  mu    = c(f.elevation, f.evi), 
  SeUp  = c(SeUp.elevation, SeUp.evi), 
  SeLo  = c(SeLo.elevation, SeLo.evi), 
  Xaxis = c(CR$elevation, CR$evi),
  ID    = factor(rep(c("Elevation smoother", "EVI smoother"), 
                     each = nrow(CR))))

#' Plot the two smoothers
ggplot() +
  xlab("Covariate") + ylab("Smoother") +
  theme(text = element_text(size = 15)) +
  geom_line(data = MyData, aes(x = Xaxis, y = mu)) +
  geom_ribbon(data = MyData, aes(x = Xaxis, ymax = SeUp, ymin = SeLo),alpha = 0.6) +
  facet_wrap(~ID, scales = "free", ncol = 2)     

XPos <- c(CR$elevation, CR$evi)
XID  <- rep(c("Elevation smoother", "EVI smoother"), each = nrow(CR) )
MyData2 <- data.frame(Xaxis = XPos,ID    = factor(XID))
ggplot() +
  xlab("Covariate") + ylab("Fitted Anthrax Cases") +
  theme(text = element_text(size = 15)) +
  geom_line(data = MyData, aes(x = Xaxis, y = mu)) +
  geom_ribbon(data = MyData, aes(x = Xaxis, ymax = SeUp, ymin = SeLo),alpha = 0.6) +
  facet_wrap(~ID, scales = "free", ncol = 2) +
  geom_rug(data = MyData2, aes(x = Xaxis), size = 0.05)





# Spatial random effects:
w.pm <- I3$summary.random$w$mean  
w.sd <- I3$summary.random$w$sd
# Plot mean
PlotField2(field = w.pm, mesh = mesh3, xlim = range(mesh3$loc[,1]), ylim = range(mesh3$loc[,2]), MyMain = "Posterior Mean SRF")
points(x = Loc[CR$cases>0,1], y = Loc[CR$cases>0,2], cex = 0.1, col = "black", pch = 16)
# Plot sd
PlotField2(field = w.sd, mesh = mesh3, xlim = range(mesh3$loc[,1]), ylim = range(mesh3$loc[,2]), MyMain = "Posterior SD SRF")
points(x = Loc[CR$cases>0,1],y = Loc[CR$cases>0,2], cex = 0.1, col = "black", pch = 16)

## Save as raster
points.em <- mesh3$loc
stepsize <- 1                           # This is given in coordinates unit (in this case this is straightforward and correspond to 150m)
east.range <- diff(range(points.em[,1]))  # calculate the length of the Easting range
north.range <- diff(range(points.em[,2])) # calculate the length of the Northing range
nxy <- round(c(east.range, north.range)/stepsize)  # Calculate the number of cells in the x and y ranges

# Project the spatial field on the mesh vertices using the inla.mesh.projector() function
projgrid <- inla.mesh.projector(mesh3, xlim = range(points.em[,1]), ylim = range(points.em[,2]), dims = nxy)
xmean <- inla.mesh.project(projgrid,I3$summary.random$w$mean)
xsd <- inla.mesh.project(projgrid, I3$summary.random$w$sd)

#We need to create spatial objects for the mean and variance of the GRF.
require(raster)
xmean2 <- t(xmean)
xmean3 <- xmean2[rev(1:length(xmean2[,1])),]
xmean_ras <- raster(xmean3,
                    xmn = range(projgrid$x)[1], xmx = range(projgrid$x)[2],
                    ymn = range(projgrid$y)[1], ymx = range(projgrid$y)[2],
                    crs = CRS("+proj=utm +zone=36N +datum=WGS84 +units=km +ellps=WGS84 +towgs84=0,0,0"))
xsd2 <- t(xsd)
xsd3 <- xsd2[rev(1:length(xsd2[,1])),]
xsd_ras <- raster(xsd3,
                  xmn = range(projgrid$x)[1], xmx =range(projgrid$x)[2],
                  ymn = range(projgrid$y)[1], ymx =range(projgrid$y)[2],
                  crs = CRS("+proj=utm +zone=36N +datum=WGS84 +units=km +ellps=WGS84 +towgs84=0,0,0"))

library(RColorBrewer)
my.palette.post <- rev(brewer.pal(n = 9, name = "YlGnBu"))
my.palette.var <- brewer.pal(n = 9, name = "BuPu")

par(mfrow = c(1,1), mar = c(2,2, 1,1))
plot(xmean_ras, asp = 1, col = my.palette.post)
points(x = Loc[,1], y = Loc[,2], cex = 0.5, col = "black", pch = 16)
plot(xsd_ras, asp = 1, col = my.palette.var)
points(x = Loc[,1],y = Loc[,2], cex = 0.5, col = "black", pch = 16)

# save
library(raster)
setwd("C:/Kenya_Spatial_Temporal_Model/Spatial Model/Prediction Rasters 1")
writeRaster(xmean_ras, filename="posterior mean of SRF_spatial model.tif", format="GTiff", overwrite=TRUE)
writeRaster(xsd_ras, filename="posterior sd of SRF_spatial model.tif", format="GTiff", overwrite=TRUE)





#################################################

### Model prediction: 

# First, predict across Kenya then check the omission rate against wildlife data

#Set the working directory and load the data
setwd("C:/Kenya_Spatial_Temporal_Model")
pred.data <- read.csv(file = "1km_grid_kenya_coords_plusvars.csv", 
                      header = TRUE,
                      na.strings = "NA",
                      stringsAsFactors = TRUE,
                      dec = ".")
head(pred.data)
pred.data$cases <- rep(NA, nrow(pred.data))
colnames(pred.data) <- c("xcoord", "ycoord", "distance_water","elevation",
                  "calcium", "pH",   "soilwater", "rainfall",  "evi", "cases")
head(pred.data)

# remove NAs from prediction data
summary(pred.data)
dim(pred.data)
pred.nas <- pred.data[is.na(pred.data$evi),]
#pred.nas2 <- pred.data[is.na(pred.data$soilwater),]
#pred.nas <- rbind (pred.nas1, pred.nas2)

pred.data <- pred.data[!is.na(pred.data$evi),]
#pred.data <- pred.data[!is.na(pred.data$soilwater),]
summary(pred.data)
dim(pred.data) # 594782 

# Grid coords that were removed due to NAs in covariates: Add these to final predicted data
colnames(pred.nas)
pred.nas <- pred.nas[,c(1,2, 10)]
head(pred.nas)
colnames(pred.nas) <- c("longitude", "latitude", "mean")
head(pred.nas)

## Clamping: 
## Set the prediction variables outside the calibration range of the 
## models to be equal to the maximum or minimum values of the training dataset
##

## Distance to water
pred.data$distance_water  <- ifelse(pred.data$distance_water < (min(CR$distance_water)), 
                                     (min(CR$distance_water)), pred.data$distance_water)
pred.data$distance_water  <- ifelse(pred.data$distance_water > (max(CR$distance_water)), 
                                     (max(CR$distance_water)), pred.data$distance_water)

## Elevation
pred.data$elevation  <- ifelse(pred.data$elevation < (min(CR$elevation)), 
                                (min(CR$elevation)), pred.data$elevation)
pred.data$elevation  <- ifelse(pred.data$elevation > (max(CR$elevation)), 
                          (max(CR$elevation)), pred.data$elevation)

## Evi
pred.data$evi  <- ifelse(pred.data$evi < (min(CR$evi)), 
                          (min(CR$evi)), pred.data$evi)
pred.data$evi  <- ifelse(pred.data$evi > (max(CR$evi)), 
                          (max(CR$evi)), pred.data$evi)

## Soil water
pred.data$soilwater  <- ifelse(pred.data$soilwater < (min(CR$soilwater)), 
                                (min(CR$soilwater)), pred.data$soilwater)
pred.data$soilwater  <- ifelse(pred.data$soilwater > (max(CR$soilwater)), 
                                (max(CR$soilwater)), pred.data$soilwater)



# Divide prediction data into small chunks to prevent INLA from crashing
list.pred <- list()
list.pred[[1]] <- pred.data[1:10000,]
list.pred[[2]] <- pred.data[10001:20000,]
list.pred[[3]] <- pred.data[20001:30000,]
list.pred[[4]] <- pred.data[30001:40000,]
list.pred[[5]] <- pred.data[40001:50000,]
list.pred[[6]] <- pred.data[50001:60000,]
list.pred[[7]] <- pred.data[60001:70000,]
list.pred[[8]] <- pred.data[70001:80000,]
list.pred[[9]] <- pred.data[80001:90000,]
list.pred[[10]] <- pred.data[90001:100000,]
list.pred[[11]] <- pred.data[100001:110000,]
list.pred[[12]] <- pred.data[110001:120000,]
list.pred[[13]] <- pred.data[120001:130000,]
list.pred[[14]] <- pred.data[130001:140000,]
list.pred[[15]] <- pred.data[140001:150000,]
list.pred[[16]] <- pred.data[150001:160000,]
list.pred[[17]] <- pred.data[160001:170000,]
list.pred[[18]] <- pred.data[170001:180000,]
list.pred[[19]] <- pred.data[180001:190000,]
list.pred[[20]] <- pred.data[190001:200000,]
list.pred[[21]] <- pred.data[200001:210000,]
list.pred[[22]] <- pred.data[210001:220000,]
list.pred[[23]] <- pred.data[220001:230000,]
list.pred[[24]] <- pred.data[230001:240000,]
list.pred[[25]] <- pred.data[240001:250000,]
list.pred[[26]] <- pred.data[250001:260000,]
list.pred[[27]] <- pred.data[260001:270000,]
list.pred[[28]] <- pred.data[270001:280000,]
list.pred[[29]] <- pred.data[280001:290000,]
list.pred[[30]] <- pred.data[290001:300000,]
list.pred[[31]] <- pred.data[300001:310000,]
list.pred[[32]] <- pred.data[310001:320000,]
list.pred[[33]] <- pred.data[320001:330000,]
list.pred[[34]] <- pred.data[330001:340000,]
list.pred[[35]] <- pred.data[340001:350000,]
list.pred[[36]] <- pred.data[350001:360000,]
list.pred[[37]] <- pred.data[360001:370000,]
list.pred[[38]] <- pred.data[370001:380000,]
list.pred[[39]] <- pred.data[380001:390000,]
list.pred[[40]] <- pred.data[390001:400000,]
list.pred[[41]] <- pred.data[400001:410000,]
list.pred[[42]] <- pred.data[410001:420000,]
list.pred[[43]] <- pred.data[420001:430000,]
list.pred[[44]] <- pred.data[430001:440000,]
list.pred[[45]] <- pred.data[440001:450000,]
list.pred[[46]] <- pred.data[450001:460000,]
list.pred[[47]] <- pred.data[460001:470000,]
list.pred[[48]] <- pred.data[470001:480000,]
list.pred[[49]] <- pred.data[480001:490000,]
list.pred[[50]] <- pred.data[490001:500000,]
list.pred[[51]] <- pred.data[500001:510000,]
list.pred[[52]] <- pred.data[510001:520000,]
list.pred[[53]] <- pred.data[520001:530000,]
list.pred[[54]] <- pred.data[530001:540000,]
list.pred[[55]] <- pred.data[540001:550000,]
list.pred[[56]] <- pred.data[550001:560000,]
list.pred[[57]] <- pred.data[560001:570000,]
list.pred[[58]] <- pred.data[570001:580000,]
list.pred[[59]] <- pred.data[580001:594782,]


list.pred.result <- list()

for(i in 1:59){
  pred.data1 <- list.pred[[i]]
  
  # Create columns for UTM in Km
  pred.data1$Xkm   <- pred.data1$xcoord #/ 1000 
  pred.data1$Ykm   <- pred.data1$ycoord #/ 1000
  
  
  # Standardize the linear covariates:
  pred.data1$distance_water.std  <- MyStd(pred.data1$distance_water)
  pred.data1$soilwater.std       <- MyStd(pred.data1$soilwater)
  
  
  # Get the basis functions for the  smoothers.
  Lego.elevation.pred        <- PredictMat(Lego.elevation, pred.data1)
  Lego.evi.pred              <- PredictMat(Lego.evi, pred.data1)
  # Give the basis functions unique names.
  colnames(Lego.elevation.pred)      <- paste("elevation", 1:(ncol(Xelevation)), sep = "")
  colnames(Lego.evi.pred)    <- paste("evi", 1:(ncol(Xevi)), sep = "")
  
  
  # dataframe of covariates
  Nx=nrow(pred.data1)
  Xm.pred <- data.frame(Intercept = rep(1,Nx),
                        distance_water.std = pred.data1$distance_water.std,
                        soilwater.std = pred.data1$soilwater.std,
                        elevation1 = Lego.elevation.pred[,"elevation1"],
                        elevation2 = Lego.elevation.pred[,"elevation2"],
                        elevation3 = Lego.elevation.pred[,"elevation3"],
                        evi1 = Lego.evi.pred[,"evi1"],
                        evi2 = Lego.evi.pred[,"evi2"],
                        evi3 = Lego.evi.pred[,"evi3"])
  
  # Here is the stack for the observed data
  N <- nrow(CR)
  StackGAM <- inla.stack(
    tag = "Fit",
    data = list(y = CR$cases),  
    A = list(1, 1, A3),                 
    effects = list(  
      Intercept = rep(1, nrow(CR)),
      Xm        = Xm[,-1],
      w         = w3.index))
  
 
  
  # The locations for which we want to do a prediction  
  LocPred <- pred.data1[,c("Xkm","Ykm")]
  LocPred <- as.matrix(LocPred)
  
  # Get the A matrix for this point:
  A.pred16 <- inla.spde.make.A(mesh = mesh3, loc = LocPred)
  dim(A.pred16)
  # Instead of predicting for only 1 spatial point (which can be anywhere in
  # the mesh), you can also select a series of points for which you want to
  # predict. Just add them to LocPred. And these do not have to be actual 
  # sampling locations, but can be anywhere in the mesh!
  # Section 6.8.1 in Blangiardo and Cameletti does this for a large grid of 
  # spatial values, which are then plotted.
  
  # Here is the stack for the prediction Note that NA for cases.
  StackPred <- inla.stack(
    tag = "Predict",
    data = list(y = NA),  
    A = list(1, 1, A.pred16),               
    effects = list(  
      Intercept   = rep(1, Nx),
      Xm.pred     = Xm.pred[,-1],    
      w           = w3.index))
  
  
    # We can combine the two stacks.         
  All.stacks <- inla.stack(StackGAM, StackPred)	              
  
  # ZIP GAM with spatial correlation 
  I1.Pred <- inla(y ~ -1 + Intercept + 
                    distance_water.std + #soilwater.std +
                    elevation1 + elevation2 + elevation3 +
                    evi1 + evi2 + evi3  +
                    f(w, model = spde3),
                  family = "zeroinflatedpoisson1",
                  data = inla.stack.data(All.stacks),
                  control.predictor = list(link = 1, A = inla.stack.A(All.stacks)))
  
  # Extract
  index.Fit  <- inla.stack.index(All.stacks, tag = "Fit")$data
  index.Pred <- inla.stack.index(All.stacks, tag = "Predict")$data
  # And we can extract the correct rows     
  Fit  <- I1.Pred$summary.fitted.values[index.Fit, c(1,2,3,5)]   #fitted values
  Pred <- I1.Pred$summary.fitted.values[index.Pred, c(1,2,3,5)]  #predicted values
  
  Pi  <- I1.Pred$summary.hyperpar["zero-probability parameter for zero-inflated poisson_1","mean"]
  Pred <- (1-Pi) * Pred
  
  list.pred.result[[i]] <- Pred
}

Pred <- rbind(list.pred.result[[1]],
              list.pred.result[[2]],
              list.pred.result[[3]],
              list.pred.result[[4]],
              list.pred.result[[5]],
              list.pred.result[[6]],
              list.pred.result[[7]],
              list.pred.result[[8]],
              list.pred.result[[9]],
              list.pred.result[[10]],
              list.pred.result[[11]],
              list.pred.result[[12]],
              list.pred.result[[13]],
              list.pred.result[[14]],
              list.pred.result[[15]],
              list.pred.result[[16]],
              list.pred.result[[17]],
              list.pred.result[[18]],
              list.pred.result[[19]],
              list.pred.result[[20]],
              list.pred.result[[21]],
              list.pred.result[[22]],
              list.pred.result[[23]],
              list.pred.result[[24]],
              list.pred.result[[25]],
              list.pred.result[[26]],
              list.pred.result[[27]],
              list.pred.result[[28]],
              list.pred.result[[29]],
              list.pred.result[[30]],
              list.pred.result[[31]],
              list.pred.result[[32]],
              list.pred.result[[33]],
              list.pred.result[[34]],
              list.pred.result[[35]],
              list.pred.result[[36]],
              list.pred.result[[37]],
              list.pred.result[[38]],
              list.pred.result[[39]],
              list.pred.result[[40]],
              list.pred.result[[41]],
              list.pred.result[[42]],
              list.pred.result[[43]],
              list.pred.result[[44]],
              list.pred.result[[45]],
              list.pred.result[[46]],
              list.pred.result[[47]],
              list.pred.result[[48]],
              list.pred.result[[49]],
              list.pred.result[[50]],
              list.pred.result[[51]],
              list.pred.result[[52]],
              list.pred.result[[53]],
              list.pred.result[[54]],
              list.pred.result[[55]],
              list.pred.result[[56]],
              list.pred.result[[57]],
              list.pred.result[[58]],
              list.pred.result[[59]]
)
pred.result <- cbind(pred.data, Pred)
head(pred.result)


### Convert to raster format
summary(pred.result)
colnames(pred.result)

pred.result$longitude <- pred.result$xcoord
pred.result$latitude <- pred.result$ycoord
colnames(pred.result)

xyz.data.mean <- pred.result[,c(15,16,11)]
xyz.data.sd <- pred.result[,c(15,16,12)]
xyz.data.low <- pred.result[,c(15,16,13)]
xyz.data.up <- pred.result[,c(15,16,14)]

head(xyz.data.mean)
head(pred.nas)

# Rbind
xyz.data.mean <- rbind(xyz.data.mean, pred.nas)

sd.nas <- pred.nas
colnames(sd.nas) <- c("longitude", "latitude", "sd")
xyz.data.sd <- rbind(xyz.data.sd, sd.nas)

low.nas <- pred.nas
colnames(low.nas) <- c("longitude", "latitude", "0.025quant")
xyz.data.low <- rbind(xyz.data.low, low.nas)

up.nas <- pred.nas
colnames(up.nas) <- c("longitude", "latitude", "0.975quant")
xyz.data.up <- rbind(xyz.data.up, up.nas)

#plot(xyz.data.mean$longitude, xyz.data.mean$latitude)
summary(xyz.data.mean)
summary(xyz.data.sd)

## Save CSVs:
setwd("C:/Kenya_Spatial_Temporal_Model/Prediction Rasters 1")
write.csv(xyz.data.mean, file = ("Anthrax_Incidence_mean.csv"))
write.csv(xyz.data.sd, file = ("Anthrax_Incidence_sd.csv"))
write.csv(xyz.data.low, file = ("Anthrax_Incidence_lower.csv"))
write.csv(xyz.data.up, file = ("Anthrax_Incidence_upper.csv"))

##########################################################################


