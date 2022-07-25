library(MASS)
library(corrplot)
library(tidyverse)
library(ggplot2)
library(GGally)
library(cowplot)
library(stringr)
library(psych)
library(Rmisc)
library(maxLik)

setwd("/home/ULreg/")

## TRAIN DATA SET
UL <- read.csv("ULReg_Train.csv", sep="")
UL[order(as.Date(UL$momento, format="%Y-%m-%d")),]
time_h=lubridate::hour(UL$momento)
time_w=lubridate::week(UL$momento)
time_m=lubridate::month(UL$momento)
time_m=as.factor(time_m)

data=cbind(UL,time_h,time_w,time_m)

colnames(data)=c("Wind",             
                 "time",
                 "Humidity",
                 "Radiation",
                 "Temperature",              
                 "time_h",             
                 "time_w",
                 "MONTH")

for(i in 1:length(data$Humidity)){
  if(data$Humidity[i]<10) data$Humidity[i]=data$Humidity[i]*10
}

# Plot
p1=ggparcoord(data[data$time_h==10,],columns = c(1,3,4,5), groupColumn = 8, 
              alphaLines = 0.03,scale="globalminmax") + 
  xlab("") + theme(
    title=element_text(size=17),
    text=element_text(size=17),
    legend.position="none")

p2=ggparcoord(data[data$time_h==10,],columns = c(1,3,4,5), groupColumn = 8,
              alphaLines = 0.3,scale="globalminmax") + 
  xlab("") + theme(
    axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 0.5,"top"),
    axis.title.y = element_blank(),
    text=element_text(size=17),
    title=element_text(size=17)) + ylim(c(-30,1610))

plot_grid(p1, p2, labels = c('', ''), rel_widths = c(0.7, 0.3))

ggsave("7a.tiff", units="in", width=12, height=5, dpi=300, compression = 'lzw')
dev.off()

## COORPLOT (ALL COVARIABLES FROM THE CHILEAN GOVERNMENT)
UL_cov <- read.csv("FULL.csv", sep="")
M<-cor(na.omit(UL_cov))
head(round(M,2))

cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}
# matrix of the p-value of the correlation
p.mat <- cor.mtest(na.omit(UL_cov))
head(p.mat[, 1:5])

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))

tiff("5.tiff", compression = 'lzw',width=600, height=480)
corrplot(M, method="color", col=col(200),  
         type="upper", order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         # Combine with significance
         p.mat = p.mat, sig.level = 0.01, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE )

dev.off()

## MULTICOLLINEARITY
model_all=lm(hr~.,data=na.omit(UL_cov))
vif_values <- car::vif(model_all)           #create vector of VIF values
barplot(vif_values,las = 2, main = "VIF Values", horiz = TRUE, col = "steelblue") #create horizontal bar chart to display each VIF value
abline(v = 5, lwd = 3, lty = 2)    #add vertical line at 5 as after 5 there is severe correlation

var_inv <- ginv(M)                                       # independent variables inverse correlation matrix 
colnames(var_inv) <- colnames(UL_cov)                      # rename the row names and column names
rownames(var_inv) <- colnames(UL_cov)
corrplot(var_inv,method='number',is.corr = F)              # visualize the multicollinearity

require(dplyr)
## CIRCULAR BARPLOT (THREE MOST IMPORTANT VARIABLES)
plot_df <- data %>%
  group_by(MONTH) %>%
  summarise(
    sum_length = mean(Radiation)/10,
    mean_gain = min(as.numeric(Humidity)),
    max_gain = max(as.numeric(Humidity)),
    Humid. = mean(as.numeric(Humidity))
  ) %>%
  mutate(mean_gain = round(mean_gain, digits = 0))

plt <- ggplot(plot_df) +
  # Make custom panel grid
  geom_hline(
    aes(yintercept = y), 
    data.frame(y = seq(0,60,20) ),
    color = "lightgrey"
  ) + 
  geom_col(
    aes(
      x = reorder(str_wrap(MONTH, 5), 1:12),
      y = sum_length,
      fill = Humid.
    ),
    position = "dodge2",
    show.legend = TRUE,
    alpha = .9
  ) +
  
  # Lollipop shaft for min gain per region
  geom_segment(
    aes(
      x = reorder(str_wrap(MONTH, 5), 1:12),
      y = 0,
      xend = reorder(str_wrap(MONTH, 5), 1:12),
      yend = 60
    ),
    linetype = "dashed",
    color = "gray12"
  ) + 
  
  # Add dots to represent the mean gain
  geom_point(
    aes(
      x = reorder(str_wrap(MONTH, 5), 1:12),
      y = mean_gain
    ),
    size = 3,
    color = "red"
  ) +
  theme(
    title=element_text(size=19),
    text=element_text(size=17))+
  xlab("")+ylab("Radiation (°/10)")+
  # Make it circular!
  coord_polar()

plt
ggsave("6c.tiff", units="in", width=12, height=5, dpi=300, compression = 'lzw')
dev.off()

## Unit-Lindley REGESSION
logvero.L=function(beta,yi,x) {
  p= exp(beta%*%t(x))/(1+exp(beta%*%t(x)))
  logvero= sum(2*log(1- p) - log(p) - 3*log(1-yi) - (yi*(1 - p)/(p*(1 - yi))))
}
set.seed(123456)
est = maxLik(logLik=logvero.L, start=rep(0,4), 
             y= as.numeric(UL$hr/100), x = cbind(rep(1, dim(UL)[1]), 
             UL[c("ffInst", "td","radiacionGlobalInst")]), method="BFGS")

summary(est)

## PHASE II - VISUALIZATION
UL2 <- read.csv("ULReg_Test.csv", sep="")
which(is.na(UL2$hr))
UL2=na.omit(UL2)
colnames(UL2)=c("Wind",             
                "Temperature",
                "Radiation",
                "time",
                "Humidity")

for(i in 1:length(UL2$Humidity)){
  if(UL2$Humidity[i]<10) UL2$Humidity[i]=UL2$Humidity[i]*10
}

## PAIRWISE
pairs.panels(UL2[, c(1,2,3,5)], 
             cex=0.5,
             method = "spearman", # correlation method
             hist.col = "#00AFBB",
             density = TRUE,  # show density plots
             ellipses = TRUE # show correlation ellipses
)


#############################
## CONTROL CHART - PHASE 2 ##
#############################
UL_f2=read.csv("muhat_UL.csv")[-1209,-1]
colnames(UL_f2)=c("mu","LI","LS")

summary(UL_f2) # CONDITIONAL PARAMETER (mu) estimation

Process_Data <- data.frame(
  Run_Number=c(1:dim(UL_f2)[1]),             #Run Order    
  Y=UL2$Humidity[-1]/100,
  Value = UL_f2 #Process A Random Data
)

Process_Data %>% mutate(Color = ifelse((Y-Value.mu) < (Value.LI-Value.mu) | (Y-Value.mu) > (Value.LS-Value.mu),  "OUT-OF-CONTROL","IN-CONTROL")) %>%
  ggplot() + #init ggplot
  geom_hline(yintercept=0)+
  geom_line(aes(x = Run_Number, y=Value.LI-Value.mu),linetype="twodash",cex=0.25)+
  geom_line(aes(x = Run_Number, y=Value.LS-Value.mu))+
  geom_point(aes(x = Run_Number, y=Y-Value.mu, color = Color),cex=0.5)+
  ylim(c(-0.5,0.5))+xlab("Time Points") + ylab("Control Value")+ 
  theme(text=element_text(size=16),
        title=element_text(size=14))

p0=ggplot(UL2,aes(x=c(1:dim(UL2)[1]),y=Wind)) + #init ggplot
  geom_smooth()+xlab("Time Points")+ylab("WIND SPEED")
p1=ggplot(UL2,aes(x=c(1:dim(UL2)[1]),y=Temperature)) + #init ggplot
  geom_step()+xlab("Time Points")+ylab(expression(paste(Delta,"TEMPERATURE")))
p2=ggplot(UL2,aes(x=c(1:dim(UL2)[1]),y=Radiation)) + #init ggplot
  geom_point()+xlab("Time Points")+ylab("SOLAR RADIATION")

multiplot(p0, p1, p2)
