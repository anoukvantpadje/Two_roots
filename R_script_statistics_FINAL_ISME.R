#R-script to analyze data belonging to manuscript:
#Temporal tracking of quantum-dot apatite across in-vitro mycorrhizal networks shows how host demand can influence fungal nutrient transfer strategies
#published in The ISME Journal
#Author: Anouk van 't Padje (anouk.vantpadje@wur.nl)

#######################################################
#######################################################
#load libraries
library(markdown)
library(lme4)
library(lmerTest)
library(Rmisc)
library(car)
require(MASS)
library(multcompView)
library(emmeans)

#the error bar function
error.bar <- function(x, y, upper, lower=upper, length=0.1, ...)
  {
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, col= "black", ...)
}

#######################################################
#######################################################
#set working directory
setwd("")

#######################################################
#######################################################
# open data
#root and hyphal anbundance
root_hyphal_abundance <-  read.csv("root_hyphal_abundance.csv")
root_hyphal_abundance$replicate <- as.factor(root_hyphal_abundance$replicate)
root_hyphal_abundance$compartment_nr <- as.factor(root_hyphal_abundance$compartment_nr)
summary(root_hyphal_abundance)

#quantum-dot
QD<- read.csv("Quantum_dot_uptake.csv")
QD$replicate <- as.factor(QD$replicate)
QD$compartment_nr <- as.factor(QD$compartment_nr)
summary(QD)

#quantum-dot ratio
ratio_QD<-read.csv("ratio_QD_young_established.csv")
ratio_QD$replicate <- as.factor(ratio_QD$replicate)
ratio_QD$compartment_nr <- as.factor(ratio_QD$compartment_nr)
summary(ratio_QD)

#RICS analysis
RICS <- read.csv("RICS.csv")
RICS$replicate <- as.factor(RICS$replicate)
RICS$color <- as.factor(RICS$color)
summary(RICS)

#P concentration of host roots
p <- read.csv("P_concentration.csv")
summary(p)

# %P from QDS
P_from_QDs <- read.csv("Quantum_dot_uptake_P_concentration.csv")
summary(P_from_QDs)

#visual hyphal colonization
colonization <- read.csv("Intraradical_Colonization.csv")
summary(colonization)

########################################################
########################################################
# Dry root weight
# omit missing values
rootweight_data <- root_hyphal_abundance[complete.cases(root_hyphal_abundance$rootweight),]

# all data together not normal:
shapiro.test(rootweight_data$rootweight.mg.)

# young and established roots sepparately are normal:
shapiro.test(rootweight_data$rootweight.mg.[rootweight_data$compartment=="established"])
shapiro.test(rootweight_data$rootweight.mg.[rootweight_data$compartment=="young"])

test_established <- (lm(rootweight.mg.~ treatment, data=rootweight_data[rootweight_data$compartment=="established",]))
test_young <- (lm(rootweight.mg.~ treatment, data=rootweight_data[rootweight_data$compartment=="young",]))

#linar model is good to study sepparately
leveneTest(test_established)
leveneTest(test_young)

qqnorm(resid(test_established))
qqnorm(resid(test_young))

#So analyse rootweight with linear model
test_all <- lmer(rootweight.mg.~ treatment*compartment+ (1|replicate), data=rootweight_data)
test_all_Levene <- lm(rootweight.mg.~ treatment*compartment, data=rootweight_data)

anova(test_all, ddf="Kenward-Roger", type="II")

#means and SE
means_root <- summarySE(rootweight_data, measurevar = "rootweight.mg.", groupvars = c("treatment", "compartment"))

################################################################################
################################################################################
#Figure 2
#quantum-dot appatite per mg

#remove missing values
QD$color <- as.factor(QD$color)
QD<- QD[complete.cases(QD$QD.nmol.),]

#all data not normally disributed
shapiro.test(QD$QD.nmol.)# NOT normal

#data transformations
QD$QD.mol.<- (1000000000*QD$QD.nmol.)
QD$QD.mol._sqrt <- (QD$QD.mol.)^(1/6)
shapiro.test(QD$QD.mol._sqrt)# NOT normal
hist(QD$QD.mol._sqrt)

#nomality better
shapiro.test(QD$QD.mol._sqrt[QD$color=="1" & QD$compartment=="young"])
shapiro.test(QD$QD.mol._sqrt[QD$color=="1" & QD$compartment=="established"])

shapiro.test(QD$QD.mol._sqrt[QD$color=="2" & QD$compartment=="young"])
shapiro.test(QD$QD.mol._sqrt[QD$color=="2" & QD$compartment=="established"])

shapiro.test(QD$QD.mol._sqrt[QD$color=="3" & QD$compartment=="young"])
shapiro.test(QD$QD.mol._sqrt[QD$color=="3" & QD$compartment=="established"])

#models
test_1_young <- lm(QD.mol._sqrt~ treatment, data=QD[QD$color=="1"& QD$compartment=="young",])
test_2_young <- lm(QD.mol._sqrt~ treatment, data=QD[QD$color=="2"& QD$compartment=="young",])
test_3_young <- lm(QD.mol._sqrt~ treatment, data=QD[QD$color=="3"& QD$compartment=="young",])

test_1_established <- lm(QD.mol._sqrt~ treatment, data=QD[QD$color=="1" & QD$compartment=="established",])
test_2_established <- lm(QD.mol._sqrt~ treatment, data=QD[QD$color=="2" & QD$compartment=="established",])
test_3_established <- lm(QD.mol._sqrt~ treatment, data=QD[QD$color=="3" & QD$compartment=="established",])

#homogeneity better
leveneTest(test_1_young)
leveneTest(test_2_young)
leveneTest(test_3_young)

leveneTest(test_1_established)
leveneTest(test_2_established)
leveneTest(test_3_established)

#normality better
qqnorm(resid(test_1_young))
qqnorm(resid(test_2_young))
qqnorm(resid(test_3_young))

qqnorm(resid(test_1_established))
qqnorm(resid(test_2_established))
qqnorm(resid(test_3_established))


#models
#use the six together in one test
test_all <- lmer(QD.mol._sqrt~ color*treatment*compartment+ (1|replicate), data=QD)
test_all_anova <-anova(test_all, ddf="Kenward-Roger", type="II")

#plot
means <- summarySE(QD, measurevar = "QD.nmol.", groupvars = c("treatment_compartment", "color"), na.rm= TRUE )

pdf("Figure_2.pdf", width = 3.385827, height = 3.385827)
op <- par(mfrow=c(1,1), mar= c(3,3,1,1), xpd=NA)

bar <- barplot(c(means$QD.nmol.[1], means$QD.nmol.[4],means$QD.nmol.[2],means$QD.nmol.[5],means$QD.nmol.[3],means$QD.nmol.[6],
                 means$QD.nmol.[7],means$QD.nmol.[10],means$QD.nmol.[8],means$QD.nmol.[11],means$QD.nmol.[9],means$QD.nmol.[12]),
               col =c("cyan3", "cyan3", "yellow2", "yellow2",  "coral1",  "coral1"),
               density=c(NA, 40), space=c(0,0,0.5,0,0.5,0,2,0,0.5,0,0.5,0), xlab=NA, ylab=NA,
               ylim=c(0,max(means$QD.nmol.+4*means$se)), axes=FALSE, cex.lab=0.7, cex.axis = 0.5, cex.names = 0.7)

axis(side=1, at=c(1,3.5, 6, 10,12.5,15), label=c("21", "14", "7", "21", "14", "7"), cex.axis=0.5, tck=c(0), padj = -4 )
axis(side=1, at=c(-0.5, 8, 16.5), label=c("", "", ""), tck=c(-0.03))
axis(side=1, at=c(2.25, 4.75, 11.25, 13.75), label=c("", "", "", ""), tck=c(-0.01))

axis(side=2, at=c(0.000,0.0002,0.0004,0.0006,0.0008,0.0010,0.0012, 0.0014),
     label=c("0.0", "2.0","4.0", "6.0",  "8.0", "10.0",  "12.0", "14.0"), cex.axis=0.5, las=1, tck=c(-0.02), hadj=0.3)
axis(side=2, at=c( 0.0001,0.0003,0.0005,0.0007,0.0009,0.0011,0.0013), label=c( "",  "",  "",  "",  "",  "",  ""), cex.axis=0.7, tck=c(-0.01))

error.bar(bar, c(means$QD.nmol.[1], means$QD.nmol.[4],means$QD.nmol.[2],means$QD.nmol.[5],means$QD.nmol.[3],means$QD.nmol.[6], 
                 means$QD.nmol.[7],means$QD.nmol.[10],means$QD.nmol.[8],means$QD.nmol.[11],means$QD.nmol.[9],means$QD.nmol.[12]),
          c(means$se[1], means$se[4],means$se[2],means$se[5],means$se[3],means$se[6], 
            means$se[7],means$se[10],means$se[8],means$se[11],means$se[9],means$se[12]), length=0.05)

text(0.5, means$QD.nmol.[1]+means$se[1]+0.0001, "*", cex=1)
text(3,  means$QD.nmol.[2]+means$se[2]+0.0001, "*" , cex=1)
text(5.5,  means$QD.nmol.[3]+means$se[3]+0.0001, "*" , cex=1)

text(9.5,  means$QD.nmol.[7]+means$se[7]+0.0001, "*" , cex=1)
text(12, means$QD.nmol.[8]+means$se[8]+0.0001, "NS", cex=0.5 )
text(14.5,  means$QD.nmol.[9]+means$se[9]+0.0001, "NS" , cex=0.5)

text(2, 0.0014, "Control", cex=0.7)
text(10, 0.0014, "Low-P", cex=0.7)

text(-3, 0.0001, expression(paste('Quantum-dot phosphorous per mg of root (x10 '^ -4*')', sep='')), srt=90, adj=0, cex=0.7)
text(8, -0.00015, "Days since quantum-dot phosphorous injection", cex=0.7)

legend(10, 0.0013, legend=c( "Established roots", "Young roots"), fill="grey70", density = c(NA,40), box.col=NA, bg=NA, cex=0.7)

par(op)
dev.off()

#Means and SE
means <- summarySE(QD, measurevar = "QD.nmol.", groupvars = c("treatment_compartment", "color"), na.rm= TRUE )
means <- summarySE(QD, measurevar = "QD.nmol.root.", groupvars = c("treatment_compartment", "color"), na.rm= TRUE )

# % difference first week
#control
(means$QD.nmol.[1]- means$QD.nmol.[4])/(means$QD.nmol.[4]+means$QD.nmol.[1])
#low P
(means$QD.nmol.[7]- means$QD.nmol.[10])/(means$QD.nmol.[7]+means$QD.nmol.[10])


#######################################################################################
#######################################################################################
#Figure 3
#Ratio of quantum dots in host (young/established)
#ratio non-injected / injected roots
ratio <- ratio_QD[complete.cases(ratio_QD$ratio_young.established),]
ratio$color <- as.factor(ratio$color)

#not normally distributed, right squwed
shapiro.test(ratio$ratio_young.established)
hist(ratio_QD$ratio_young.established)

#transform data
ratio$ratio_young.established_sqrt <- (ratio$ratio_young.established)^(1/16)
hist(ratio$ratio_young.established_sqrt)
shapiro.test(ratio$ratio_young.established_sqrt)

#now normal, so lmer
model_ratios <- lmer(ratio_young.established_sqrt ~ color*treatment + (1|replicate), data=ratio)
leveneTest(model_ratios)
qqnorm(resid(model_ratios))

anova(model_ratios, ddf="Kenward-Roger", type="II")

t.test(ratio$ratio_young.established_sqrt[ratio$treatment=="control"& ratio$color=="1"], ratio$ratio_young.established_sqrt[ratio$treatment=="lowP"& ratio$color=="1"])
t.test(ratio$ratio_young.established_sqrt[ratio$treatment=="control"& ratio$color=="2"],ratio$ratio_young.established_sqrt[ratio$treatment=="lowP"& ratio$color=="2"] )
t.test(ratio$ratio_young.established_sqrt[ratio$treatment=="control"& ratio$color=="3"], ratio$ratio_young.established_sqrt[ratio$treatment=="lowP"& ratio$color=="3"])


#plot
means<- summarySE(data=ratio, measurevar = "ratio_young.established", groupvars = c("treatment", "color"), na.rm= TRUE )

pdf("Figure_3.pdf", width = 3.385827, height = 3.385827)
op <- par(mfrow=c(1,1), mar= c(2.0,2.5,1,1), xpd=NA)

plot <- plot(c(1,2,3,1,2,3), c(means$ratio_young.established),
             ylab = "",  xlab = "", ylim=c(0,3), axes=FALSE, cex.lab=0.7, cex.axis = 0.5)

lines(x=c(1,2,3), c(means$ratio_young.established[1],means$ratio_young.established[2], means$ratio_young.established[3]), lty=1)
lines(x=c(1,2,3), c(means$ratio_young.established[4],means$ratio_young.established[5], means$ratio_young.established[6]), lty=2)

arrows(x0=1, y0=means$ratio_young.established[1]-means$se[1], x1=1, y1=means$ratio_young.established[1]+means$se[1], angle=90, code=3, length=0.05)
arrows(x0=2, y0=means$ratio_young.established[2]-means$se[2], x1=2, y1=means$ratio_young.established[2]+means$se[2], angle=90, code=3, length=0.05)
arrows(x0=3, y0=means$ratio_young.established[3]-means$se[3], x1=3, y1=means$ratio_young.established[3]+means$se[3], angle=90, code=3, length=0.05)
arrows(x0=1, y0=means$ratio_young.established[4]-means$se[4], x1=1, y1=means$ratio_young.established[4]+means$se[4], angle=90, code=3, length=0.05)
arrows(x0=2, y0=means$ratio_young.established[5]-means$se[5], x1=2, y1=means$ratio_young.established[5]+means$se[5], angle=90, code=3, length=0.05)
arrows(x0=3, y0=means$ratio_young.established[6]-means$se[6], x1=3, y1=means$ratio_young.established[6]+means$se[6], angle=90, code=3, length=0.05)

axis(side=1, at=c(1,2,3), label=c( "21", "14", "7"), cex.axis=0.5, tck=c(-0.02), padj = -3)
axis(side=2, at=c(0, 0.5, 1, 1.5, 2.0, 2.5, 3), label=c("0.0", "0.5",  "1",  "1.5",  "2",  "2.5",  "3"),  cex.axis=0.5, las=1, tck=c(-0.02), hadj=0.3)
axis(side=2, at=c(0, 0.25, 0.75, 1.25, 1.75, 2.25, 2.75), label=c("", "", "",  "", "",  "",  ""), cex.axis=0.7, tck=c(-0.01))

text(1, means$ratio_young.established[1]+means$se[1]+0.15, "NS", cex=0.5)
text(2, means$ratio_young.established[5]+means$se[5]+0.15, "NS" , cex=0.5)
text(3, means$ratio_young.established[6]+means$se[6]+0.15, "*" , cex=1)

text(0.6, 1.5, "Ratio of quantum-dot phosphorous \n in young/established roots", srt=90, adj=0.5, cex=0.7)
text(2, -0.45, "Days since quantum-dot phosphorous injection", cex=0.7)

legend("topleft", legend= c("Low-P", "Control"), lty=c(2,1), bty="n", cex=0.7)

par(op)
dev.off()

#Supplementary means and SE
means<- summarySE(data=ratio, measurevar = "ratio_young.established", groupvars = c("treatment", "color"), na.rm= TRUE )


############################################################
############################################################
# Figure 4 
#P content & P from QDS

summary(p)

#P per mg of root
shapiro.test(p$X.ng.P.mg.root.)
barplot(p$X.ng.P.mg.root.)
lm <- lm (X.ng.P.mg.root.~treatment*compartment, data= p)
anova(lm)

boxplot(X.ng.P.mg.root.~treatment_compartment, data= p)
means_p<- summarySE(data= p, measurevar = "conc..nmol.P.mg.root.", groupvars = "treatment_compartment")

# %P from QDS
summary(P_from_QDs)
P_from_QDs <- P_from_QDs[complete.cases(P_from_QDs$X.PfromQD),]

#total p from QDs
means_P_from_QD<- summarySE(data= P_from_QDs, measurevar = "X.PfromQD", groupvars = "treatment")
shapiro.test(P_from_QDs$X.PfromQD)
hist(P_from_QDs$X.PfromQD)
P_from_QDs$X.PfromQD_transformed <- log10(P_from_QDs$X.PfromQD)
shapiro.test(P_from_QDs$X.PfromQD_transformed)
lm <- lm (X.PfromQD_transformed~treatment*compartment, data= P_from_QDs)
anova(lm)

#%P from 1st t QD-apatite addition (21 days)
means_1<- summarySE(data= P_from_QDs, measurevar = "X.Pfrom1st", groupvars = "treatment_compartment")
shapiro.test(P_from_QDs$X.Pfrom1st)
hist(P_from_QDs$X.Pfrom1st)
P_from_QDs$X.Pfrom1st_transformed <- log10(P_from_QDs$X.Pfrom1st)
shapiro.test(P_from_QDs$X.Pfrom1st_transformed)

lm <- lm (X.Pfrom1st_transformed~treatment*compartment, data= P_from_QDs)
anova(lm)

#%P from 2nd t QD-apatite addition (14 days)
means_2<- summarySE(data= P_from_QDs, measurevar = "X.Pfrom2nd", groupvars = "treatment_compartment")
shapiro.test(P_from_QDs$X.Pfrom2nd)
hist(P_from_QDs$X.Pfrom2nd)

P_from_QDs$X.Pfrom2nd_transformed <- log10(P_from_QDs$X.Pfrom2nd)
shapiro.test(P_from_QDs$X.Pfrom2nd_transformed)

lm <- lm (X.Pfrom2nd_transformed~treatment*compartment, data= P_from_QDs)
anova(lm)

#%P from 3de QD-apatite addition (7 days)
means_3<- summarySE(data= P_from_QDs, measurevar = "X.Pfrom3rd", groupvars = "treatment_compartment")
shapiro.test(P_from_QDs$X.Pfrom3rd)
hist(P_from_QDs$X.Pfrom3rd)

P_from_QDs$X.Pfrom3rd_transformed <- log10(P_from_QDs$X.Pfrom3rd)
shapiro.test(P_from_QDs$X.Pfrom3rd_transformed)
hist(P_from_QDs$X.Pfrom3rd_transformed)

lm <- lm(X.Pfrom3rd_transformed~treatment*compartment, data= P_from_QDs)
qqnorm(resid(lm))
Anova(lm, test="F")


#plot
pdf("Figure_4.pdf", width = 6.66, height = 6.66)

op <- par(mfrow=c(2,1),mar=c(4,5,1,1), xpd=NA)

# Figure 4A: %P in host roots
#total %P from QDs
bar <- barplot(c(means_P_from_QD$X.PfromQD), space=c(0,0, 0.5,0), col = c("grey", "grey90"),
               xlab= "Treatment", ylab="P from Quantum-dot appatite", ylim=c(0,20), axes=FALSE)
error.bar(bar, c(means_P_from_QD$X.PfromQD), c(means_P_from_QD$se), length=0.05)

axis(side=1, at=c( 1,3.5), label=c("Control","Low P"), cex.axis=1, tck=0, adj=0.5)
axis(side=1, at=c(-0.3,2.25,4.6), label=c("", "", ""), tck=-0.03)
axis(side=2, at=c(0,5,10,15), label=c("0","5", "10","15"), cex.axis=1, las=1, tck=c(-0.02), hadj=0.3)
axis(side=2, at=c(0, 2.5, 7.5, 12.5), label=c("", "", "",""), cex.axis=1, tck=c(-0.01))

legend("topleft", legend=c( "Established roots", "Young roots"), fill=c("grey", "grey90"), box.col=NA, bg=NA, cex=1)

# %P from 1st to third addition
bar <- barplot(c(means_1$X.Pfrom1st[1], means_1$X.Pfrom1st[2],
                 means_2$X.Pfrom2nd[1],means_2$X.Pfrom2nd[2],
                 means_3$X.Pfrom3rd[1],means_3$X.Pfrom3rd[2],
                 means_1$X.Pfrom1st[3],means_1$X.Pfrom1st[4],
                 means_2$X.Pfrom2nd[3],means_2$X.Pfrom2nd[4],
                 means_3$X.Pfrom3rd[3],means_3$X.Pfrom3rd[4]),
               col =c("cyan3", "cyan1", "yellow3", "yellow1",  "coral3",  "coral1"),
               space=c(0,0,0.5,0,0.5,0,2,0,0.5,0,0.5,0), xlab= "Days since quantum-dot phosphorous injection", ylab="P from Quantum-dot appatite",
               ylim=c(0,15), axes=FALSE, cex.lab=1, cex.axis = 1, cex.names = 1)

axis(side=1, at=c(1,3.5, 6, 10,12.5,15), label=c("21", "14", "7", "21", "14", "7"), cex.axis=1, tck=c(0), padj = -1 )
axis(side=1, at=c(-0.5, 8, 16.5), label=c("", "", ""), tck=c(-0.03))
axis(side=1, at=c(2.25, 4.75, 11.25, 13.75), label=c("", "", "", ""), tck=c(-0.01))

axis(side=2, at=c(0,5,10,15), label=c("0","5", "10","15"), cex.axis=1, las=1, tck=c(-0.02), hadj=0.3)
axis(side=2, at=c(0, 2.5, 7.5, 12.5), label=c("", "", "",""), cex.axis=1, tck=c(-0.01))

error.bar(bar, c(means_1$X.Pfrom1st[1], means_1$X.Pfrom1st[2],
                 means_2$X.Pfrom2nd[1],means_2$X.Pfrom2nd[2],
                 means_3$X.Pfrom3rd[1],means_3$X.Pfrom3rd[2],
                 means_1$X.Pfrom1st[3],means_1$X.Pfrom1st[4],
                 means_2$X.Pfrom2nd[3],means_2$X.Pfrom2nd[4],
                 means_3$X.Pfrom3rd[3],means_3$X.Pfrom3rd[4]),
          c(means_1$se[1], means_1$se[2],
            means_2$se[1],means_2$se[2],
            means_3$se[1],means_3$se[2],
            means_1$se[3],means_1$se[4],
            means_2$se[3],means_2$se[4],
            means_3$se[3],means_3$se[4]), length=0.05)


text(2, 15, "Control", cex=1)
text(10, 15, "Low-P", cex=1)

legend("topleft", legend=c( "Established roots", "Young roots"), fill=c("grey", "grey90"), box.col=NA, bg=NA, cex=1)

close.screen(all.screens = TRUE)
par(op)
dev.off()


######################################################################################
######################################################################################
# Figure 5
#Quantum dot ratio retention
ratio_RICS <- RICS
ratio_RICS$color <- as.factor(ratio_RICS$color)

ratio_RICS$ratio_log <- log10(ratio_RICS$ratio)

shapiro.test(ratio_RICS$ratio_log) # normal

lm <- lm(ratio_log ~treatment*compartment*color, data= ratio_RICS)

leveneTest(lm)
qqnorm(resid(lm))

Anova(lm, test.statistic = "F")


#plot
means <- summarySE(measurevar ="ratio", groupvars = c("treatment", "color", "compartment"), data=ratio_RICS)

pdf("Figure_4.pdf", width = 3.385827, height = 3.385827)
op <- par(mfrow=c(1,1), mar= c(3,3,1,1), xpd=NA)

bar <- barplot(means$ratio[c(1,7, 2,8, 5,11, 6,12)], ylab="", xlab="", ylim=c(0,1.65), axes=FALSE, 
               space=c(0,0,0.5,0,2,0,0.5,0), 
               col=c("cyan4", "cyan3", "cyan4","cyan3", "coral3", "coral", "coral3", "coral"),  density=c(NA, NA))

error.bar(bar, means$ratio[c(1,7, 2,8, 5,11, 6,12)], means$se[c(1,7, 2,8, 5,11, 6,12)], length = 0.03)
axis(1,c(-0.3, 5.5, 11.3), label=c("", "", ""), tck=-0.03)
axis(1,c(2.25, 8.75), label=c("", ""), tck=-0.01)
axis(1,c(1,3.5,7.5,10), label=c("Established", "Young", "Established", "Young"), cex.axis=0.5, tck=0, padj = -3, adj=0.5)
axis(2,c(0,0.2,0.4,0.6,0.8,1.0, 1.2, 1.4, 1.6), 
     label=c("0.0", "0.20", "0.40","0.60","0.80","1.00","1.20","1.40","1.60"), cex.axis=0.5, las=1, tck=c(-0.02), hadj=0.3)
axis(2,c(0.1,0.3,0.5,0.7,0.9,1.1, 1.3, 1.5), label=c("","", "","","","","",""), cex.axis=0.7, tck=c(-0.01))

lines(x=c(-0.1,11.1), y=c(1,1), col="red", lty=2)
text(1.5,1.60, "21 days", cex=0.7)
text(7.5,1.60, "7 days", cex=0.7)
text(-2.2, 0.8, "Ratio of quantum-dot phosphorous per biovolume \n retained in hyphae/roots ", cex=0.7, srt=90, adj=0.5)
text(5.5, -0.19, "Root sytem", cex=0.7, adj=0.5)

legend(8, 1.50, legend=c( "Control", "Low-P"), fill=c("grey50", "grey70"), density = c(NA,NA), box.col=NA, bg=NA, cex=0.7)

par(op)
dev.off()

#supplementary material: means and SE
means <- summarySE(measurevar ="ratio", groupvars = c("treatment", "compartment", "color"), data=ratio_RICS)


############################################################################################
############################################################################################
#Figure 6
# Fungal abundance

#a) EXTRARADICAL ABDUNCANCE
extr_abundance <- root_hyphal_abundance[complete.cases(root_hyphal_abundance$hyphae_weight.mg.),]

#not normal
shapiro.test(extr_abundance$hyphae_weight.mg.)
hist(extr_abundance$hyphae_weight.mg.)

#right skewed, take sqaure root to correct
extr_abundance$hyphae_weight.mg._sqrt <- (extr_abundance$hyphae_weight.mg.)^(1/2)
hist(extr_abundance$hyphae_weight.mg._sqrt)
shapiro.test(extr_abundance$hyphae_weight.mg._sqrt)

#analyse with linear model
model_extra <- lmer(hyphae_weight.mg._sqrt~treatment*compartment +(1|replicate), data=extr_abundance)
leveneTest(model_extra)
qqnorm(resid(model_extra))

anova(model_extra, ddf="Kenward-Roger", type="II")

post <- lsmeans(model_extra, ~ treatment*compartment)
cld(post, alpha=0.05, Letters=letters, adjust="tukey")


#plot
means <- summarySE(extr_abundance, measurevar = "hyphae_weight.mg.", groupvars = c("treatment", "compartment"), na.rm= TRUE )

pdf("Figure_5.pdf", width = 3.385827, height = 3.385827)
op <- par(mfrow=c(1,1), mar= c(2.0,2.5,1,1), xpd=NA)

bar <- barplot(c(means$hyphae_weight.mg.), col=c("black", "white", "black", "white"), space=c(0,0, 0.5,0),
               xlab= "", ylab="", ylim=c(0,max(means$hyphae_weight.mg.+4*means$se)), axes=FALSE)
error.bar(bar, c(means$hyphae_weight.mg.), c(means$se), length=0.05)

axis(side=1, at=c( 1,3.5), label=c("Control","Low P"), cex.axis=0.5, tck=0, padj = -3, adj=0.5)
axis(side=1, at=c(-0.1,2.25,4.6), label=c("", "", ""), tck=-0.03)
axis(side=2, at=c(0,1,2,3,4,5,6,7), label=c("0.0","1.0","2.0","3.0","4.0","5.0","6.0","7.0"), cex.axis=0.5, las=1, tck=c(-0.02), hadj=0.3)
axis(side=2, at=c(0, 0.5, 1.5, 2.5,3.5, 4.5, 5.5, 6.5), label=c("", "", "","", "",  "",  "",  ""), cex.axis=0.7, tck=c(-0.01))

legend("topright", legend=c( "Established root compartment", "Young root compartment"), fill=c("black", "white"), density = c(NA,40), box.col=NA, bg=NA, cex=0.7)
text(-1, 4, "Extra-radical hyphal biomass (mg)", cex=0.7, srt=90, adj=0.5)
text(2.25, -0.75, "Treatment", cex=0.7, adj=0.5)

text(0.5,6.5, "NS", cex=0.5)
text(3,3.75, "NS", cex=0.5)

par(op)
dev.off()

#means and SE for supplementary material
means <- summarySE(extr_abundance, measurevar = "hyphae_weight.mg.", groupvars = c("treatment", "compartment"), na.rm= TRUE )

#B) INTRARADICAL HYPHAL ABUNDANCE OR HYPHAL COLONIZATION
#remove missing data
intra_abundance <- root_hyphal_abundance[complete.cases(root_hyphal_abundance$log_copynr_intraradical_hyphae),]
summary(intra_abundance)

#caculate total copy numbers and take the log again
intra_abundance$total_log_copynr_intraradical <- (log10((10^intra_abundance$log_copynr_intraradical_hyphae) * intra_abundance$rootweight.mg.))

#remove missing data
intra_abundance <- intra_abundance[complete.cases(intra_abundance$total_log_copynr_intraradical),]

#check data > not normal
shapiro.test(intra_abundance$total_log_copynr_intraradical)
hist(intra_abundance$total_log_copynr_intraradical)

#transform. Still not normal, but bimodal
intra_abundance$total_log_copynr_intraradical_log <- (intra_abundance$total_log_copynr_intraradical)^(1/3)
hist(intra_abundance$total_log_copynr_intraradical_log)
shapiro.test(intra_abundance$total_log_copynr_intraradical_log)

#check normality of the two main groups:
shapiro.test(intra_abundance$total_log_copynr_intraradical_log[intra_abundance$compartment=="established"])
shapiro.test(intra_abundance$total_log_copynr_intraradical_log[intra_abundance$compartment=="young"])

#two distributions
model_intra_young <- glm(total_log_copynr_intraradical_log~treatment, data=intra_abundance[intra_abundance$compartment=="young",])
model_intra_establisehd <- glm(total_log_copynr_intraradical_log~treatment, data=intra_abundance[intra_abundance$compartment=="established",])

leveneTest(model_intra_young)
leveneTest(model_intra_establisehd)

qqnorm(resid(model_intra_young))
qqnorm(resid(model_intra_establisehd))

# complete model
model_intra <- lmer(total_log_copynr_intraradical_log~treatment*compartment +(1|replicate), data=intra_abundance)
anova(lmer(total_log_copynr_intraradical_log~treatment*compartment +(1|replicate), data=intra_abundance), ddf="Kenward-Roger", type="II")

post <- lsmeans(model_intra, ~ treatment*compartment)
cld(post, alpha=0.05, Letters=letters, adjust="tukey")


#plot 
means <- summarySE(intra_abundance, measurevar = "log_copynr_intraradical_hyphae", groupvars = c("treatment", "compartment"), na.rm= T )

pdf("Figure_5B.pdf", width = 3.385827, height = 3.385827)
op <- par(mfrow=c(1,1), mar= c(3,3,1,1), xpd=NA)

bar <- barplot(c(means$log_copynr_intraradical_hyphae), space=c(0,0, 0.5,0), col = c("grey", "grey90"),
               xlab= "", ylab="", ylim=c(0,5.2), axes=FALSE)
error.bar(bar, c(means$log_copynr_intraradical_hyphae), c(means$se), length=0.05)

axis(side=1, at=c( 1,3.5), label=c("Control","Low P"), cex.axis=0.5, tck=0, padj = -3, adj=0.5)
axis(side=1, at=c(-0.1,2.25,4.6), label=c("", "", ""), tck=-0.03)
axis(side=2, at=c(0,1,2,3,4,5), label=c("0.0","1.0","2.0","3.0","4.0","5.0"), cex.axis=0.5, las=1, tck=c(-0.02), hadj=0.3)
axis(side=2, at=c(0, 0.5, 1.5, 2.5,3.5, 4.5), label=c("", "", "","", "",""), cex.axis=0.7, tck=c(-0.01))

legend("topright", legend=c( "Established roots", "Young roots"), fill=c("grey", "grey90"), box.col=NA, bg=NA, cex=0.7)
text(-1, 2.5, "Intra-radical colonization (log(copy #))", cex=0.7, srt=90, adj=0.5)
text(2.25, -0.5, "Treatment", cex=0.7, adj=0.5)

par(op)
dev.off()

#descriptive
#caculate absolute copy numbers per mg
intra_abundance$abs_copynr_intraradical <- (10^intra_abundance$log_copynr_intraradical_hyphae)

#control vs low P
means_treatment <- summarySE( data = intra_abundance, measurevar = "abs_copynr_intraradical", groupvars = ("treatment"))
means_treatment[1,3]/means_treatment[2,3]

#young vs old
means_compartment <- summarySE( data = intra_abundance, measurevar = "abs_copynr_intraradical", groupvars = ("compartment"))
means_compartment[2,3]/means_compartment[1,3]

#supplementary material: means and SE
means <- summarySE(intra_abundance, measurevar = "abs_copynr_intraradical", groupvars = c("treatment", "compartment"), na.rm= T )

#############################################################################
##############################################################################
#Figure S1
#P concentration of host roots
hist(p$conc..nmol.P.mg.root.)
shapiro.test(p$conc..nmol.P.mg.root.[p$compartment=="established"])
shapiro.test(p$conc..nmol.P.mg.root.[p$compartment=="young"])

p$conc..nmol.P.mg.root.transformed <- log10(p$conc..nmol.P.mg.root.)
shapiro.test(p$conc..nmol.P.mg.root.transformed)
shapiro.test(p$conc..nmol.P.mg.root.transformed[p$compartment=="established"])
shapiro.test(p$conc..nmol.P.mg.root.transformed[p$compartment=="young"])
hist(p$conc..nmol.P.mg.root.transformed)


lm <- glm(conc..nmol.P.mg.root.~treatment*compartment, family= "gaussian", data= p)
qqnorm(lm)
leveneTest(lm)

anova(lm, test= "F")

#plot
pdf("Figure_S1.pdf", width = 3.33, height = 3.333)
op <- par(mfrow=c(1,1), mar= c(2,5,1,0), xpd=NA)

means_p<- summarySE(data= p, measurevar = "conc..nmol.P.mg.root.", groupvars = "treatment_compartment")
bar <- barplot(c(means_p$conc..nmol.P.mg.root.), space=c(0,0, 0.5,0), col = c("grey", "grey90"),
               xlab= "Treatment", ylab="P concentration (nmol/mg root)", ylim=c(0,30), axes=FALSE)
error.bar(bar, c(means_p$conc..nmol.P.mg.root.), c(means_p$se), length=0.05)

axis(side=1, at=c( 1,3.5), label=c("Control","Low P"), cex.axis=0.7, tck=0, adj=0.5)
axis(side=1, at=c(-0.3,2.25,4.6), label=c("", "", ""), tck=-0.03)
axis(side=2, at=c(0,10,20,30), label=c("0.0","10","20","30"), cex.axis=0.7, las=1, tck=c(-0.02), hadj=0.3)
axis(side=2, at=c(0, 5, 15, 25), label=c("", "", "",""), cex.axis=1, tck=c(-0.01))

legend("topright", legend=c( "Established roots", "Young roots"), fill=c("grey", "grey90"), box.col=NA, bg=NA, cex=0.7)

par(op)
dev.off()

# % difference first week
#control
(means_p$conc..nmol.P.mg.root.[1]- means_p$conc..nmol.P.mg.root.[2])/
  (means_p$conc..nmol.P.mg.root.[1]+means_p$conc..nmol.P.mg.root.[2])

#low P
(means_p$conc..nmol.P.mg.root.[3]- means_p$conc..nmol.P.mg.root.[4])/
  (means_p$conc..nmol.P.mg.root.[3]+means_p$conc..nmol.P.mg.root.[4])


#################################################################
#################################################################
# Figure s3
#% colonization

boxplot(hyphal.colonization..~treatment_compartment, data= colonization)
boxplot(arbuscules..~treatment_compartment, data= colonization)
boxplot(vesicles..~treatment_compartment, data= colonization)

means_colonization <- summarySE(data= colonization, measurevar = "hyphal.colonization..", groupvars = "treatment_compartment")
means_arbuscules <- summarySE(data= colonization, measurevar = "arbuscules..", groupvars = "treatment_compartment")
means_vesicles <- summarySE(data= colonization, measurevar = "vesicles..", groupvars = "treatment_compartment")

shapiro.test(colonization$hyphal.colonization..)
lm <- lm (hyphal.colonization..~treatment*compartment, data= colonization)
anova(lm)

shapiro.test(colonization$arbuscules..)
lm <- lm (arbuscules..~treatment*compartment, data= colonization)
anova(lm)

shapiro.test(colonization$vesicles..)
lm <- lm (vesicles..~treatment*compartment, data= colonization)
anova(lm)

#plot
#a) hyphal colonization
pdf("Figure_S3ABC.pdf", width = 10, height = 3.385827)
op <- par(mfrow=c(1,3), mar= c(5,5,5,5), xpd=NA)

bar <- barplot(c(means_colonization$hyphal.colonization..), space=c(0,0, 0.5,0), col = c("grey", "grey90"),
               xlab= "Treatment", ylab="Hyphal colonization (%)", ylim=c(0,50), axes=FALSE)
error.bar(bar, c(means_colonization$hyphal.colonization..), c(means_colonization$se), length=0.05)

axis(side=1, at=c( 1,3.5), label=c("Control","Low P"), cex.axis=1, tck=0, adj=0.5)
axis(side=1, at=c(-0.1,2.25,4.6), label=c("", "", ""), tck=-0.03)
axis(side=2, at=c(0,10,20,30,40,50), label=c("0.0","10","20","30","40","50"), cex.axis=1, las=1, tck=c(-0.02), hadj=0.3)
axis(side=2, at=c(0, 5, 15, 25,25, 45), label=c("", "", "","", "",""), cex.axis=1, tck=c(-0.01))

legend("topright", legend=c( "Established roots", "Young roots"), fill=c("grey", "grey90"), box.col=NA, bg=NA, cex=1)


#b) Arbuscules
bar <- barplot(c(means_arbuscules$arbuscules..), space=c(0,0, 0.5,0), col = c("grey", "grey90"),
               xlab= "Treatment", ylab="Arbuscules (%)", ylim=c(0,50), axes=FALSE)
error.bar(bar, c(means_arbuscules$arbuscules..), c(means_arbuscules$se), length=0.05)

axis(side=1, at=c( 1,3.5), label=c("Control","Low P"), cex.axis=1, tck=0, adj=0.5)
axis(side=1, at=c(-0.3,2.25,4.6), label=c("", "", ""), tck=-0.03)
axis(side=2, at=c(0,10,20,30,40,50), label=c("0.0","10","20","30","40","50"), cex.axis=1, las=1, tck=c(-0.02), hadj=0.3)
axis(side=2, at=c(0, 5, 15, 25,25, 45), label=c("", "", "","", "",""), cex.axis=1, tck=c(-0.01))

legend("topright", legend=c( "Established roots", "Young roots"), fill=c("grey", "grey90"), box.col=NA, bg=NA, cex=1)

#c) vesicles
bar <- barplot(c(means_vesicles$vesicles..), space=c(0,0, 0.5,0), col = c("grey", "grey90"),
               xlab= "Treatment", ylab="Vesicles (%)", ylim=c(0,50), axes=FALSE)
error.bar(bar, c(means_vesicles$vesicles..), c(means_arbuscules$se), length=0.05)

axis(side=1, at=c( 1,3.5), label=c("Control","Low P"), cex.axis=1, tck=0, adj=0.5)
axis(side=1, at=c(-0.3,2.25,4.6), label=c("", "", ""), tck=-0.03)
axis(side=2, at=c(0,10,20,30,40,50), label=c("0.0","10","20","30","40","50"), cex.axis=1, las=1, tck=c(-0.02), hadj=0.3)
axis(side=2, at=c(0, 5, 15, 25,25, 45), label=c("", "", "","", "",""), cex.axis=1, tck=c(-0.01))

legend("topright", legend=c( "Established roots", "Young roots"), fill=c("grey", "grey90"), box.col=NA, bg=NA, cex=1)

par(op)
dev.off()


