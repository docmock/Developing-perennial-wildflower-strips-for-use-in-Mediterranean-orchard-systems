library(mvabund)
library(dplyr)
library(tidyr)
library(lme4)
library(MASS)
library(multcomp)
library(car)
library(vegan)
library(bibtex)
library(ggplot2)
library(gllvm)
library(boral)
library(emmeans)

# functions
Indicators<-function(ANOVA_OBJ, VARIABLE, INDICATOR_Num){ 
  # where ANOVA_OBJ is an anova object, 
  # VARIABLE (numeric) variable to test from anova object and refers to the row number it appears in the anova table,
  # INDICATOR_Num (numeric) is the number of indicator species to output
  SORTEDANOV = sort((ANOVA_OBJ)$uni.test[VARIABLE,],decreasing=TRUE, index.return=TRUE)
  TEST=(SORTEDANOV$x[1:INDICATOR_Num]/(ANOVA_OBJ)$table[VARIABLE,3])*100
  P = sort((ANOVA_OBJ)$uni.p[VARIABLE,], index.return=TRUE)
  PVAL = P$x
  list(TEST, PVAL)
}

repro_DF<-function(Input_DF){
  veg<-dplyr::select(Input_DF, ends_with("vegetative"))
  bud<-dplyr::select(Input_DF, ends_with("budding"))
  flo<-dplyr::select(Input_DF, ends_with("flowering"))
  see<-dplyr::select(Input_DF, ends_with("seedset"))
  x<-transmute(Input_DF,
               Vegetative = rowSums(veg[,1:ncol(veg)]),
               Budding = rowSums(bud[,1:ncol(bud)]),
               flowering = rowSums(flo[,1:ncol(flo)]),
               Seedset = rowSums(see[,1:ncol(see)]))
}

#######################################################################################################
# sward height and structural heterogeneity:
swardheight<- read.csv(file="Sward_Height_Data.csv", header=T)
str(swardheight)

## Pre-modelling data manipulation
swardheight$Year<-as.factor(swardheight$Year)
swardheight$Site<-as.factor(swardheight$Site)
swardheight$Treatment<-as.factor(swardheight$Treatment)
swardheight$Date<-as.factor(swardheight$Date)
swardheight.17 <- filter(swardheight, Year == "2017")
swardheight.17$Treatment<-gsub("amwt", "ews", (gsub("smwt", "ews", swardheight.17$Treatment)))
swardheight.17$Treatment<-as.factor(swardheight.17$Treatment)
swardheight.18.19<-swardheight %>% filter(Year == "2018" | Year == "2019" )

cov.df.17<-swardheight.17 %>% group_by(Site, Date, Treatment, Alleyway) %>% 
  summarise(Count = n(), 
            mean.height = mean(Height, na.rm = TRUE),
            SD = sd(Height, na.rm = TRUE), 
            SE = sd(Height, na.rm = TRUE)/sqrt(Count)) %>%
  as.data.frame()
cov.df.17$CoV <- (cov.df.17$SD/cov.df.17$mean.height)*100
cov.df.17 <- na.omit(cov.df.17)
cov.df.17 <- unite(cov.df.17, "Alleyway", c(Site, Treatment, Alleyway), sep = "_", remove = F)
cov.df.17$Treatment <- as.factor(cov.df.17$Treatment)
cov.df.17$Alleyway <- as.factor(cov.df.17$Alleyway)
cov.df.17$Site <- as.factor(cov.df.17$Site)
cov.df.17$Treatment <- relevel(cov.df.17$Treatment, ref = "control")
cov.df.17$Treatment<-droplevels(cov.df.17$Treatment)

cov.df.18.19<-swardheight.18.19 %>% group_by(Year, Date, Site, Treatment, Alleyway) %>% 
  summarize(Count = n(), 
         mean.height = mean(Height, na.rm = TRUE),
         SD = sd(Height, na.rm = TRUE), 
         SE = sd(Height, na.rm = TRUE)/sqrt(Count)) %>%
  as.data.frame()
cov.df.18.19$CoV <- (cov.df.18.19$SD/cov.df.18.19$mean.height)*100
cov.df.18.19 <- na.omit(cov.df.18.19)
cov.df.18.19 <- unite(cov.df.18.19, "Alleyway", c(Site, Treatment, Alleyway), sep = "_", remove = F)
cov.df.18.19$Treatment <- as.factor(cov.df.18.19$Treatment)
cov.df.18.19$Alleyway <- as.factor(cov.df.18.19$Alleyway)
cov.df.18.19$Site <- as.factor(cov.df.18.19$Site)
cov.df.18.19$Treatment <- relevel(cov.df.18.19$Treatment, ref = "control")
str(cov.df.18.19)
cov.df.18.19$Year<-droplevels(cov.df.18.19$Year)
cov.df.18.19$Treatment <- relevel(cov.df.18.19$Treatment, ref = "control")

### 2017
#### Eyeball data:
str(cov.df.17)
par(mfrow=c(2,2))
hist(cov.df.17$CoV) # transform?
cov.df.17$ln_Cov<-log(cov.df.17$CoV+1)
hist(cov.df.17$ln_Cov) # GREAT!
boxplot(cov.df.17$ln_Cov~cov.df.17$Treatment) 
boxplot(cov.df.17$ln_Cov~cov.df.17$Site) 

hist(cov.df.17$mean.height) # transform?
cov.df.17$ln_mean.height<-log(cov.df.17$mean.height+1)
hist(cov.df.17$ln_mean.height)
boxplot(cov.df.17$ln_mean.height~cov.df.17$Treatment) 
boxplot(cov.df.17$ln_mean.height~cov.df.17$Site) 

#### calculate proportional difference in height
hght.sum.17<-cov.df.17 %>% group_by(Treatment) %>% 
  summarize(Count = n(), 
            Mean = mean(mean.height, na.rm = TRUE),
            SD = sd(mean.height, na.rm = TRUE), 
            SE = sd(mean.height, na.rm = TRUE)/sqrt(Count)) %>%
  as.data.frame()
(hght.sum.17[2,3]-hght.sum.17[1,3])*(100/hght.sum.17[1,3])

hght.sum.18.19<-cov.df.18.19 %>% group_by(Treatment) %>% 
  summarize(Count = n(), 
            Mean = mean(mean.height, na.rm = TRUE),
            SD = sd(mean.height, na.rm = TRUE), 
            SE = sd(mean.height, na.rm = TRUE)/sqrt(Count)) %>%
  as.data.frame()
(hght.sum.18.19[2,3]-hght.sum.18.19[1,3])*(100/hght.sum.18.19[1,3])
(hght.sum.18.19[3,3]-hght.sum.18.19[1,3])*(100/hght.sum.18.19[1,3])

#### calculate proportional difference in COV
cov.sum.17<-cov.df.17 %>% group_by(Treatment) %>% 
  summarize(Count = n(), 
            Mean = mean(CoV, na.rm = TRUE),
            SD = sd(CoV, na.rm = TRUE), 
            SE = sd(CoV, na.rm = TRUE)/sqrt(Count)) %>%
  as.data.frame()
(cov.sum.17[2,3]-cov.sum.17[1,3])*(100/cov.sum.17[1,3])

cov.sum.18.19<-cov.df.18.19 %>% group_by(Treatment) %>% 
  summarize(Count = n(), 
            Mean = mean(CoV, na.rm = TRUE),
            SD = sd(CoV, na.rm = TRUE), 
            SE = sd(CoV, na.rm = TRUE)/sqrt(Count)) %>%
  as.data.frame()
(cov.sum.18.19[2,3]-cov.sum.18.19[1,3])*(100/cov.sum.18.19[1,3])
(cov.sum.18.19[3,3]-cov.sum.18.19[1,3])*(100/cov.sum.18.19[1,3])

## Model creation
### 2017 Sward Height
height.17 <- na.omit(cov.df.17)
height17.max <- lmer(ln_mean.height ~ Treatment
                           + (1|Site)+ (1|Date),  data=cov.df.17)
isSingular(height17.max)
summary(height17.max)
Anova(height17.max)
AIC(height17.max)
prwshght17<-emmeans(height17.max, "Treatment")
pairs(prwshght17)

### 2017 structural heterogeneity
CoV17.max <- lmer(ln_Cov ~ Treatment
                  + (1|Site) + (1|Date),  data=cov.df.17)
isSingular(CoV17.max)
summary(CoV17.max)
Anova(CoV17.max)
AIC(CoV17.max)
plot(CoV17.max)
prwscov17<-emmeans(CoV17.max, "Treatment")
pairs(prwscov17)

### 2018 and 2019
par(mfrow=c(2,2))
hist(cov.df.18.19$mean.height) 
cov.df.18.19$ln_Mean<-log(cov.df.18.19$mean.height+1)
hist(cov.df.18.19$ln_Mean)
boxplot(cov.df.18.19$ln_Mean~cov.df.18.19$Treatment)
boxplot(cov.df.18.19$ln_Mean~cov.df.18.19$Site) 

### 2018/19 Sward Height
cov.df.18.19<- na.omit(cov.df.18.19)
Height18.19.max <- lmer(ln_Mean ~ Treatment * Year 
                  + (1|Site) + (1|Date),  data=cov.df.18.19)
isSingular(Height18.19.max)
summary(Height18.19.max)
Anova(Height18.19.max)
#### reduce 
Height18.19.alt <- lmer(ln_Mean ~ Treatment + Year 
                        + (1|Site) + (1|Date),  data=cov.df.18.19)
isSingular(Height18.19.alt)
summary(Height18.19.alt)
Anova(Height18.19.alt)
#### reduce
Height18.19.alt1 <- lmer(ln_Mean ~ Treatment 
                        + (1|Site) + (1|Date),  data=cov.df.18.19)
isSingular(Height18.19.alt1)
summary(Height18.19.alt1)
Anova(Height18.19.alt1)
plot(Height18.19.alt1)
AIC(Height18.19.alt1)
prwshght1819<-emmeans(Height18.19.alt1, "Treatment")
pairs(prwshght1819)

### 2018/19 sward height
hist(cov.df.18.19$CoV)
cov.df.18.19$ln_Cov<-log(cov.df.18.19$CoV+1)
hist(cov.df.18.19$ln_Cov) 
boxplot(cov.df.18.19$ln_Cov~cov.df.18.19$Treatment) 
boxplot(cov.df.18.19$ln_Cov~cov.df.18.19$Site) 
CoV18.19.max <- lmer(ln_Cov ~ Treatment*Year
                               + (1|Site) + (1|Date),  data=cov.df.18.19)
isSingular(CoV18.19.max)
summary(CoV18.19.max)
Anova(CoV18.19.max)
AIC(CoV18.19.max) 
plot(CoV18.19.max)

CoV18.19.alt <- lmer(ln_Cov ~ Treatment + Year 
                     + (1|Site) + (1|Date),  data=cov.df.18.19)
isSingular(CoV18.19.alt)
summary(CoV18.19.alt)
Anova(CoV18.19.alt)
plot(CoV18.19.alt)
AIC(CoV18.19.alt)

CoV18.19.alt1 <- lmer(ln_Cov ~ Treatment 
                     + (1|Site) + (1|Date),  data=cov.df.18.19)
isSingular(CoV18.19.alt1)
summary(CoV18.19.alt1)
Anova(CoV18.19.alt1)
plot(CoV18.19.alt1)
AIC(CoV18.19.alt1) 
prwscov1819<-emmeans(CoV18.19.alt1, "Treatment")
pairs(prwscov1819)

####################################################################################################################################

#  1)	To determine whether the sowing of a habitat composed of 12 forbs species and two grasses can increase resource richness

spp_abundance_data<- read.csv(file="Botanical_Surveys_Data.csv", header =T)
spp.abundance.data<-as.data.frame(spp_abundance_data)
str(spp.abundance.data)

## Pre-modelling data manipulation
spp.abundance.data$Treatment<- factor(spp.abundance.data$Treatment, levels = c("control", "smwt", "amwt"))
spp.abundance.data$Year<-as.factor(spp.abundance.data$Year)
spp.abundance.no.nas<-spp.abundance.data %>% drop_na() # drop rows which total 0
str(spp.abundance.no.nas)

### remove all columns which sum to 0
spp.abundance.no.nas<-spp.abundance.no.nas[, colSums(spp.abundance.no.nas != 0) > 0]
colSums(spp.abundance.no.nas[, 8:307]) # check

### Create environmental datasheet
environmental.data<-spp.abundance.no.nas[, 1:5]
# select only percentage cover data
percent.cover.no.nas<-dplyr::select(spp.abundance.no.nas, ends_with("_PC"))
colSums(percent.cover.no.nas) # check

### Create variable richness
require(vegan)
percent.cover.no.nas$richness<-specnumber(percent.cover.no.nas)

## summarise and eyeball data:
### Summarise by treatment
group_by(percent.cover.no.nas, environmental.data$Treatment) %>%
  summarise(count = n(),
            mean = mean(richness, na.rm = TRUE),
            sd = sd(richness, na.rm = TRUE),
            se = sd(richness, na.rm = TRUE)/sqrt(count),
            median = median(richness, na.rm = TRUE),
            IQR = IQR(richness, na.rm = TRUE))
###  Summarise May by year
group_by(percent.cover.no.nas, environmental.data$Year) %>%
  summarise(count = n(),
            mean = mean(richness, na.rm = TRUE),
            sd = sd(richness, na.rm = TRUE),
            se = sd(richness, na.rm = TRUE)/sqrt(count),
            median = median(richness, na.rm = TRUE),
            IQR = IQR(richness, na.rm = TRUE))
###  Summarise May by treatment*Year 
group_by(percent.cover.no.nas, environmental.data$Treatment, environmental.data$Year) %>%
  summarise(count = n(),
            mean = mean(richness, na.rm = TRUE),
            sd = sd(richness, na.rm = TRUE),
            se = sd(richness, na.rm = TRUE)/sqrt(count),
            median = median(richness, na.rm = TRUE),
            IQR = IQR(richness, na.rm = TRUE))
###  Summarise May by replicate block
group_by(percent.cover.no.nas, environmental.data$Site) %>%
  summarise(count = n(),
            mean = mean(richness, na.rm = TRUE),
            sd = sd(richness, na.rm = TRUE),
            se = sd(richness, na.rm = TRUE)/sqrt(count),
            median = median(richness, na.rm = TRUE),
            IQR = IQR(richness, na.rm = TRUE))

### summarise by species 
summary.spp.abund.tot <- group_by(percent.cover.no.nas, environmental.data$Treatment) %>%
  summarise(across(1:83, mean)) %>%
  as.data.frame()



### Eyeball data:
par(mfrow=c(2,2))
hist(percent.cover.no.nas$richness) # zero inflated
boxplot(percent.cover.no.nas$richness~environmental.data$Treatment)
boxplot(percent.cover.no.nas$richness~environmental.data$Site)
boxplot(percent.cover.no.nas$richness~environmental.data$Year)

## Model creation
### Guassian mod:
richness.mod.pois<-lmer(percent.cover.no.nas$richness~
                         Treatment*Year
                       +(1|Year/Site), data=environmental.data)
summary(richness.mod.pois)
Anova(richness.mod.pois)
deviance(richness.mod.pois)/df.residual(richness.mod.pois) # overdispersed

### Negative Binomial Max mod
nb.richness.max<-glmer.nb(percent.cover.no.nas$richness~
                            Treatment*Year+(1|Site), 
                          data=environmental.data)
summary(nb.richness.max)
plot(nb.richness.max)
Anova(nb.richness.max) 
lsmeans(nb.richness.max, pairwise ~ Treatment)
lsmeans(nb.richness.max, pairwise ~ Year)
lsmeans(nb.richness.max, pairwise ~ Treatment | Year)
lsmeans(nb.richness.max, pairwise ~ Year | Treatment)

########################################################################################################################################
#  2)	To determine how alleyway management affects plant community composition  

str(spp_abundance_data)
## Pre-modelling data manipulation
### select only 2018 and 2019 as requires balance
percent.cover.data.df<-dplyr::filter(spp_abundance_data, Year != "2017")

### Create environmental data sheet
environmental.cover.data<-percent.cover.data.df[, 1:5]

### select only percentage cover data
percent.cover.data<-dplyr::select(percent.cover.data.df, ends_with("_PC"))

### Remove cols which total 0
percent.cover.data<-percent.cover.data[, colSums(percent.cover.data != 0, na.rm=T) > 0]
colSums(percent.cover.data, na.rm=T) # check

environmental.cover.data$Site<-as.factor(environmental.cover.data$Site)
environmental.cover.data$Obs_unit<-as.factor(environmental.cover.data$Obs_unit)
environmental.cover.data<-tidyr::unite(environmental.cover.data, "Site_id", c(Year, Site, Treatment), sep="_", remove=F)
environmental.cover.data$Site_id<-as.factor(environmental.cover.data$Site_id)

# group by site id
spp.cover.grouped <- percent.cover.data%>%
  group_by(environmental.cover.data$Site_id) %>%
  summarise(across(everything(), sum)) %>%
  as.data.frame()
spp.cover.grouped<-as.data.frame(spp.cover.grouped)
environmental.grouped<-as.data.frame(cbind(site_id=levels(environmental.cover.data$site_id), 
                             Year=c( rep("2018", 12), rep("2019", 12)), 
                             Site=rep(c(rep("CAL", 3), rep("MDA", 3), rep("MTP1", 3), rep("MTP2", 3)), 2), 
                             Treatment=rep(c( "amwt", "smwt", "control"), 8)))
environmental.grouped$Treatment<-environmental.grouped$Treatment %>%
  as.factor() %>%
  relevel(ref = "control")
spp.cover.long<-gather(cbind(environmental.grouped, spp.cover.grouped), species, abundance, 6:75)
spp.cover.summary<-group_by(spp.cover.long, species, Treatment) %>%
  summarise(count = n(),
            mean = mean(abundance, na.rm = TRUE),
            sd = sd(abundance, na.rm = TRUE),
            se = sd(abundance, na.rm = TRUE)/sqrt(count),
            median = median(abundance, na.rm = TRUE),
            IQR = IQR(abundance, na.rm = TRUE))
write.table(spp.cover.summary, file="spp_cover_summary.csv", sep=",")

## Model creation
### Create matrix for ordination
spp.cover.matrix=as.matrix(spp.cover.grouped[, 2:72])

### Fit compositional model via manyglm
mv.spp.cover.grouped=mvabund(spp.cover.grouped[, 2:72])
par(mar = c(2, 10, 2, 2)) 
boxplot(spp.cover.grouped[, 2:72], horizontal = TRUE, las = 2, main = "Abundance")
meanvar.plot(mv.spp.cover.grouped)
plot(mv.spp.cover.grouped~as.factor(environmental.grouped$Treatment), cex.axis = 0.8, cex = 0.8)

environmental.grouped<-unite(environmental.grouped, "year.x.farm", c(Year, Site), sep="_", remove=F)
CTRL = how(block = environmental.grouped$year.x.farm)
PermIDs = shuffleSet(length(environmental.grouped$year.x.farm),nset=999,control=CTRL)
head(PermIDs)

par(mfrow=c(1,1))
par(mar=c(2, 2, 2, 0))
par(xpd=F)

ft.1=manyglm(mv.spp.cover.grouped~Treatment*Year, data=environmental.grouped, family="negative.binomial")
plot(ft.1, which = 1:2) # Not a great fit

ft.2=manyglm(mv.spp.cover.grouped~Treatment*Year+Site, data=environmental.grouped, family="negative.binomial")
plot(ft.2, which = 1:2) # plots look much better
anv.ft2<-anova(ft.2, bootID=PermIDs, show.time="all")

ft.3=manyglm(mv.spp.cover.grouped~Treatment+Site+Year, data=environmental.grouped, family="negative.binomial")
plot(ft.3, which=1:2) # this perhaps looks the best
anv.ft3<-anova(ft.3, bootID=PermIDs, show.time="all")
summary(ft.3)

### Is treatment significant
null=manyglm(mv.spp.cover.grouped~Site+Year, data=environmental.grouped, family="negative.binomial")
anv.ft3.uni<-anova(null, ft.3, bootID=PermIDs, p.uni="adjusted", resamp="monte-carlo", show.time="all")
anv.ft3.indicators<-Indicators(anv.ft3.uni, 2, 20)
sort((anv.ft3.uni)$uni.test[2,],decreasing=TRUE, index.return=TRUE)
write.table(anv.ft3.indicators[1], file="Plant_Composition_Test_Stats.csv", sep=",")
write.table(anv.ft3.indicators[2], file="Plant_Composition_Pvalues.csv", sep=",")

### Is treatment*Year significant
anv.ft4.uni<-anova(ft.3, ft.2, bootID=PermIDs,  p.uni="adjusted", show.time="all")

composition.summary.T<-group_by(percent.cover.data, environmental.cover.data$Treatment) %>%
  summarise(count = n(),
            P_lanceolata.mean = mean(SF_P_lanceolata_PC, na.rm = TRUE),
            P_lanceolata.sd = sd(SF_P_lanceolata_PC, na.rm = TRUE),
            P_lanceolata.se = sd(SF_P_lanceolata_PC, na.rm = TRUE)/sqrt(count),
            D_glomerata.mean = mean(SG_D_glomerata_PC, na.rm = TRUE),
            D_glomerata.sd = sd(SG_D_glomerata_PC, na.rm = TRUE),
            D_glomerata.se = sd(SG_D_glomerata_PC, na.rm = TRUE)/sqrt(count),
            A_millefolium.mean = mean(SF_A_millefolium_PC, na.rm = TRUE),
            A_millefolium.sd = sd(SF_A_millefolium_PC, na.rm = TRUE),
            A_millefolium.se = sd(SF_A_millefolium_PC, na.rm = TRUE)/sqrt(count),
            S_verbenaca.mean = mean(SF_S_verbenaca_PC, na.rm = TRUE),
            S_verbenaca.sd = sd(SF_S_verbenaca_PC, na.rm = TRUE),
            S_verbenaca.se = sd(SF_S_verbenaca_PC, na.rm = TRUE)/sqrt(count),
            C_capillaris.mean = mean(UF_C_capillaris_PC, na.rm = TRUE),
            C_capillaris.sd = sd(UF_C_capillaris_PC, na.rm = TRUE),
            C_capillaris.se = sd(UF_C_capillaris_PC, na.rm = TRUE)/sqrt(count),
            T_campestre.mean = mean(UF_T_campestre_PC, na.rm = TRUE),
            T_campestre.sd = sd(UF_T_campestre_PC, na.rm = TRUE),
            T_campestre.se = sd(UF_T_campestre_PC, na.rm = TRUE)/sqrt(count))
options(tibble.print_max = Inf)
write.table(composition.summary.T, file="composition.summary.csv", sep=",")

composition.summary.all<-group_by(percent.cover.data, environmental.cover.data$Treatment) %>%
  summarise(across(everything(), mean))  %>%
  as.data.frame()
options(tibble.print_max = Inf)
write.table(composition.summary.T, file="All_species_summary.csv", sep=",")

## visualise using Boral
Tot_Percent_Matrix=as.matrix(spp.cover.grouped[, 2:72])
rowIDs<-cbind(1:24, rep(1:8,each=3))
Tot_Percent.ord.1 = boral(Tot_Percent_Matrix, lv.control = list(num.lv = 2), family="negative.binomial", 
                          row.eff = "random",
                          row.ids = rowIDs, 
                          save.model=TRUE)
par(mfrow=c(2,2))
plot(Tot_Percent.ord.1, est = "median", include.ranef = TRUE, jitter = FALSE)

####################################################################################################################################
#  3) How does management strategies applied to sown habitats affect the performance of sown species?

str(spp_abundance_data)
## Pre-modelling data manipulation
### Create environmental data sheet
environmental.species.data<-spp_abundance_data[, 1:5]

### select only percentage cover data
species.cover.data<-dplyr::select(spp_abundance_data, ends_with("_PC"))

environmental.species.data$Site<-as.factor(environmental.species.data$Site)
environmental.species.data$Year<-as.factor(environmental.species.data$Year)
environmental.species.data$Treatment<-as.character(environmental.species.data$Treatment)
environmental.species.data$Obs_unit<-as.factor(environmental.species.data$Obs_unit)
environmental.species.data<-unite(environmental.species.data, "SITE_ID", c(Year, Site, Treatment), sep="_", remove=F)
environmental.species.data$SITE_ID<-as.factor(environmental.species.data$SITE_ID)

### filter To only include sown species
species.cover.data<-dplyr::select(species.cover.data, starts_with(c("SG_", "SF_")))
sown.species<-cbind(environmental.species.data, species.cover.data)
### filter To only include sown treatments: Passive management and active management
sown.species<-filter(sown.species, Treatment!="control")
sown.species<-dplyr::rename(sown.species, 
                               A_azurea = SF_A_azurea_PC,
                               A_millefolium = SF_A_millefolium_PC,
                               C_intybus = SF_C_intybus_PC,
                               D_glomerata = SG_D_glomerata_PC,
                               H_perforatum = SF_H_perforatum_PC,
                               M_suaveolens = SF_M_suaveolens_PC,
                               M_vulgare = SF_M_vulgare_PC, 
                               O_natrix = SF_O_natrix_PC,
                               P_bituminosa = SF_P_bituminosa_PC, 
                               P_lanceolata =  SF_P_lanceolata_PC, 
                               S_arundinaceaus = SG_S_arundinaceaus_PC,
                               S_verbenaca = SF_S_verbenaca_PC)
sown.species.ordered<-sown.species[,c("A_azurea", "A_millefolium", "C_intybus", "D_glomerata", 
                       "H_perforatum", "M_suaveolens", "M_vulgare", "O_natrix", "P_bituminosa", "P_lanceolata",
                        "S_arundinaceaus", "S_verbenaca")]
sown.species<-cbind(sown.species[,1:6], sown.species.ordered)                  
sown.species$Treatment<-as.factor(sown.species$Treatment)
sown.species$Treatment<-relevel(sown.species$Treatment, ref = "smwt")

## Model creation
mv.sown.species<-mvabund(sown.species[,7:18])
par(xpd=TRUE)
par(mfrow=c(1,1))
par(mar=c(2,10,1,1))
boxplot(sown.species[,7:18], horizontal = T,las=2, main="mean Abundance")
par(mar=c(2,2,2,2)) 
meanvar.plot(mv.sown.species~sown.species$Treatment)
abline(a=0,b=1,col="green")
plot(mv.sown.species~as.factor(sown.species$Treatment), cex.axis=0.8, cex=0.8)

sown_CTRL = how(block = sown.species$Site)
Sown_permIDs = shuffleSet(length(sown.species$Site),nset=999,control=sown_CTRL)
head(Sown_permIDs)

sown.species.max<-manyglm(mv.sown.species~Treatment*Year,
                                data=sown.species)
plot(sown.species.max, which=1)
plot(sown.species.max, which=2)
anv.sown.max<-anova(sown.species.max,  
                     resamp="monte-carlo",
                     p.uni="adjusted",
                     show.time = "all")
anv.sown.max
sum.sown.max<-summary(sown.species.max,  
                    resamp="monte-carlo",
                    p.uni="adjusted",
                    show.time = "all")
sort((anv.sown.max)$uni.test[4,],decreasing=TRUE, index.return=TRUE)
inicators.sown.treatment<-Indicators(anv.sown.max, 2, 10)
inicators.sown.year<-Indicators(anv.sown.max, 3, 10)
inicators.sown.treatmentXyear<-Indicators(anv.sown.max, 4, 12)
write.table(inicators.sown.treatment[1], file="SownSpecies_Indicators_test_stat_Treatment.csv", sep=",")
write.table(inicators.sown.treatment[2], file="SownSpecies_Indicators_P_values_Treatment.csv", sep=",")
write.table(inicators.sown.year[1], file="SownSpecies_Indicators_test_stat_Year.csv", sep=",")
write.table(inicators.sown.year[2], file="SownSpecies_Indicators_P_values_Year.csv", sep=",")
write.table(inicators.sown.treatmentXyear[1], file="SownSpecies_Indicators_test_stat_TreatmentxYear.csv", sep=",")
write.table(inicators.sown.treatmentXyear[2], file="SownSpecies_Indicators_P_values_TreatmentxYear.csv", sep=",")

### summarise
sown.species.long<-gather(sown.species, species, abundance, 7:18)
sown.species.summary<-group_by(sown.species.long, species, Treatment) %>%
  summarise(count = n(),
            mean = mean(abundance, na.rm = TRUE),
            sd = sd(abundance, na.rm = TRUE),
            se = sd(abundance, na.rm = TRUE)/sqrt(count),
            median = median(abundance, na.rm = TRUE),
            IQR = IQR(abundance, na.rm = TRUE))
write.table(sown.species.summary, file="sown_species_summary.csv", sep=",")

sown.species.Yearsummary<-group_by(sown.species.long, species, Year) %>%
  summarise(count = n(),
            mean = mean(abundance, na.rm = TRUE),
            sd = sd(abundance, na.rm = TRUE),
            se = sd(abundance, na.rm = TRUE)/sqrt(count),
            median = median(abundance, na.rm = TRUE),
            IQR = IQR(abundance, na.rm = TRUE))
write.table(sown.species.Yearsummary, file="sown_species_Yearsummary.csv", sep=",")

sown.species.TxYsummary<-group_by(sown.species.long, species, Year, Treatment) %>%
  summarise(count = n(),
            mean = mean(abundance, na.rm = TRUE),
            sd = sd(abundance, na.rm = TRUE),
            se = sd(abundance, na.rm = TRUE)/sqrt(count),
            median = median(abundance, na.rm = TRUE),
            IQR = IQR(abundance, na.rm = TRUE))
write.table(sown.species.TxYsummary, file="sown_species_TxYsummary.csv", sep=",")



######################################### Differences is resource type ####################################################################
# Classified as:
# Bare soil, Leaf litter, 
# vegetative grasses, budding  grasses, flowering  grasses, seed-set grasses, 
# vegetative forbs, budding  forbs, flowering  forbs, seed-set forbs, 

str(spp_abundance_data)

## Pre-modelling data manipulation
### select only those which are in the seed set stage to identify which are contributing to the control
seedset.data<-dplyr::select(spp.abundance.no.nas, ends_with("_seedset"))
environmental.seed<-spp.abundance.no.nas[, 1:5]
seedset.long<-gather(cbind(environmental.seed, seedset.data), species, abundance, 6:45)
seedset.summary<-group_by(seedset.long, species, Treatment) %>%
  summarise(count = n(),
            mean = mean(abundance, na.rm = TRUE),
            sd = sd(abundance, na.rm = TRUE),
            se = sd(abundance, na.rm = TRUE)/sqrt(count),
            median = median(abundance, na.rm = TRUE),
            IQR = IQR(abundance, na.rm = TRUE))
write.table(seedset.summary, file="seedset_summary.csv", sep=",")

### Create environmental data sheet
environmental.resource.data<-spp_abundance_data[, 1:5]

environmental.resource.data$Site<-as.factor(environmental.resource.data$Site)
environmental.resource.data$Year<-as.factor(environmental.resource.data$Year)
environmental.resource.data$Obs_unit<-as.factor(environmental.resource.data$Obs_unit)
environmental.resource.data<-unite(environmental.resource.data, "SITE_ID", c(Year, Site, Treatment), sep="_", remove=F)
environmental.resource.data$SITE_ID<-as.factor(environmental.resource.data$SITE_ID)

## Total_cover 
### filter by unsown grasses, unsown forbs, sown grasses and sown forbs
unsown.forbs<-dplyr::select(spp_abundance_data, starts_with("UF_"))
unsown.grasses<-dplyr::select(spp_abundance_data, starts_with("UG_"))
sown.forbs<-dplyr::select(spp_abundance_data, starts_with("SF_"))
sown.grasses<-dplyr::select(spp_abundance_data, starts_with("SG_"))

### combine all grasses and all forbs
grasses<-cbind(unsown.grasses, sown.grasses)
forbs<-cbind(unsown.forbs, sown.forbs)

### create new dataframe
forbs.df<-repro_DF(forbs)
grasses.df<-repro_DF(grasses)
forbs.df<-dplyr::rename(forbs.df, 
                        c_Vegetative_Forbs = Vegetative,
                        e_Budding_Forbs = Budding,
                        g_Flowering_Forbs = flowering,
                        i_Seeding_Forbs = Seedset)
grasses.df<-dplyr::rename(grasses.df, 
                          d_Vegetative_Grasses = Vegetative,
                          f_Budding_Grasses = Budding,
                          h_Flowering_Grasses = flowering,
                          j_Seeding_Grasses = Seedset)

a_Bare_soil<-(spp_abundance_data$Soil)
b_Leaf_litter<-(spp_abundance_data$Leaf)
resource.df<-cbind(spp_abundance_data[, 1:4],  a_Bare_soil,  b_Leaf_litter, forbs.df, grasses.df)
resource.df<-mutate(resource.df, Total_cover =
                        c_Vegetative_Forbs + e_Budding_Forbs + g_Flowering_Forbs + i_Seeding_Forbs + 
                        d_Vegetative_Grasses + f_Budding_Grasses + h_Flowering_Grasses + j_Seeding_Grasses +
                  a_Bare_soil + b_Leaf_litter)
resource.df<-mutate(resource.df, 
                      a_Bare_soil_prop = (a_Bare_soil/Total_cover)*100, 
                      b_Leaf_litter_prop = (a_Bare_soil/Total_cover)*100, 
                      c_Vegetative_Forbs_prop = (c_Vegetative_Forbs/Total_cover)*100, 
                      d_Vegetative_Grasses_prop = (d_Vegetative_Grasses/Total_cover)*100, 
                      e_Budding_Forbs_prop = (e_Budding_Forbs/Total_cover)*100, 
                      f_Budding_Grasses_prop = (f_Budding_Grasses/Total_cover)*100, 
                      g_Flowering_Forbs_prop = (g_Flowering_Forbs/Total_cover)*100, 
                      h_Flowering_Grasses_prop = (h_Flowering_Grasses/Total_cover)*100, 
                      i_Seeding_Forbs_prop = (i_Seeding_Forbs/Total_cover)*100, 
                      j_Seeding_Grasses_prop = (j_Seeding_Grasses/Total_cover)*100)
str(resource.df)
resource.df$Treatment<-as.factor(resource.df$Treatment)
resource.df$Year<-as.factor(resource.df$Year)

resource.df$Treatment<-relevel(resource.df$Treatment, ref = "control")

options(pillar.sigfig=5)
group_by(resource.df, Treatment) %>%
  summarise(count = n(),
            mean = mean(a_Bare_soil, na.rm = TRUE),
            se = sd(a_Bare_soil, na.rm = TRUE)/sqrt(count),
            median = median(a_Bare_soil, na.rm = TRUE))
### reshape format from wide to long
resource.long <- gather(resource.df, cover_type, cover, 5:14)
resource.long<-mutate(resource.long, Other_Cover=Total_cover-cover, 
                       covertype_prop=cover/Total_cover)
str(resource.long)
resource.long$cover_type<-as.factor(resource.long$cover_type)
resource.long$LOGIT_Prop_Cover<-logit(resource.long$covertype_prop)
resource.long$LOGIT_Cover<-logit(resource.long$cover)
### summarise
cover_summary<-resource.long  %>% group_by(cover_type, Treatment, .add = TRUE) %>%
  summarise(count = n(),
            mean = mean(cover, na.rm = TRUE),
            se = sd(cover, na.rm = TRUE)/sqrt(count),
            median = median(cover, na.rm = TRUE))
write.table(cover_summary, "resource_cover_summary_all.csv", sep=",") 

# remove 2017 to model separately
resource.long.17<-dplyr::filter(resource.long, Year == "2017")
resource.long.17$Treatment<-gsub("amwt", "EWS", resource.long.17$Treatment)
resource.long.17$Treatment<-gsub("smwt", "EWS", resource.long.17$Treatment)
resource.long.17$Treatment<-as.factor(resource.long.17$Treatment)
resource.long.18.19<-dplyr::filter(resource.long, Year == "2018" | Year == "2019")

options(tibble.print_max = Inf)
cover_summary<-resource.long.18.19  %>% group_by(cover_type, Treatment, .add = TRUE) %>%
  summarise(count = n(),
            mean = mean(cover, na.rm = TRUE),
            se = sd(cover, na.rm = TRUE)/sqrt(count),
            median = median(cover, na.rm = TRUE))
write.table(cover_summary, "resource_cover_summary18and19.csv", sep=",")
cover_summary<-resource.long.17  %>% group_by(cover_type, Treatment, Year, .add = TRUE) %>%
  summarise(count = n(),
            mean = mean(cover, na.rm = TRUE),
            se = sd(cover, na.rm = TRUE)/sqrt(count),
            median = median(cover, na.rm = TRUE))
write.table(cover_summary, "resource_cover_summary17.csv", sep=",")

## Model creation
### Modelling establishment year separately
#### 2017
resource.max.17 = lmer(LOGIT_Cover~cover_type*Treatment+(1|Site), 
                          data=resource.long.17)
summary(resource.max.17) 
Anova(resource.max.17)
resource.null.17 = lmer(LOGIT_Cover~cover_type+Treatment+(1|Site), 
                           data=resource.long.17)
summary(resource.null.17) 
anova(resource.max.17, resource.null.17)

library(lsmeans)
pairwise.test<-lsmeans(resource.max.17, pairwise ~ Treatment | cover_type, pbkrtest.limit = 3840) 
x<-as.data.frame(summary(pairwise.test))
write.table(x, "pairwise_resource17.csv", sep=",")

# 2018 & 2019
resource.max.18.19 = lmer(LOGIT_Cover~cover_type*Treatment*Year+(1|Site), 
                    data=resource.long.18.19)
summary(resource.max.18.19) 
resource.null.18.19 = lmer(LOGIT_Cover~cover_type+Treatment*Year+(1|Site), 
                          data=resource.long.18.19)
summary(resource.null.18.19) 
Anova(resource.null.18.19)
anova(resource.max.18.19, resource.null.18.19)
resource.null.18.19.2 = lmer(LOGIT_Cover~cover_type*Treatment+Year+(1|Site), 
                           data=resource.long.18.19)
summary(resource.null.18.19.2) 
anova(resource.max.18.19, resource.null.18.19.2)

pairwise.test<-lsmeans(resource.max.18.19, pairwise ~ Treatment | cover_type, pbkrtest.limit = 3840) 
x<-as.data.frame(summary(pairwise.test))
write.table(x, "pairwise_resource18and19.csv", sep=",")
pairwise.test<-lsmeans(resource.max.18.19, pairwise ~ Treatment | cover_type | Year, pbkrtest.limit = 3840) 
x<-as.data.frame(summary(pairwise.test))
write.table(x, "pairwise_resource18_19.csv", sep=",")
pairwise.test<-lsmeans(resource.max.18.19, pairwise ~ Year | cover_type | Treatment, pbkrtest.limit = 3840) 
x<-as.data.frame(summary(pairwise.test))
write.table(x, "pairwise_resource18_19Year.csv", sep=",")

resource.TxYsummary<-group_by(resource.long, cover_type, Year, Treatment) %>%
  summarise(count = n(),
            mean = mean(cover, na.rm = TRUE),
            sd = sd(cover, na.rm = TRUE),
            se = sd(cover, na.rm = TRUE)/sqrt(count),
            median = median(cover, na.rm = TRUE),
            IQR = IQR(cover, na.rm = TRUE))
write.table(resource.TxYsummary, file="resource_TxYsummary.csv", sep=",")



# Export references
setwd("C:/Users/moca1/OneDrive - University of Worcester/SpanishOrangeGrovesV1/papers/Chapter 3 - Establishment and performance of wildflower strips")
write.bib(citation("mvabund"), file = "mvabund.bib", append = FALSE, verbose = TRUE)
write.bib(citation("multcomp"), file = "multcomp.bib", append = FALSE, verbose = TRUE)
write.bib(citation("car"), file = "car.bib", append = FALSE, verbose = TRUE)
write.bib(citation("lme4"), file = "lme4.bib", append = FALSE, verbose = TRUE)
write.bib(citation("permute"), file = "R.bib", append = FALSE, verbose = TRUE)
write.bib(citation("ggplot2"), file = "R.bib", append = FALSE, verbose = TRUE)
write.bib(citation("lsmeans"), file = "R.bib", append = FALSE, verbose = TRUE)
write.bib(citation("RColorBrewer"), file = "R.bib", append = FALSE, verbose = TRUE)
write.bib(citation("emmeans"), file = "R.bib", append = FALSE, verbose = TRUE)
