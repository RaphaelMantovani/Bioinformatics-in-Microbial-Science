library(dplyr)
library(ggplot2)
library(readxl)
library(writexl)

uw_samples <- read_excel("C:/Users/Raphael/Documents/suen_lab/cow-data/feTraitsUW.xlsx")

# filtering for outliers  & adding reference columns

uw_samples$dateSampling <- as.Date(as.character(uw_samples$dateSampling))
uw_samples$dateTrial <- as.Date(as.character(uw_samples$dateTrial))


uw_samples2 <- mutate(uw_samples, rfi_scale = rfi <= 0, 
         DMI_zscore = (DMI - mean(DMI))/sd(DMI), 
         milkE_zscore = (milkE - mean(milkE))/sd(milkE)) %>% 
  filter(!(abs(milkE_zscore) > 3 | abs(DMI_zscore) > 3)) %>% 
  group_by(RFID) %>%
  slice(which.max(dateSampling))
  

attach(uw_samples2)

# recoding the rfi_scale to rfi classification

uw_samples2$rfi_scale <- recode(as.factor(rfi_scale), 'TRUE' = 'negative or null', 'FALSE' = 'positive')

# plotting with original rfi classification

ggplot(uw_samples2, aes(milkE, DMI, color = rfi_scale))+
  geom_point(alpha = .5, size = .8)+
  geom_smooth(method = lm, se = F, color = 'black', size = .6)+
  theme_bw()+
  theme(aspect.ratio = 1)+
  ggtitle('DMI x milkE for UW samples')

# finding the regression line's coefficients and residuals

model <- lm(formula = DMI ~ milkE, data = uw_samples2)

uw_samples2['residuals'] <- model$residuals

# 2.552  is the residual standard error checked with summary(model)

attach(uw_samples2)

# adding efficiency classification 
uw_samples2$classification = ifelse(test = residuals < -2.552, 
    yes = 'highly efficient',
     no = ifelse(test = residuals >= -2.552 & residuals  < 0,
                 yes = 'mildly efficient',
                  no = ifelse(test = 0 <= residuals & residuals <= 2.552,
                             yes = 'mildly inefficient',
                              no = ifelse(residuals > 2.552,
                                         yes = 'highly inefficient',
                                         no = NA))))

# plotting efficiency graph

ggplot(uw_samples2, aes(milkE, DMI, color = classification))+
  geom_point(alpha = .5, size = .8)+
  geom_smooth(method = lm, se = F, color = 'black', size = .6)+
  theme_bw()+
  theme(aspect.ratio = 1)+
  ggtitle('DMI x milkE for UW samples')      


# saving new data frame and exporting to excel 

uw_samples3 <- mutate(uw_samples, rfi_scale = rfi <= 0, 
         DMI_zscore = (DMI - mean(DMI))/sd(DMI), 
         milkE_zscore = (milkE - mean(milkE))/sd(milkE)) %>% 
  filter(!(abs(milkE_zscore) > 3 | abs(DMI_zscore) > 3))

uw_samples3 <- select(uw_samples3, !c(rfi_scale, DMI_zscore, milkE_zscore))
uw_samples3$classification <- uw_samples2$classification

write_xlsx(uw_samples3, "C:\\Users\\Raphael\\Documents\\suen_lab\\cow-data\\feTraitsUW2.xlsx")
