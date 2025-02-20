##Load necessary R packages##
library(haven) #Load sav file
library(intsvy) #For PISA/PV Data
library(tidyverse)
library(dplyr)
library(knitr) #Base for kableExtra
library(kableExtra) #Produce neat tables
library(car) #For VIF calculation
library(corrplot) #For plotting correlation
library(ggplot2)
library(patchwork) #Combine plots
library(quantreg) #Quantile Regression
##Edit Functions from intsvy Package##
#Modify Function to calculate correlation - SE is taken out from the output to allow corrplot
pisa.rho_modified <-
  function(variables, by, data, export=FALSE, name= "output", folder=getwd()) {
    intsvy.rho_modified(variables=variables, by=by, data=data, export=export, name=name,
                        folder=folder, config = pisa_conf)
  }
intsvy.rho_modified <-
  function(variables, by, data, export=FALSE, name= "output", folder=getwd(), config) {
    rho.input <- function(variables, data, config) {
      # BRR / JK
      if (config$parameters$weights == "BRR") {
        # balanced repeated replication
        # Replicate weighted %s (sampling error)
        # in PISA / PIAAC
        weights <- grep("^W_.*[0-9]+$", names(data), value = TRUE)
        data <- na.omit(data[c(variables, config$variables$weightFinal,
                               grep(config$variables$weightBRR, names(data), value=TRUE))])
        # Fifth element is correlation matrix
        rhorp <- lapply(1:config$parameters$BRRreps, function(i) cov.wt(data[variables], wt=
                                                                          data[[weights[i]]], cor = TRUE)[[5]])
        rhotot <- cov.wt(data[variables], wt=data[[config$variables$weightFinal]],
                         cor=TRUE)[[5]]
        # SE formula
        #SE IS DELETED TO ALLOW CORRPLOT
        # Standard error (sampling eror)
        colnames(rhotot) <- unlist(lapply(1:length(variables), function(x)
          c(paste(variables, "Rho", sep=" ")[x])))
        return(round(rhotot, 6))
        
        52
        
      }
      if (config$parameters$weights == "mixed_piaac") {
        # mixed design, different for different coutnries
        # PIAAC
        stop("Not implemented yet")
      }
    }
    # If by no supplied, calculate for the complete sample
    if (missing(by)) {
      output <- rho.input(variables=variables, data=data, config=config)
    } else {
      for (i in by) {
        data[[c(i)]] <- as.factor(data[[c(i)]])
      }
      output <- ddply(data, by, function(x) rho.input(data=x, variables=variables,
                                                      config=config))
    }
    if (export) {
      write.csv(output, file=file.path(folder, paste(name, ".csv", sep="")))
    }
    class(output) <- c("intsvy.rho", class(output))
    return(output)
  }
#Modify regression function to allow for interaction
pisa.reg.pv_modified <-
  function(x, formula, pvlabel, by, data, export=FALSE, name= "output", folder=getwd(),
           std=FALSE) {
    intsvy.reg.pv_modified(x=x, formula = formula, pvnames = pvlabel, by=by, data=data,
                           std=std, export=export,
                           name= name, folder=folder, config=pisa_conf)
  }
#Replaced the lapply function by predetermined formula to include interaction term
intsvy.reg.pv_modified <-
  function(x, formula, pvnames, by, data, std=FALSE, export=FALSE, name= "output",
           folder=getwd(), config) {
    # Remove missing data in IVs
    data <- data[complete.cases(data[, x]), ]
    reg.pv.input <- function(x, pvnames, data, std, config) {
      if (any(sapply(data[x], function(i) all(duplicated(i))))) {
        results <- list("replicates"=NA, "residuals"= NA, "var.w"=NA, "var.b"=NA, "reg"=NA)
        return(results)
      }
      
      53
      
      # BRR / JK
      if (config$parameters$weights == "BRR") {
        # balanced repeated replication
        # Replicate weighted %s (sampling error)
        # in PISA
        #pvnames <- paste0(pvnames, ".*[0-9]|[0-9].*", pvnames)
        #pvnames <- grep(pvnames, names(data), value = TRUE)
        weights <- grep(paste0("^", config$variables$weightBRR , ".*[0-9]+$"),
                        names(data), value = TRUE)
        # remove missings in pvalues and weights
        data <- data[complete.cases(data[, c(pvnames[1], weights[1],
                                             config$variables$weightFinal)]), ]
        # List of formulas for each PV
        ###I CHANGED THIS PART ###
        regform <- lapply(paste0("PV",1:10,"MATH"), function(i) paste(i, "~", formula))
        # Standardise IV and DV variables
        if(std) {
          data <- cbind(scale(data[c(pvnames, x)]), data[!names(data) %in% c(pvnames, x)])
        }
        # Replicate weighted coefficients for sampling error (PVs)
        reg.rep <- lapply(regform, function(pv) lapply(1:length(weights), function(rep)
          summary(lm(formula=as.formula(pv), data=data, weights=data[[weights[rep]]]))))
        # Combining coefficients and R-squared replicates
        coe.rep <- lapply(1:length(pvnames), function(pv) sapply(1:length(weights),
                                                                 function(rep)
                                                                   c(reg.rep[[pv]][[rep]]$coefficients[,1], "R-squared"= reg.rep[[pv]][[rep]]$r.squared)))
        resid <- lapply(1:length(pvnames), function(pv)
          sapply(1:length(weights),
                 function(rep) reg.rep[[pv]][[rep]]$residuals))
        # Total weighted coefficient for each PV for imputation (between) error
        reg.pv <- lapply(regform, function(pv)
          summary(lm(formula=as.formula(pv), data=data,
                     weights=data[[config$variables$weightFinal]])))
        coe.tot <- sapply(1:length(pvnames), function(pv)
          c(reg.pv[[pv]]$coefficients[, 1], "R-squared" = reg.pv[[pv]]$r.squared))
        
        # Mean total coefficients (across PVs)
        stat.tot <- apply(coe.tot, 1, mean)
        
        54
        
        # Sampling error (variance within)
        cc = 1/(length(weights)*(1-0.5)^2)
        var.w <- apply(cc*sapply(lapply(1:length(pvnames), function(pv)
          (coe.rep[[pv]]-coe.tot[,pv])^2), function(e) apply(e, 1, sum)), 1, mean)
        # Imputation error (variance between)
        var.b <- (1/(length(pvnames)-1))*apply(sapply(1:length(pvnames), function(pv)
          (coe.tot[, pv] - stat.tot)^2), 1, sum)
        stat.se <- (var.w +(1+1/length(pvnames))*var.b)^(1/2)
        stat.t <- stat.tot/stat.se
        # Reg Table
        reg.tab <- data.frame("Estimate"=stat.tot, "Std. Error"=stat.se, "t value"=stat.t,
                              check.names=F)
        results <- list("replicates"=lapply(coe.rep, t), "residuals"= resid, "var.w"=var.w,
                        "var.b"=var.b, "reg"=reg.tab)
        return(results)
      }
    }
    # If by no supplied, calculate for the complete sample
    if (missing(by)) {
      output <- reg.pv.input(x=x, pvnames=pvnames, data=data, std=std, config=config)
    } else {
      output <- lapply(split(data, droplevels(data[by])), function(i)
        reg.pv.input(x=x, pvnames=pvnames, data=i, std=std, config=config))
    }
    if (export) {
      write.csv(do.call(rbind, lapply(output, function(x) x$reg)), file=file.path(folder,
                                                                                  paste(name, ".csv", sep="")))
    }
    class(output) <- "intsvy.reg"
    return(output)
  }
##Load data##
student_data <- read_sav("CY08MSP_STU_QQQ.sav")
school_data <- read_sav("CY08MSP_SCH_QQQ.sav")
##Data filtering and cleaning##
var_student <-
  c("CNT","CNTSCHID","CNTSTUID","ST004D01T","FISCED","MISCED","ESCS","ICTS
CH","ICTHOME","ICTRES","ICTQUAL","W_FSTUWT","SENWT")
school_data["Prop_Girl"] = school_data["SC002Q02TA"] / (school_data["SC002Q02TA"] +
                                                          school_data["SC002Q01TA"])

55

var_school <-
  c("CNTSCHID","MCLSIZE","STRATIO","MTTRAIN","PROATCE","Prop_Girl","SC013
Q01TA")
rep_weights <- paste0("W_FSTURWT",1:80)
math_pv <- paste0("PV",1:10,"MATH")
variables <- c(var_student,rep_weights,math_pv)
student_subset <- student_data[variables]
school_subset <- school_data[var_school]
subset_data <- merge(student_subset,school_subset,by = "CNTSCHID")
subset_data <- subset_data %>%
  rename(
    School_ID = CNTSCHID,
    Country = CNT,
    Student_ID = CNTSTUID,
    Gender = ST004D01T,
    FatherEdu = FISCED,
    MotherEdu = MISCED,
    Class_Size = MCLSIZE,
    Math_Training = MTTRAIN,
    Certified_Teacher = PROATCE,
    Ownership = SC013Q01TA
  )
subset_data$Gender <- subset_data$Gender - 1
subset_data$Ownership <- subset_data$Ownership - 1
data <- subset_data
##Data Imputation##
#Flag missing values
data$Gender_flag <- ifelse(is.na(data$Gender),data$Gender_flag <- 1,data$Gender_flag <-
                             0)
data$FatherEdu_flag <- ifelse(is.na(data$FatherEdu),data$FatherEdu_flag <-
                                1,data$FatherEdu_flag <- 0)
data$MotherEdu_flag <- ifelse(is.na(data$MotherEdu),data$MotherEdu_flag <-
                                1,data$MotherEdu_flag <- 0)
data$ESCS_flag <- ifelse(is.na(data$ESCS),data$ESCS_flag <- 1,data$ESCS_flag <- 0)
data$ICTHOME_flag <- ifelse(is.na(data$ICTHOME),data$ICTHOME_flag <-
                              1,data$ICTHOME_flag <- 0)
data$ICTSCH_flag <- ifelse(is.na(data$ICTSCH),data$ICTSCH_flag <-
                             1,data$ICTSCH_flag <- 0)
data$ICTRES_flag <- ifelse(is.na(data$ICTRES),data$ICTRES_flag <- 1,data$ICTRES_flag
                           <- 0)
data$ICTQUAL_flag <- ifelse(is.na(data$ICTQUAL),data$ICTQUAL_flag <-
                              1,data$ICTQUAL_flag <- 0)
data$Class_Size_flag <- ifelse(is.na(data$Class_Size),data$Class_Size_flag <-
                                 1,data$Class_Size_flag <- 0)
data$STRATIO_flag <- ifelse(is.na(data$STRATIO),data$STRATIO_flag <-
                              1,data$STRATIO_flag <- 0)

56

data$Math_Training_flag <- ifelse(is.na(data$Math_Training),data$Math_Training_flag <-
                                    1,data$Math_Training_flag <- 0)
data$Certified_Teacher_flag <-
  ifelse(is.na(data$Certified_Teacher),data$Certified_Teacher_flag <-
           1,data$Certified_Teacher_flag <- 0)
data$Prop_Girl_flag <- ifelse(is.na(data$Prop_Girl),data$Prop_Girl_flag <-
                                1,data$Prop_Girl_flag <- 0)
data$Ownership_flag <- ifelse(is.na(data$Ownership),data$Ownership_flag <-
                                1,data$Ownership_flag <- 0)
#Replace NA in Categorical
data$Gender[is.na(data$Gender)] <- 99
data$Ownership[is.na(data$Ownership)] <- 99
#Replace NA for Continuous
escs_impute <- pisa.mean(variable = "ESCS", data = data, by = "Country") %>%
  select(Country, Mean)
names(escs_impute) <- c("Country","ESCS_mean")
father_impute <- pisa.mean(variable = "FatherEdu", data = data, by = "Country") %>%
  select(Country, Mean)
names(father_impute) <- c("Country","Father_mean")
mother_impute <- pisa.mean(variable = "MotherEdu", data = data, by = "Country") %>%
  select(Country, Mean)
names(mother_impute) <- c("Country","Mother_mean")
ICTHOME_impute <- pisa.mean(variable = "ICTHOME", data = data, by =
                              "Country") %>%
  select(Country, Mean)
names(ICTHOME_impute) <- c("Country","ICTHOME_mean")
ICTSCH_impute <- pisa.mean(variable = "ICTSCH", data = data, by = "Country") %>%
  select(Country, Mean)
names(ICTSCH_impute) <- c("Country","ICTSCH_mean")
ICTRES_impute <- pisa.mean(variable = "ICTRES", data = data, by = "Country") %>%
  select(Country, Mean)
names(ICTRES_impute) <- c("Country","ICTRES_mean")
ICTQUAL_impute <- pisa.mean(variable = "ICTQUAL", data = data, by = "Country") %>%
  select(Country, Mean)
names(ICTQUAL_impute) <- c("Country","ICTQUAL_mean")
Class_Size_impute <- pisa.mean(variable = "Class_Size", data = data, by = "Country") %>%
  select(Country, Mean)
names(Class_Size_impute) <- c("Country","Class_Size_mean")
STRATIO_impute <- pisa.mean(variable = "STRATIO", data = data, by = "Country") %>%
  select(Country, Mean)

57

names(STRATIO_impute) <- c("Country","STRATIO_mean")
Math_Training_impute <- pisa.mean(variable = "Math_Training", data = data, by =
                                    "Country") %>%
  select(Country, Mean)
names(Math_Training_impute) <- c("Country","Math_Training_mean")
Certified_Teacher_impute <- pisa.mean(variable = "Certified_Teacher", data = data, by =
                                        "Country") %>%
  select(Country, Mean)
names(Certified_Teacher_impute) <- c("Country","Certified_Teacher_mean")
Prop_Girl_impute <- pisa.mean(variable = "Prop_Girl", data = data, by = "Country") %>%
  select(Country, Mean)
names(Prop_Girl_impute) <- c("Country","Prop_Girl_mean")
data <- merge(data, escs_impute, by="Country")
data <- merge(data, father_impute, by="Country")
data <- merge(data, mother_impute, by="Country")
data <- merge(data, ICTHOME_impute, by="Country")
data <- merge(data, ICTSCH_impute, by="Country")
data <- merge(data, ICTRES_impute, by="Country")
data <- merge(data, ICTQUAL_impute, by="Country")
data <- merge(data, Class_Size_impute, by="Country")
data <- merge(data, STRATIO_impute, by="Country")
data <- merge(data, Math_Training_impute, by="Country")
data <- merge(data, Certified_Teacher_impute, by="Country")
data <- merge(data, Prop_Girl_impute, by="Country")
data <- data %>%
  mutate(ESCS = ifelse(is.na(ESCS), ESCS_mean, ESCS)) %>%
  mutate(FatherEdu = ifelse(is.na(FatherEdu), Father_mean, FatherEdu)) %>%
  mutate(MotherEdu = ifelse(is.na(MotherEdu), Mother_mean, MotherEdu)) %>%
  mutate(ICTHOME = ifelse(is.na(ICTHOME), ICTHOME_mean, ICTHOME)) %>%
  mutate(ICTSCH = ifelse(is.na(ICTSCH), ICTSCH_mean, ICTSCH)) %>%
  mutate(ICTRES = ifelse(is.na(ICTRES), ICTRES_mean, ICTRES)) %>%
  mutate(ICTQUAL = ifelse(is.na(ICTQUAL), ICTQUAL_mean, ICTQUAL)) %>%
  mutate(Class_Size = ifelse(is.na(Class_Size), Class_Size_mean, Class_Size)) %>%
  mutate(STRATIO = ifelse(is.na(STRATIO), STRATIO_mean, STRATIO)) %>%
  mutate(Math_Training = ifelse(is.na(Math_Training), Math_Training_mean,
                                Math_Training)) %>%
  mutate(Certified_Teacher = ifelse(is.na(Certified_Teacher), Certified_Teacher_mean,
                                    Certified_Teacher)) %>%
  mutate(Prop_Girl = ifelse(is.na(Prop_Girl), Prop_Girl_mean, Prop_Girl)) %>%
  select(-ESCS_mean, -Father_mean, -Mother_mean, -ICTHOME_mean, -ICTSCH_mean, -
           ICTRES_mean, -ICTQUAL_mean, -Class_Size_mean, -STRATIO_mean, -
           Math_Training_mean, -Certified_Teacher_mean, -Prop_Girl_mean)
##Omit rows with NAs
complete_data <- na.omit(data)

58

##Descriptive Stats##
#Calculate mean
math_PV <- pisa.mean.pv(pvlabel = paste0("PV",1:10,"MATH"),data = complete_data)
gender_mean <- pisa.mean(variable = "Gender", data =
                           complete_data[complete_data["Gender"]!=99,])
father_mean <- pisa.mean(variable = "FatherEdu", data = complete_data)
mother_mean <- pisa.mean(variable = "MotherEdu", data = complete_data)
escs_mean <- pisa.mean(variable = "ESCS", data = complete_data)
icthome_mean <- pisa.mean(variable = "ICTHOME", data = complete_data)
ictsch_mean <- pisa.mean(variable = "ICTSCH", data = complete_data)
ictres_mean <- pisa.mean(variable = "ICTRES", data = complete_data)
ictqual_mean <- pisa.mean(variable = "ICTQUAL", data = complete_data)
classsize_mean <- pisa.mean(variable = "Class_Size", data = complete_data)
stratio_mean <- pisa.mean(variable = "STRATIO", data = complete_data)
mathtrain_mean <- pisa.mean(variable = "Math_Training", data = complete_data)
certify_mean <- pisa.mean(variable = "Certified_Teacher", data = complete_data)
girlprop_mean <- pisa.mean(variable = "Prop_Girl", data = complete_data)
ownership_mean <- pisa.mean(variable = "Ownership", data =
                              complete_data[complete_data["Ownership"]!=99,])
#Combine means together to form one table
table_comb <- round(rbind(math_PV,gender_mean, father_mean, mother_mean,
                          escs_mean,classsize_mean,stratio_mean,mathtrain_mean,certify_mean,girlprop_mean,owner
                          ship_mean),2)
#Change table labels
rownames(table_comb) <- c("Math Performance","Gender", "Father's Education","Mother's
Education","Economic, Social and Cultural Status","Maths Class Size", "Student-Teacher
Ratio", "Maths Professional Training", "Certified Teacher","Gender
Distribution","Ownership")
names(table_comb) <- c("Frequency","Mean","SE of Mean","Standard Deviation","SE of
Standard Deviation")
table_comb <- data.frame(Rowname = rownames(table_comb),table_comb)
#Calculate max and min for categorical variables
continuous_subset <- complete_data[,c("ESCS","FatherEdu","MotherEdu","Class_Size",
                                      "STRATIO", "Math_Training", "Certified_Teacher","Prop_Girl")]
cont_max <- data.frame(sapply(continuous_subset, max, na.rm = T))
cont_min <- data.frame(sapply(continuous_subset, min, na.rm = T))
continuous_stats <- cbind(cont_min,cont_max)
row.names(continuous_stats) <- c("Economic, Social and Cultural Status","Father's
Education","Mother's Education","Maths Class Size", "Student-Teacher Ratio", "Maths
Professional Training", "Certified Teacher","Gender Distribution")
names(continuous_stats) <- c("Min","Max")
continuous_stats <- data.frame(Rowname = rownames(continuous_stats),continuous_stats)
stats_table <- merge(table_comb, continuous_stats, by = "Rowname", all = TRUE)
flag_subset <-
  complete_data[,c("Gender_flag","ESCS_flag","FatherEdu_flag","MotherEdu_flag","Class_S
ize_flag", "STRATIO_flag", "Math_Training_flag",
                   "Certified_Teacher_flag","Prop_Girl_flag","Ownership_flag")]
flag_count <- data.frame(sapply(flag_subset, sum, na.rm = T))

59

rownames(flag_count) <- c("Gender", "Father's Education","Mother's Education","Economic,
Social and Cultural Status","Maths Class Size", "Student-Teacher Ratio", "Maths
Professional Training", "Certified Teacher","Gender Distribution","Ownership")
names(flag_count) <- c("Number of Imputed Values")
flag_count <- data.frame(Rowname = rownames(flag_count),flag_count)
stats_table <- merge(stats_table, flag_count, by = "Rowname", all = TRUE)
#Reorder table
row_order <- c(6,4,2,3,9,7,11,8,1,5,10)
stats_table <- stats_table[row_order,]
#Calculate max and min for PV
stats_table[1,7] <-
  min(apply(complete_data[,c("PV1MATH","PV2MATH","PV3MATH","PV4MATH","PV5
MATH","PV6MATH","PV7MATH","PV8MATH","PV9MATH","PV10MATH")],1,min))
stats_table[1,8] <-
  max(apply(complete_data[,c("PV1MATH","PV2MATH","PV3MATH","PV4MATH","PV5
MATH","PV6MATH","PV7MATH","PV8MATH","PV9MATH","PV10MATH")],1,max))
#Fill min/max with 0/1 for binary variables
stats_table[2,7] <- 0
stats_table[2,8] <- 1
stats_table[11,7] <- 0
stats_table[11,8] <- 1
rownames(stats_table) <- NULL
names(stats_table) <- c("Variable","Frequency","Mean","SE of Mean","Standard
Deviation","SE of Standard Deviation","Min","Max","Number of Imputed Values")
stats_table[c("Min","Max")] <- round(stats_table[c("Min","Max")],2)
statsvar_subset <- c("Variable","Frequency","Mean","Standard
Deviation","Min","Max","Number of Imputed Values")
stats_table <- stats_table[statsvar_subset]
stats_table %>%
  kbl() %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = F)
kable(stats_table, format = "latex", booktabs = TRUE) %>% kable_styling(latex_options =
                                                                          c("striped","scale_down"))
##Generate Plot##
#Graph for ESCS
escs_plot <- ggplot(complete_data, aes(x = ESCS, weight = W_FSTUWT)) +
  geom_density(fill = "blue", alpha = 0.2) +
  geom_vline(xintercept = -0.49, color = "red", linetype = "dashed", size = 1) +
  labs(x = "Economic, Social and Cultural Status", y = "Density") +
  theme_minimal()
#Graph for Father's Education
father_plot <- ggplot(complete_data, aes(x = FatherEdu, weight = W_FSTUWT)) +
  geom_density(fill = "blue", alpha = 0.2) +
  geom_vline(xintercept = 5.79, color = "red", linetype = "dashed", size = 1) +
  labs(x = "Father's Education", y = "Density") +
  
  60

scale_x_continuous(limits = c(1,10) )
theme(plot.title = element_text(hjust = 0.5))
#Graph for Father's Education
mother_plot <- ggplot(complete_data, aes(x = MotherEdu, weight = W_FSTUWT)) +
  geom_density(fill = "blue", alpha = 0.2) +
  geom_vline(xintercept = 5.93, color = "red", linetype = "dashed", size = 1) +
  labs(x = "Mother's Education", y = "Density") +
  scale_x_continuous(limits = c(1,10) )
theme(plot.title = element_text(hjust = 0.5))
#Graph for Maths Class Size
class_plot <- ggplot(complete_data, aes(x = Class_Size, weight = W_FSTUWT)) +
  geom_density(fill = "blue", alpha = 0.2) +
  geom_vline(xintercept = 29.34, color = "red", linetype = "dashed", size = 1) +
  labs(x = "Maths Class Size", y = "Density") +
  scale_x_continuous(limits = c(10,55) )
theme(plot.title = element_text(hjust = 0.5))
#Graph for Student-Teacher Ratio
stratio_plot <- ggplot(complete_data, aes(x = STRATIO, weight = W_FSTUWT)) +
  geom_density(fill = "blue", alpha = 0.2) +
  geom_vline(xintercept = 17.56, color = "red", linetype = "dashed", size = 1) +
  labs(x = "Student-Teacher Ratio", y = "Density") +
  scale_x_continuous(limits = c(1,100) )
theme(plot.title = element_text(hjust = 0.5))
#Graph for Maths Professional Training
train_plot <- ggplot(complete_data, aes(x = Math_Training, weight = W_FSTUWT)) +
  geom_density(fill = "blue", alpha = 0.2) +
  geom_vline(xintercept = 0.29, color = "red", linetype = "dashed", size = 1) +
  labs(x = "Maths Professional Training", y = "Density") +
  scale_x_continuous(limits = c(-1.8,1.5) )
theme(plot.title = element_text(hjust = 0.5))
#Graph for Certified Teacher
certified_plot <- ggplot(complete_data, aes(x = Certified_Teacher, weight = W_FSTUWT))
+ geom_density(fill = "blue", alpha = 0.2) +
  geom_vline(xintercept = 0.86, color = "red", linetype = "dashed", size = 1) +
  labs(x = "Proportion of Certified Teacher", y = "Density") +
  scale_x_continuous(limits = c(0,1) )
theme(plot.title = element_text(hjust = 0.5))
#Graph for Girl's Proportion
girl_plot <- ggplot(complete_data, aes(x = Prop_Girl, weight = W_FSTUWT)) +
  geom_density(fill = "blue", alpha = 0.2) +
  geom_vline(xintercept = 0.5, color = "red", linetype = "dashed", size = 1) +
  labs(x = "Proportion of Girls", y = "Density") +
  scale_x_continuous(limits = c(0,1) )
theme(plot.title = element_text(hjust = 0.5))

61

#Graph for ICTHOME
icthome_plot <- ggplot(complete_data, aes(x = ICTHOME, weight = W_FSTUWT)) +
  geom_density(fill = "blue", alpha = 0.2) +
  geom_vline(xintercept = -0.08156, color = "red", linetype = "dashed", size = 1) +
  labs(x = "ICT at Home", y = "Density") +
  scale_x_continuous(limits = c(-6.5, 0.5)) +
  theme(plot.title = element_text(hjust = 0.5))
#Graph for ICTSCH
ictsch_plot <- ggplot(complete_data, aes(x = ICTSCH, weight = W_FSTUWT)) +
  geom_density(fill = "blue", alpha = 0.2) +
  geom_vline(xintercept = -0.045230, color = "red", linetype = "dashed", size = 1) +
  labs(x = "ICT at School", y = "Density") +
  scale_x_continuous(limits = c(-5.7, 0.65)) +
  theme(plot.title = element_text(hjust = 0.5))
#Graph for ICTQUAL
ictqual_plot <- ggplot(complete_data, aes(x = ICTQUAL, weight = W_FSTUWT)) +
  geom_density(fill = "blue", alpha = 0.2) +
  geom_vline(xintercept = -0.09313, color = "red", linetype = "dashed", size = 1) +
  labs(x = "ICT Quality", y = "Density") +
  scale_x_continuous(limits = c(-3, 3)) +
  theme(plot.title = element_text(hjust = 0.5))
#Graph for ICTRES
ictres_plot <-ggplot(complete_data, aes(x = ICTRES, weight = W_FSTUWT)) +
  geom_density(fill = "blue", alpha = 0.2) +
  geom_vline(xintercept = -0.49373, color = "red", linetype = "dashed", size = 1) +
  labs(x = "ICT Resource", y = "Density") +
  scale_x_continuous(limits = c(-6.2, 5.4)) +
  theme(plot.title = element_text(hjust = 0.5))
#Combine density plots
escs_plot
main_plot <- father_plot + mother_plot + plot_layout(ncol = 2)
main_plot
ict_plot <- icthome_plot + ictsch_plot + ictqual_plot + ictres_plot + plot_layout(ncol = 2)
ict_plot
school_plot <- class_plot + stratio_plot + train_plot + certified_plot + girl_plot +
  plot_layout(ncol = 2)
school_plot
#Calculate mean PISA score by country and plot the top and bottom 10
pmeans <- pisa.mean.pv(pvlabel=paste0("PV",1:10,"MATH"), by="Country",
                       data=complete_data, export=FALSE)
pmeans <-sort(pmeans,by = "Mean",ascending=FALSE)
plot(pmeans,sort=TRUE) + theme(plot.title = element_text(hjust = 0.5))
##Create Categorical ICT Variables##

62

#Create Binary Variables based on proportion for interaction
summary(complete_data$ICTHOME)
quantile(complete_data$ICTHOME,p=c(0.25,0.5,0.75))
#ICT Home
sum(complete_data$ICTHOME>=0.3346,na.rm = TRUE) /
  (sum(complete_data$ICTHOME>=0.3346,na.rm = TRUE) +
     sum(complete_data$ICTHOME<0.3346,na.rm = TRUE))
#0.6687492
quantile(complete_data$ICTHOME,p=c(0.3))
#0.0978 - 30 percentile
complete_data$ICTHOME_binary <- ifelse(complete_data$ICTHOME >= 0.0978, 1, 0)
complete_data$ICTHOME_binary <- ifelse(complete_data$ICTHOME_flag == 1, 99,
                                       complete_data$ICTHOME_binary)
#ICT School
summary(complete_data$ICTSCH)
quantile(complete_data$ICTSCH,p=c(0.25,0.5,0.75))
sum(complete_data$ICTSCH>=0.4062,na.rm = TRUE) /
  (sum(complete_data$ICTSCH>=0.4062,na.rm = TRUE) +
     sum(complete_data$ICTSCH<0.4062,na.rm = TRUE))
#0.6268301
quantile(complete_data$ICTSCH,p=c(0.3))
#0.04752493
complete_data$ICTSCH_binary <- ifelse(complete_data$ICTSCH >= 0.04752493, 1, 0)
complete_data$ICTSCH_binary <- ifelse(complete_data$ICTSCH_flag == 1, 99,
                                      complete_data$ICTSCH_binary)
#ICT Quality
quantile(complete_data$ICTQUAL,p=c(0.333,0.667))
#-0.4775,0.1515
complete_data$ICTQUAL_cat <- ifelse(complete_data$ICTQUAL >= 0.1515, 2,
                                    ifelse(complete_data$ICTQUAL < -0.4775,0,1))
complete_data$ICTQUAL_cat <- ifelse(complete_data$ICTQUAL_flag == 1, 99,
                                    complete_data$ICTQUAL_cat)
#ICT Resource
quantile(complete_data$ICTRES,p=c(0.333,0.667))
#-1.0153,0.1020
complete_data$ICTRES_cat <- ifelse(complete_data$ICTRES >= 0.1020, 2,
                                   ifelse(complete_data$ICTRES < -1.0153,0,1))
complete_data$ICTRES_cat <- ifelse(complete_data$ICTRES_flag == 1, 99,
                                   complete_data$ICTRES_cat)
##Main Regression##
#Change variables to factors
complete_data$Gender <- factor(complete_data$Gender)
complete_data$ICTHOME_binary <- factor(complete_data$ICTHOME_binary)
complete_data$ICTSCH_binary <- factor(complete_data$ICTSCH_binary)
complete_data$ICTQUAL_cat <- factor(complete_data$ICTQUAL_cat)

63

complete_data$ICTRES_cat <- factor(complete_data$ICTRES_cat)
complete_data$Ownership <- factor(complete_data$Ownership)
complete_data$Gender_flag <- factor(complete_data$Gender_flag)
complete_data$FatherEdu_flag <- factor(complete_data$FatherEdu_flag)
complete_data$MotherEdu_flag <- factor(complete_data$MotherEdu_flag)
complete_data$ESCS_flag <- factor(complete_data$ESCS_flag)
complete_data$Class_Size_flag <- factor(complete_data$Class_Size_flag)
complete_data$STRATIO_flag <- factor(complete_data$STRATIO_flag)
complete_data$Math_Training_flag <- factor(complete_data$Math_Training_flag)
complete_data$Certified_Teacher_flag <- factor(complete_data$Certified_Teacher_flag)
complete_data$Prop_Girl_flag <- factor(complete_data$Prop_Girl_flag)
#Descrptive Results
bivariate_HOME <- pisa.reg.pv_modified(x=c("ICTHOME_binary"), formula =
                                         "ICTHOME_binary", pvlabel = paste0("PV",1:10,"MATH"),data = complete_data)
bivariate_SCH <- pisa.reg.pv_modified(x=c("ICTSCH_binary"), formula =
                                        "ICTSCH_binary", pvlabel = paste0("PV",1:10,"MATH"),data = complete_data)
bivariate_resource <- pisa.reg.pv_modified(x=c("ICTRES_cat"), formula = "ICTRES_cat",
                                           pvlabel = paste0("PV",1:10,"MATH"),data = complete_data)
bivariate_quality <- pisa.reg.pv_modified(x=c("ICTQUAL_cat"), formula =
                                            "ICTQUAL_cat", pvlabel = paste0("PV",1:10,"MATH"),data = complete_data)
#Visualising Descriptive Results
#bivariate_HOME
bivariate_HOME_df <- data.frame(round(bivariate_HOME$reg,3))
bivariate_HOME_row_order <- c(2,3,1,4)
rownames(bivariate_HOME_df) <- c("Intercept","ICT at Home (High)", "ICT at Home
(Missing)","R-Squared")
bivariate_HOME_df <- bivariate_HOME_df[bivariate_HOME_row_order,]
names(bivariate_HOME_df) <- c("Estimated Coefficient","Standard Error","T Value")
#bivariate_SCH
bivariate_SCH_df <- data.frame(round(bivariate_SCH$reg,3))
bivariate_SCH_row_order <- c(2,3,1,4)
rownames(bivariate_SCH_df) <- c("Intercept","ICT at School (High)", "ICT at School
(Missing)","R-Squared")
bivariate_SCH_df <- bivariate_SCH_df[bivariate_SCH_row_order,]
names(bivariate_SCH_df) <- c("Estimated Coefficient","Standard Error","T Value")
#bivariate_resource
bivariate_resource_df <- data.frame(round(bivariate_resource$reg,3))
bivariate_resource_row_order <- c(2,3,4,1,5)
rownames(bivariate_resource_df) <- c("Intercept","ICT Resource (Medium)","ICT Resource
(High)", "ICT Resource (Missing)","R-Squared")
bivariate_resource_df <- bivariate_resource_df[bivariate_resource_row_order,]
names(bivariate_resource_df) <- c("Estimated Coefficient","Standard Error","T Value")
#bivariate_quality
bivariate_quality_df <- data.frame(round(bivariate_quality$reg,3))
bivariate_quality_row_order <- c(2,3,4,1,5)
rownames(bivariate_quality_df) <- c("Intercept","ICT Quality (Medium)","ICT Quality
(High)", "ICT Quality (Missing)","R-Squared")
bivariate_quality_df <- bivariate_quality_df[bivariate_quality_row_order,]

64

names(bivariate_quality_df) <- c("Estimated Coefficient","Standard Error","T Value")
#Tidy Regression Models
for (i in 1:nrow(bivariate_HOME_df)){
  bivariate_HOME_df$Home_coef[i] <- paste(toString(bivariate_HOME_df$`Estimated
                                                   Coefficient`[i]),"(",toString(bivariate_HOME_df$`Standard Error`[i]),")")
  bivariate_HOME_df$Home_coef[i] <- ifelse((1.645 >= bivariate_HOME_df$`T Value`[i])
                                           & (bivariate_HOME_df$`T Value`[i] >= -1.645), bivariate_HOME_df$Home_coef[i],
                                           paste(bivariate_HOME_df$Home_coef[i],"*"))
  bivariate_HOME_df$Home_coef[i] <- ifelse((-1.96 >= bivariate_HOME_df$`T Value`[i]) |
                                             (bivariate_HOME_df$`T Value`[i] >= 1.96), paste(bivariate_HOME_df$Home_coef[i],"*"),
                                           bivariate_HOME_df$Home_coef[i])
  bivariate_HOME_df$Home_coef[i] <- ifelse((-2.576 >= bivariate_HOME_df$`T Value`[i])
                                           | (bivariate_HOME_df$`T Value`[i] >= 2.576),
                                           paste(bivariate_HOME_df$Home_coef[i],"*"), bivariate_HOME_df$Home_coef[i])
}
for (i in 1:nrow(bivariate_SCH_df)){
  bivariate_SCH_df$School_coef[i] <- paste(toString(bivariate_SCH_df$`Estimated
                                                    Coefficient`[i]),"(",toString(bivariate_SCH_df$`Standard Error`[i]),")")
  bivariate_SCH_df$School_coef[i] <- ifelse((1.645 >= bivariate_SCH_df$`T Value`[i]) &
                                              (bivariate_SCH_df$`T Value`[i] >= -1.645), bivariate_SCH_df$School_coef[i],
                                            paste(bivariate_SCH_df$School_coef[i],"*"))
  bivariate_SCH_df$School_coef[i] <- ifelse((-1.96 >= bivariate_SCH_df$`T Value`[i]) |
                                              (bivariate_SCH_df$`T Value`[i] >= 1.96), paste(bivariate_SCH_df$School_coef[i],"*"),
                                            bivariate_SCH_df$School_coef[i])
  bivariate_SCH_df$School_coef[i] <- ifelse((-2.576 >= bivariate_SCH_df$`T Value`[i]) |
                                              (bivariate_SCH_df$`T Value`[i] >= 2.576), paste(bivariate_SCH_df$School_coef[i],"*"),
                                            bivariate_SCH_df$School_coef[i])
}
for (i in 1:nrow(bivariate_resource_df)){
  bivariate_resource_df$Resource_coef[i] <- paste(toString(bivariate_resource_df$`Estimated
                                                           Coefficient`[i]),"(",toString(bivariate_resource_df$`Standard Error`[i]),")")
  bivariate_resource_df$Resource_coef[i] <- ifelse((1.645 >= bivariate_resource_df$`T
                                                    Value`[i]) & (bivariate_resource_df$`T Value`[i] >= -1.645),
                                                   bivariate_resource_df$Resource_coef[i], paste(bivariate_resource_df$Resource_coef[i],"*"))
  bivariate_resource_df$Resource_coef[i] <- ifelse((-1.96 >= bivariate_resource_df$`T
                                                    Value`[i]) | (bivariate_resource_df$`T Value`[i] >= 1.96),
                                                   paste(bivariate_resource_df$Resource_coef[i],"*"), bivariate_resource_df$Resource_coef[i])
  bivariate_resource_df$Resource_coef[i] <- ifelse((-2.576 >= bivariate_resource_df$`T
                                                    Value`[i]) | (bivariate_resource_df$`T Value`[i] >= 2.576),
                                                   paste(bivariate_resource_df$Resource_coef[i],"*"), bivariate_resource_df$Resource_coef[i])
}
for (i in 1:nrow(bivariate_quality_df)){
  bivariate_quality_df$Quality_coef[i] <- paste(toString(bivariate_quality_df$`Estimated
                                                         Coefficient`[i]),"(",toString(bivariate_quality_df$`Standard Error`[i]),")")
  bivariate_quality_df$Quality_coef[i] <- ifelse((1.645 >= bivariate_quality_df$`T Value`[i])
                                                 & (bivariate_quality_df$`T Value`[i] >= -1.645), bivariate_quality_df$Quality_coef[i],
                                                 paste(bivariate_quality_df$Quality_coef[i],"*"))
  
  65
  
  bivariate_quality_df$Quality_coef[i] <- ifelse((-1.96 >= bivariate_quality_df$`T Value`[i]) |
                                                   (bivariate_quality_df$`T Value`[i] >= 1.96), paste(bivariate_quality_df$Quality_coef[i],"*"),
                                                 bivariate_quality_df$Quality_coef[i])
  bivariate_quality_df$Quality_coef[i] <- ifelse((-2.576 >= bivariate_quality_df$`T Value`[i])
                                                 | (bivariate_quality_df$`T Value`[i] >= 2.576),
                                                 paste(bivariate_quality_df$Quality_coef[i],"*"), bivariate_quality_df$Quality_coef[i])
}
merge1 <- merge(bivariate_HOME_df, bivariate_SCH_df, by = "row.names",all = TRUE)
merge2 <- merge(bivariate_resource_df,bivariate_quality_df,by = "row.names",all = TRUE)
descriptive_df <- merge(merge1,merge2,by = "Row.names", all = TRUE)
descriptive_coefs <-
  select(descriptive_df,Row.names,Home_coef,School_coef,Resource_coef,Quality_coef)
names(descriptive_coefs) <- c("Coefficients","ICT Usage at Home","ICT Usage at
School","ICT Resources","ICT Quality")
descriptive_coefs <- descriptive_coefs[c(1,3,6,5,9,8,11,12),]
kable(descriptive_coefs, format = "latex", booktabs = TRUE) %>%
  kable_styling(latex_options = c("striped","scale_down"))
#Main Model 1
model_1 <-
  pisa.reg.pv_modified(x=c("Gender","ESCS","ICTHOME_binary","ICTSCH_binary",
                           "ICTQUAL_cat","ICTRES_cat","Class_Size","STRATIO","Math_Training","Certified_Teac
her","Prop_Girl","Ownership","ESCS_flag","Class_Size_flag","STRATIO_flag","Math_Trai
ning_flag","Certified_Teacher_flag","Prop_Girl_flag"), formula = "Gender + ESCS +
ICTHOME_binary + ICTSCH_binary + ICTQUAL_cat + ICTRES_cat + Class_Size +
STRATIO + Math_Training + Certified_Teacher + Prop_Girl + Ownership + ESCS_flag +
Class_Size_flag + STRATIO_flag + Math_Training_flag + Certified_Teacher_flag +
Prop_Girl_flag", pvlabel = paste0("PV",1:10,"MATH"),data = complete_data)
#Main Model 2
model_2 <-
  pisa.reg.pv_modified(x=c("Gender","ESCS","ICTHOME_binary","ICTSCH_binary",
                           "ICTQUAL_cat","ICTRES_cat","Class_Size","STRATIO","Math_Training","Certified_Teac
her","Prop_Girl","Ownership","ESCS_flag","Class_Size_flag","STRATIO_flag","Math_Trai
ning_flag","Certified_Teacher_flag","Prop_Girl_flag"), formula = "Gender + ESCS +
ICTHOME_binary * ICTSCH_binary + ICTQUAL_cat + ICTRES_cat + Class_Size +
STRATIO + Math_Training + Certified_Teacher + Prop_Girl + Ownership + ESCS_flag +
Class_Size_flag + STRATIO_flag + Math_Training_flag + Certified_Teacher_flag +
Prop_Girl_flag", pvlabel = paste0("PV",1:10,"MATH"),data = complete_data)
#Main Model 3
model_3 <-
  pisa.reg.pv_modified(x=c("Gender","ESCS","ICTHOME_binary","ICTSCH_binary",
                           "ICTQUAL_cat","ICTRES_cat","Class_Size","STRATIO","Math_Training","Certified_Teac
her","Prop_Girl","Ownership","ESCS_flag","Class_Size_flag","STRATIO_flag","Math_Trai
ning_flag","Certified_Teacher_flag","Prop_Girl_flag"), formula = "Gender + ESCS +
ICTHOME_binary * ICTRES_cat + ICTSCH_binary + ICTQUAL_cat + Class_Size +
STRATIO + Math_Training + Certified_Teacher + Prop_Girl + Ownership + ESCS_flag +

66

Class_Size_flag + STRATIO_flag + Math_Training_flag + Certified_Teacher_flag +
Prop_Girl_flag", pvlabel = paste0("PV",1:10,"MATH"),data = complete_data)
#Main Model 4
model_4 <-
  pisa.reg.pv_modified(x=c("Gender","ESCS","ICTHOME_binary","ICTSCH_binary",
                           "ICTQUAL_cat","ICTRES_cat","Class_Size","STRATIO","Math_Training","Certified_Teac
her","Prop_Girl","Ownership","ESCS_flag","Class_Size_flag","STRATIO_flag","Math_Trai
ning_flag","Certified_Teacher_flag","Prop_Girl_flag"), formula = "Gender + ESCS +
ICTSCH_binary * ICTQUAL_cat + ICTHOME_binary + ICTRES_cat + Class_Size +
STRATIO + Math_Training + Certified_Teacher + Prop_Girl + Ownership + ESCS_flag +
Class_Size_flag + STRATIO_flag + Math_Training_flag + Certified_Teacher_flag +
Prop_Girl_flag", pvlabel = paste0("PV",1:10,"MATH"),data = complete_data)
#Visualising Regression Result Table
#Model 1
model_1_df <- data.frame(round(model_1$reg,3))
names(model_1_df) <- c("Estimated Coefficient","Standard Error","T Value")
model_1_row_order <-
  c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,1,28)
rownames(model_1_df) <- c("Intercept","Gender (Male)","Gender (Missing)","Economic
Social and Cultural Status","ICT at Home (High)", "ICT at Home (Missing)","ICT at School
(High)","ICT at School (Missing)","ICT Quality (Medium)","ICT Quality (High)","ICT
Quality (Missing)","ICT Resource (Medium)","ICT Resource (High)","ICT Resource
(Missing)","Maths Class Size","Student-Teacher Ratio","Maths Professional
Training","Proportion of Certified Teacher","Proportion of Girls","Ownership
(Private)","Ownership (Missing)","ESCS Missing Flag","Maths Class Size Missing
Flag","Student-Teacher Ratio Missing Flag","Maths Professional Training Missing

Flag","Proportion of Certified Teacher Missing Flag","Proportion of Girls Missing Flag","R-
Squared")

model_1_df <- model_1_df[model_1_row_order,]
names(model_1) <- c("Estimated Coefficient","Standard Error","T Value")
#Model 2
model_2_df <- data.frame(round(model_2$reg,3))
names(model_2_df) <- c("Estimated Coefficient","Standard Error","T Value")
model_2_row_order <-
  c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,1,32)
rownames(model_2_df) <- c("Intercept","Gender (Male)","Gender (Missing)","Economic
Social and Cultural Status","ICT at Home (High)", "ICT at Home (Missing)","ICT at School
(High)","ICT at School (Missing)","ICT Quality (Medium)","ICT Quality (High)","ICT
Quality (Missing)","ICT Resource (Medium)","ICT Resource (High)","ICT Resource
(Missing)","Maths Class Size","Student-Teacher Ratio","Maths Professional
Training","Proportion of Certified Teacher","Proportion of Girls","Ownership
(Private)","Ownership (Missing)","ESCS Missing Flag","Maths Class Size Missing
Flag","Student-Teacher Ratio Missing Flag","Maths Professional Training Missing
Flag","Proportion of Certified Teacher Missing Flag","Proportion of Girls Missing
Flag","ICT at Home (High) : ICT at School (High)","ICT at Home (Missing) : ICT at School
(High)","ICT at Home (High) : ICT at School (Missing)","ICT at Home (Missing) : ICT at
School (Missing)","R-Squared")

67

model_2_df <- model_2_df[model_2_row_order,]
names(model_2) <- c("Estimated Coefficient","Standard Error","T Value")
#Model 3
model_3_df <- data.frame(round(model_3$reg,3))
names(model_3_df) <- c("Estimated Coefficient","Standard Error","T Value")
model_3_row_order <-
  c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,1,
    34)
rownames(model_3_df) <- c("Intercept","Gender (Male)","Gender (Missing)","Economic
Social and Cultural Status","ICT at Home (High)", "ICT at Home (Missing)","ICT Resource
(Medium)","ICT Resource (High)","ICT Resource (Missing)","ICT at School (High)","ICT
at School (Missing)","ICT Quality (Medium)","ICT Quality (High)","ICT Quality
(Missing)","Maths Class Size","Student-Teacher Ratio","Maths Professional
Training","Proportion of Certified Teacher","Proportion of Girls","Ownership
(Private)","Ownership (Missing)","ESCS Missing Flag","Maths Class Size Missing
Flag","Student-Teacher Ratio Missing Flag","Maths Professional Training Missing
Flag","Proportion of Certified Teacher Missing Flag","Proportion of Girls Missing
Flag","ICT at Home (High) : ICT Resource (Medium)","ICT at Home (Missing) : ICT
Resource (Medium)","ICT at Home (High) : ICT Resource (High)","ICT at Home (Missing) :
ICT Resource (High)","ICT at Home (High) : ICT Resource (Missing)","ICT at Home
(Missing) : ICT Resource (Missing)","R-Squared")
model_3_df <- model_3_df[model_3_row_order,]
names(model_3) <- c("Estimated Coefficient","Standard Error","T Value")
#Model 4
model_4_df <- data.frame(round(model_4$reg,3))
names(model_4_df) <- c("Estimated Coefficient","Standard Error","T Value")
model_4_row_order <-
  c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,1,
    34)
rownames(model_4_df) <- c("Intercept","Gender (Male)","Gender (Missing)","Economic
Social and Cultural Status","ICT at School (High)","ICT at School (Missing)","ICT Quality
(Medium)","ICT Quality (High)","ICT Quality (Missing)","ICT at Home (High)", "ICT at
Home (Missing)","ICT Resource (Medium)","ICT Resource (High)","ICT Resource
(Missing)","Maths Class Size","Student-Teacher Ratio","Maths Professional
Training","Proportion of Certified Teacher","Proportion of Girls","Ownership
(Private)","Ownership (Missing)","ESCS Missing Flag","Maths Class Size Missing
Flag","Student-Teacher Ratio Missing Flag","Maths Professional Training Missing
Flag","Proportion of Certified Teacher Missing Flag","Proportion of Girls Missing
Flag","ICT at School (High) : ICT Quality (Medium)","ICT at School (Missing) : ICT
Quality (Medium)","ICT at School (High) : ICT Quality (High)","ICT at School (Missing) :
ICT Quality (High)","ICT at School (High) : ICT Quality (Missing)","ICT at School
(Missing) : ICT Quality (Missing)","R-Squared")
model_4_df <- model_4_df[model_4_row_order,]
names(model_4) <- c("Estimated Coefficient","Standard Error","T Value")
##Test Regression Assumption##
#Multicollinearity - VIF

68

#Dependent variable doesn't matter for VIF
vif_model <- lm(PV3MATH ~ Gender + ESCS + ICTHOME_binary + ICTSCH_binary +
                  ICTRES_cat + ICTQUAL_cat + Class_Size + STRATIO + Math_Training +
                  Certified_Teacher + Prop_Girl + Ownership,data=complete_data)
vif(vif_model)
#Correlation Matrix with Modified Function
correlation_matrix <-
  pisa.rho_modified(variables=c("ESCS","Class_Size","STRATIO","Math_Training","Certifie
d_Teacher","Prop_Girl"), data=complete_data)
corrplot(correlation_matrix, type = "upper",order = "hclust", addCoef.col = 'black', tl.col =
           "black", tl.srt = 45)
##Robustness Check 1: Replace ESCS with Parental Education##
#Replace ESCS with Parental Education
robust_full_model <-
  pisa.reg.pv_modified(x=c("Gender","FatherEdu","MotherEdu","ICTHOME_binary","ICTSC
H_binary",
                           "ICTQUAL_cat","ICTRES_cat","Class_Size","STRATIO","Math_Training","Certified_Teac
her","Prop_Girl","Ownership","ESCS_flag","Class_Size_flag","STRATIO_flag","Math_Trai
ning_flag","Certified_Teacher_flag","Prop_Girl_flag"), formula = "Gender + FatherEdu +
MotherEdu + ICTHOME_binary + ICTSCH_binary + ICTQUAL_cat + ICTRES_cat +
Class_Size + STRATIO + Math_Training + Certified_Teacher + Prop_Girl + Ownership +
FatherEdu_flag + MotherEdu_flag + Class_Size_flag + STRATIO_flag +
Math_Training_flag + Certified_Teacher_flag + Prop_Girl_flag", pvlabel =
                         paste0("PV",1:10,"MATH"),data = complete_data)
#Robust Full
robust_full_df <- data.frame(round(robust_full_model$reg,3))
names(robust_full_df) <- c("Estimated Coefficient","Standard Error","T Value")
robust_row_order <-
  c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,1,30)
rownames(robust_full_df) <- c("Intercept","Gender (Male)","Gender (Missing)","Father's
Education","Mother's Education","ICT at Home (High)", "ICT at Home (Missing)","ICT at
School (High)","ICT at School (Missing)","ICT Quality (Medium)","ICT Quality
(High)","ICT Quality (Missing)","ICT Resource (Medium)","ICT Resource (High)","ICT
Resource (Missing)","Maths Class Size","Student-Teacher Ratio","Maths Professional
Training","Proportion of Certified Teacher","Proportion of Girls","Ownership
(Private)","Ownership (Missing)","Father's Education Missing Flag","Mother's Education
Missing Flag","Maths Class Size Missing Flag","Student-Teacher Ratio Missing
Flag","Maths Professional Training Missing Flag","Proportion of Certified Teacher Missing
Flag","Proportion of Girls Missing Flag","R-Squared")
robust_full_df <- robust_full_df[robust_row_order,]
names(robust_full_df) <- c("Estimated Coefficient","Standard Error","T Value")
##Robustness Check 2: Quantile Regression##
#Quantile Regression

69

library(quantreg)
#Quantile 1 - 25%
taus <- seq(0.1,0.9,by = 0.1)
quant_seq <- rq(PV1MATH ~ Gender + ESCS + ICTHOME_binary + ICTSCH_binary +
                  ICTQUAL_cat + ICTRES_cat + Class_Size + STRATIO + Math_Training +
                  Certified_Teacher + Prop_Girl + Ownership + ESCS_flag + Class_Size_flag +
                  STRATIO_flag + Math_Training_flag + Certified_Teacher_flag + Prop_Girl_flag, weights =
                  W_FSTUWT, data = complete_data,tau = taus)
quant_coef_df <- data.frame(round(quant_seq$coefficients,3))
names(quant_coef_df) <- c("10th Percentile","20th Percentile","30th Percentile","40th
Percentile","50th Percentile","60th Percentile","70th Percentile","80th Percentile","90th
Percentile")
rownames(quant_coef_df) <- c("Intercept","Gender (Male)","Gender (Missing)","Economic
Social and Cultural Status","ICT at Home (High)", "ICT at Home (Missing)","ICT at School
(High)","ICT at School (Missing)","ICT Quality (Medium)","ICT Quality (High)","ICT
Quality (Missing)","ICT Resource (Medium)","ICT Resource (High)","ICT Resource
(Missing)","Maths Class Size","Student-Teacher Ratio","Maths Professional
Training","Proportion of Certified Teacher","Proportion of Girls","Ownership
(Private)","Ownership (Missing)","ESCS Missing Flag","Maths Class Size Missing
Flag","Student-Teacher Ratio Missing Flag","Maths Professional Training Missing
Flag","Proportion of Certified Teacher Missing Flag","Proportion of Girls Missing Flag")
par(mfrow = c(2, 3))
#Quantile Regression Result of ICT at Home (High)
plot(x = seq(0.1,0.9,by = 0.1), y = quant_coef_df[5,],type = "o",xlab = "Percentile",ylab =
       "Coefficients",main = "ICT at Home (High)")
abline(h = 3.494, col = "red", lwd = 2)
#Quantile Regression Result of ICT at School (High)
plot(x = seq(0.1,0.9,by = 0.1), y = quant_coef_df[7,],type = "o",xlab = "Percentile",ylab =
       "Coefficients",main = "ICT at School (High)")
abline(h = -3, col = "red", lwd = 2)
#Quantile Regression Result of ICT Quality (Medium)
plot(x = seq(0.1,0.9,by = 0.1), y = quant_coef_df[9,],type = "o",xlab = "Percentile",ylab =
       "Coefficients",main = "ICT Quality (Medium)")
abline(h = 16.202, col = "red", lwd = 2)
#Quantile Regression Result of ICT Quality (High)
plot(x = seq(0.1,0.9,by = 0.1), y = quant_coef_df[10,],type = "o",xlab = "Percentile",ylab =
       "Coefficients",main = "ICT Quality (High)")
abline(h = 8, col = "red", lwd = 2)
#Quantile Regression Result of ICT Resource (Medium)
plot(x = seq(0.1,0.9,by = 0.1), y = quant_coef_df[12,],type = "o",xlab = "Percentile",ylab =
       "Coefficients",main = "ICT Resource (Medium)")
abline(h = 32.631, col = "red", lwd = 2)
#Quantile Regression Result of ICT Resource (High)
plot(x = seq(0.1,0.9,by = 0.1), y = quant_coef_df[13,],type = "o",xlab = "Percentile",ylab =
       "Coefficients",main = "ICT Resource (High)")
abline(h = 37.814, col = "red", lwd = 2)
##Tidy All Regression Results##

70

for (i in 1:nrow(model_1_df)){
  model_1_df$Model_1_coef[i] <- paste(toString(model_1_df$`Estimated
                                               Coefficient`[i]),"(",toString(model_1_df$`Standard Error`[i]),")")
  model_1_df$Model_1_coef[i] <- ifelse((1.645 >= model_1_df$`T Value`[i]) &
                                         (model_1_df$`T Value`[i] >= -1.645), model_1_df$Model_1_coef[i],
                                       paste(model_1_df$Model_1_coef[i],"*"))
  model_1_df$Model_1_coef[i] <- ifelse((-1.96 >= model_1_df$`T Value`[i]) |
                                         (model_1_df$`T Value`[i] >= 1.96), paste(model_1_df$Model_1_coef[i],"*"),
                                       model_1_df$Model_1_coef[i])
  model_1_df$Model_1_coef[i] <- ifelse((-2.576 >= model_1_df$`T Value`[i]) |
                                         (model_1_df$`T Value`[i] >= 2.576), paste(model_1_df$Model_1_coef[i],"*"),
                                       model_1_df$Model_1_coef[i])
}
for (i in 1:nrow(model_2_df)){
  model_2_df$Model_2_coef[i] <- paste(toString(model_2_df$`Estimated
                                               Coefficient`[i]),"(",toString(model_2_df$`Standard Error`[i]),")")
  model_2_df$Model_2_coef[i] <- ifelse((1.645 >= model_2_df$`T Value`[i]) &
                                         (model_2_df$`T Value`[i] >= -1.645), model_2_df$Model_2_coef[i],
                                       paste(model_2_df$Model_2_coef[i],"*"))
  model_2_df$Model_2_coef[i] <- ifelse((-1.96 >= model_2_df$`T Value`[i]) |
                                         (model_2_df$`T Value`[i] >= 1.96), paste(model_2_df$Model_2_coef[i],"*"),
                                       model_2_df$Model_2_coef[i])
  model_2_df$Model_2_coef[i] <- ifelse((-2.576 >= model_2_df$`T Value`[i]) |
                                         (model_2_df$`T Value`[i] >= 2.576), paste(model_2_df$Model_2_coef[i],"*"),
                                       model_2_df$Model_2_coef[i])
}
for (i in 1:nrow(model_3_df)){
  model_3_df$Model_3_coef[i] <- paste(toString(model_3_df$`Estimated
                                               Coefficient`[i]),"(",toString(model_3_df$`Standard Error`[i]),")")
  model_3_df$Model_3_coef[i] <- ifelse((1.645 >= model_3_df$`T Value`[i]) &
                                         (model_3_df$`T Value`[i] >= -1.645), model_3_df$Model_3_coef[i],
                                       paste(model_3_df$Model_3_coef[i],"*"))
  model_3_df$Model_3_coef[i] <- ifelse((-1.96 >= model_3_df$`T Value`[i]) |
                                         (model_3_df$`T Value`[i] >= 1.96), paste(model_3_df$Model_3_coef[i],"*"),
                                       model_3_df$Model_3_coef[i])
  model_3_df$Model_3_coef[i] <- ifelse((-2.576 >= model_3_df$`T Value`[i]) |
                                         (model_3_df$`T Value`[i] >= 2.576), paste(model_3_df$Model_3_coef[i],"*"),
                                       model_3_df$Model_3_coef[i])
}
for (i in 1:nrow(model_4_df)){
  model_4_df$Model_4_coef[i] <- paste(toString(model_4_df$`Estimated
                                               Coefficient`[i]),"(",toString(model_4_df$`Standard Error`[i]),")")
  model_4_df$Model_4_coef[i] <- ifelse((1.645 >= model_4_df$`T Value`[i]) &
                                         (model_4_df$`T Value`[i] >= -1.645), model_4_df$Model_4_coef[i],
                                       paste(model_4_df$Model_4_coef[i],"*"))
  model_4_df$Model_4_coef[i] <- ifelse((-1.96 >= model_4_df$`T Value`[i]) |
                                         (model_4_df$`T Value`[i] >= 1.96), paste(model_4_df$Model_4_coef[i],"*"),
                                       model_4_df$Model_4_coef[i])
  
  71
  
  model_4_df$Model_4_coef[i] <- ifelse((-2.576 >= model_4_df$`T Value`[i]) |
                                         (model_4_df$`T Value`[i] >= 2.576), paste(model_4_df$Model_4_coef[i],"*"),
                                       model_4_df$Model_4_coef[i])
}
for (i in 1:nrow(robust_full_df)){
  robust_full_df$Robust_coef[i] <- paste(toString(robust_full_df$`Estimated
                                                  Coefficient`[i]),"(",toString(robust_full_df$`Standard Error`[i]),")")
  robust_full_df$Robust_coef[i] <- ifelse((1.645 >= robust_full_df$`T Value`[i]) &
                                            (robust_full_df$`T Value`[i] >= -1.645), robust_full_df$Robust_coef[i],
                                          paste(robust_full_df$Robust_coef[i],"*"))
  robust_full_df$Robust_coef[i] <- ifelse((-1.96 >= robust_full_df$`T Value`[i]) |
                                            (robust_full_df$`T Value`[i] >= 1.96), paste(robust_full_df$Robust_coef[i],"*"),
                                          robust_full_df$Robust_coef[i])
  robust_full_df$Robust_coef[i] <- ifelse((-2.576 >= robust_full_df$`T Value`[i]) |
                                            (robust_full_df$`T Value`[i] >= 2.576), paste(robust_full_df$Robust_coef[i],"*"),
                                          robust_full_df$Robust_coef[i])
}
merge_12 <- merge(model_1_df, model_2_df, by = "row.names",all = TRUE)
merge_34 <- merge(model_3_df,model_4_df,by = "row.names",all = TRUE)
combined_result <- merge(merge_12,merge_34,by = "Row.names", all = TRUE)
combined_robust <- merge(model_1_df,robust_full_df,by = "row.names", all = TRUE)
result_coefs <-
  select(combined_result,Row.names,Model_1_coef,Model_2_coef,Model_3_coef,Model_4_c
         oef)
robust_coefs <- select(combined_robust,Row.names,Model_1_coef,Robust_coef)
names(result_coefs) <- c("Coefficients","Model 1","Model 2","Model 3","Model 4")
names(robust_coefs) <- c("Coefficients","Model 1","Robustness Check Model")
#Main Regression Coefficient DF
result_main <- result_coefs[c(1,3,5,17,29,28,26,25,6,9,8,19,18,31,42),]
rownames(result_main) <- 1:nrow(result_main)
#School Variables Coefficient DF
sch_var_df <- result_coefs[c(32,33,34,35,36,37,38,39,40,41,44,31,42),]
rownames(sch_var_df) <- 1:nrow(sch_var_df)
#Main Regression Coefficient DF with Missing Flag
result_main_flag <-
  result_coefs[c(1,2,3,4,5,11,17,21,29,28,30,26,25,27,6,7,12,13,9,8,10,15,14,16,19,18,20,23,22
                 ,24,31,42),]
rownames(result_main_flag) <- 1:nrow(result_main_flag)
#Robustness Regression Coefficient DF
result_robust <- robust_coefs[c(1,3,22,5,7,9,12,11,15,14,17,30),]
rownames(result_robust) <- 1:nrow(result_robust)
#Quantile Regression Coefficient DF
result_quantile <- quant_coef_df[c("Economic Social and Cultural Status","Gender
(Male)","ICT at Home (High)","ICT at School (High)","ICT Quality (Medium)","ICT
Quality (High)","ICT Resource (Medium)","ICT Resource (High)","Intercept"),]
kable(result_main, format = "latex", booktabs = TRUE) %>% kable_styling(latex_options =
                                                                          c("striped","scale_down"))

72

kable(sch_var_df, format = "latex", booktabs = TRUE) %>% kable_styling(latex_options =
                                                                         c("striped","scale_down"))
kable(result_main_flag, format = "latex", booktabs = TRUE) %>%
  kable_styling(latex_options = c("striped","scale_down"))
kable(result_robust, format = "latex", booktabs = TRUE) %>% kable_styling(latex_options =
                                                                            c("striped","scale_down"))
kable(result_quantile, format = "latex", booktabs = TRUE) %>% kable_styling(latex_options
                                                                            = c("striped","scale_down"))