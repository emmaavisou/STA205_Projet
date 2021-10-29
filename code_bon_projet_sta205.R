
# ======================================== #
# ============ PROJET STA 205 ============ #
# ======================================== #



#rm(list = ls())

# ======= Librairies et données ======= #

library(epiR)
library(epiDisplay)
library(survival)
library(cmprsk)

projet <- read.delim("projet1.txt")


# ======= Modification des types de variables ======= #

projet$MRI_LD_Baseline <- as.numeric(projet$MRI_LD_Baseline)
projet$MRI_LD_Int_Reg <- as.numeric(projet$MRI_LD_Int_Reg)
projet$MRI_LD_PreSurg <- as.numeric(projet$MRI_LD_PreSurg)
projet$id <- as.factor(projet$id)
projet$erp <- as.factor(projet$erp)
projet$prp <- as.factor(projet$prp)
projet$hrp <- as.factor(projet$hrp)
projet$bilat <- as.factor(projet$bilat)
projet$statut <- as.numeric(projet$statut)


# ======== Analyse descriptives des variables ==== #
# Variables quantitatives
# age
epi.descriptives(projet$age)
# delai
epi.descriptives(projet$delai)
# MRI_LD_Baseline
epi.descriptives(projet$MRI_LD_Baseline)
# MRI_LD_Int_Reg
epi.descriptives(projet$MRI_LD_Int_Reg)
# MRI_LD_PreSurg
epi.descriptives(projet$MRI_LD_PreSurg)


# Variables qualitatives
tab1(projet$erp)
tab1(projet$prp)
tab1(projet$hrp)
tab1(projet$bilat)
tab1(projet$statut)

# ======= Création de variables ======= #

# évolution de la taille de la tumeur entre avant chimiothérapie et entre les 2 cycles de traitement
projet$m12 <- ((projet$MRI_LD_Int_Reg - projet$MRI_LD_Baseline) / projet$MRI_LD_Baseline)
# évolution de la taille de la tumeur entre avant chimiothérapie et entre avant l'opération
projet$m13 <- ((projet$MRI_LD_PreSurg - projet$MRI_LD_Baseline) / projet$MRI_LD_Baseline)
# évolution de la taille de la tumeur entre les 2 cycles de traitement et entre avant l'opération
projet$m23 <- ((projet$MRI_LD_PreSurg - projet$MRI_LD_Int_Reg) / projet$MRI_LD_Int_Reg)

# Classes (augmentation ou diminution)
projet$cm12 <- ifelse(projet$m12 >= 0, 1, 0) # 1 si augmentation entre 1 et 2, 0 sinon
projet$cm13 <- ifelse(projet$m13 >= 0, 1, 0) # 1 si augmentation entre 1 et 3, 0 sinon
projet$cm23 <- ifelse(projet$m23 >= 0, 1, 0) # 1 si augmentation entre 2 et 3, 0 sinon

tab1(projet$cm12)
tab1(projet$cm13)
tab1(projet$cm23)

# Classes (1 si inférieur à -0.5 ; 0 si entre -0.5 et 0 ; 2 si plus de 0)
projet$cm12bis <- ifelse(projet$m12 <= 0 & projet$m12 >= -0.5, 0, ifelse(projet$m12 <= -0.5, 1, 2))
projet$cm13bis <- ifelse(projet$m13 <= 0 & projet$m13 >= -0.5, 0, ifelse(projet$m13 <= -0.5, 1, 2))
projet$cm23bis <- ifelse(projet$m23 <= 0 & projet$m23 >= -0.5, 0, ifelse(projet$m23 <= -0.5, 1, 2))

# Succès chimio (pour régression logistique)
projet$succes <- ifelse(projet$MRI_LD_PreSurg == 0, 1, 0) # 1 si taille tumeur = 0 avant opération, 0 sinon
tab1(projet$succes)
tab1(projet$chormones)
# Classes de succès
projet$csucces <- ifelse(projet$MRI_LD_PreSurg >= 45, 3, ifelse(projet$MRI_LD_PreSurg >= 22, 2, ifelse(projet$MRI_LD_PreSurg >= 4.75, 1, 0)))
projet$csucces <- as.factor(projet$csucces)

for (i in c(1:nrow(projet))){
  if(is.na(projet$m23[i]) == TRUE){
    projet$m23[i] <- 0
  }
}

# Classes de ERP/PRP/HRP
projet$chormones <- ifelse((projet$erp==1&projet$prp==0), 0, ifelse((projet$prp==1&projet$erp==1), 1, ifelse((projet$prp == 1& projet$erp ==0), 2, 3)))

# ======= Modèles de Cox ======= #

# Modèle de Cox avec les tailles de tumeur
cox1 <- coxph(Surv(delai,statut)~age+erp+prp+MRI_LD_Baseline+MRI_LD_Int_Reg+MRI_LD_PreSurg,data=projet)
summary(cox1)
prop1 <- cox.zph(cox1, transform="identity")
print(prop1)


# Modèle de Cox avec les évolutions de taille de tumeur
cox2 <- coxph(Surv(delai,statut)~age+erp+prp+m12+m13+m23,data=projet)
summary(cox2)
drop1(cox2, test="Chisq")
prop2 <- cox.zph(cox2, transform="identity")
print(prop2)

# Modèle de Cox avec les classes des évolutions de taille de tumeur
cox3 <- coxph(Surv(delai,statut) ~ age + erp + prp + cm12 + cm13 + cm23, data = projet)
summary(cox3)

# Modèle de Cox : il est BIEN ce modèle !
cox4 <- coxph(Surv(delai,statut) ~ age + erp + prp + MRI_LD_PreSurg + factor(cm12), data = projet)
summary(cox4)
# __ Vérification proportionnalité des risques
prop4 <- cox.zph(cox4, transform="identity")
print(prop4)

# Modèle TOP TOP TOP
cox5 <- coxph(Surv(delai, statut)~age +factor(chormones,c(1,0,2,3))+MRI_LD_PreSurg+factor(cm12), data = projet)
summary(cox5)

cox6 <- coxph(Surv(delai,statut)~age + MRI_LD_PreSurg+factor(cm12), data=projet)
summary(cox6)
anova(cox5)

# __ Vérification proportionnalité des risques
prop5 <- cox.zph(cox5, transform="identity")
print(prop5)

# Pour avoir la p-valeur globale des récepteurs hormonaux (variables chormones)
vrai <- -2*(cox5$loglik[2] - cox6$loglik[2])
vrai
1-pchisq(16.6047,df=3)


# ========= Beau Kaplan Meier =============#
install.packages("survminer")
library(survminer)
require("survival")

fits<-list()
# Kaplan- Meier de la variable hrp : significatif (à mettre dans résultats) --> mais dire qu'on ne la garde pas dans le modèle car elle est trop générale (non distinstion de ERP et PRP)
fit <- survfit(Surv(delai,statut) ~ hrp, data = projet)
fits[[1]] <- ggsurvplot(fit, 
                     pval = TRUE, conf.int = FALSE,
                     pval.size = 4, pval.coord = c(2000, .25),
                     risk.table = TRUE, 
                     risk.table.col = "strata", 
                     linetype = "solid", 
                     censor=FALSE,
                     palette = c("#B8255F","#FF9933"),
                     xlab="Temps (en jour)",
                     ylab="Probabilité de survie",
                     legend.title = "Statut récepteur hormonal HRP",
                     legend.labs = c("négatif", "positif"),
                     legend = c(.3,.3),
                     ggtheme = theme_classic2(base_size=12, base_family = "Arial"),
                     risk.table.title = "Nombre de personnes à risque",
                     risk.table.fontsize = 4,
                     risk.table.y.text = FALSE
)
fits[[1]]$table <- fits[[1]]$table + theme(plot.title = element_text(size = 12, hjust = 0, face = "bold")) + labs(y = NULL)
fits[[1]]$plot <- fits[[1]]$plot + theme(plot.title = element_text(hjust = .5)) + labs(title = "A")
fits[[1]]

# Kaplan- Meier de la variable ERP : significatif (à mettre dans résultats)
fit0 <- survfit(Surv(delai,statut) ~ erp, data = projet)
fits[[2]] <- ggsurvplot(fit0, 
                     pval = TRUE, conf.int = FALSE,
                     pval.size = 4, pval.coord = c(2000, .25),
                     risk.table = TRUE, 
                     risk.table.col = "strata", 
                     linetype = "solid", 
                     censor=FALSE,
                     palette = c("#B8255F","#FF9933"),
                     xlab="Temps (en jour)",
                     ylab="Probabilité de survie",
                     legend.title = "Statut récepteur hormonal ERP",
                     legend.labs = c("négatif", "positif"),
                     legend = c(.3,.3),
                     ggtheme = theme_classic2(base_size=12, base_family = "Arial"),
                     risk.table.title = "Nombre de personnes à risque",
                     risk.table.fontsize = 4,
                     risk.table.y.text = FALSE
)
fits[[2]]$table <- fits[[2]]$table + theme(plot.title = element_text(size = 12, hjust = 0, face = "bold")) + labs(y = NULL)
fits[[2]]$plot <- fits[[2]]$plot + theme(plot.title = element_text(hjust = .5)) + labs(title = "B")
fits[[2]]

arrange_ggsurvplots(fits, print=TRUE, ncol=2, nrow=1, risk.table.height = 0.25, surv.plot.height = 0.1)
# Kaplan- Meier de la variable PRP : non significatif (à mettre dans en annexe) --> même si on le met en annexe on explique que comme HRP est significatifs on crée la variable chormones et on fait le fit7
fit1 <- survfit(Surv(delai,statut) ~ prp, data = projet)
ggsurv <- ggsurvplot(fit1, 
                     pval = TRUE, conf.int = FALSE,
                     pval.size = 4, pval.coord = c(2000, .25),
                     risk.table = TRUE, 
                     risk.table.col = "strata", 
                     linetype = "solid", 
                     censor=FALSE,
                     palette = c("#B8255F","#FF9933"),
                     xlab="Temps (en jour)",
                     ylab="Probabilité de survie",
                     legend.title = "Statut récepteur hormonal PRP",
                     legend.labs = c("négatif", "positif"),
                     legend = c(.4,.3),
                     ggtheme = theme_classic2(base_size=12, base_family = "Arial"),
                     risk.table.title = "Nombre de personnes à risque",
                     risk.table.fontsize = 4,
                     risk.table.y.text = FALSE
)
ggsurv$table <- ggsurv$table + theme(plot.title = element_text(size = 12, hjust = 0, face = "bold")) + labs(y = NULL)
ggsurv$plot <- ggsurv$plot + theme(plot.title = element_text(hjust = .5))
ggsurv

# Kaplan- Meier de la variation de la taille de la tumeur durant le traitement (entre baseline et avant opération) : significatif (à mettre dans résultats) --> On l'enlève dans le modèle de cox car une fois ajustée elle n'est plus significative, l'efficacité du traitement se voit durant la première phase uniquement (voir fit3 et fit4)

fits2<-list()
fit2 <- survfit(Surv(delai,statut) ~ cm13, data = projet)
fits2[[1]] <- ggsurvplot(fit2, 
                     pval = TRUE, conf.int = FALSE,
                     pval.size = 4, pval.coord = c(2000, .25),
                     risk.table = TRUE, 
                     risk.table.col = "strata", 
                     linetype = "solid", 
                     censor=FALSE,
                     palette = c("#72B000","#009999"),
                     xlab="Temps (en jour)",
                     ylab="Probabilité de survie",
                     legend.title = "Variation taille de la tumeur",
                     legend.labs = c("diminution", "augmentation"),
                     legend = c(.3,.3),
                     ggtheme = theme_classic2(base_size=10, base_family = "Arial"),
                     risk.table.title = "Nombre de personnes à risque",
                     risk.table.fontsize = 4,
                     risk.table.y.text = FALSE
)
fits2[[1]]$table <- fits2[[1]]$table + theme(plot.title = element_text(size = 10, hjust = 0, face = "bold")) + labs(y = NULL)
fits2[[1]]$plot <- fits2[[1]]$plot + theme(plot.title = element_text(hjust = .5))+ labs(title = "A")
fits2[[1]]

# Kaplan- Meier de la variation de la taille de la tumeur durant première phase de traitement (entre baseline et entre deux cycle de traitements) : significatif (à mettre dans résultats)
fit3 <- survfit(Surv(delai,statut) ~ cm23, data = projet)
fits2[[3]] <- ggsurvplot(fit3, 
                     pval = TRUE, conf.int = FALSE,
                     pval.size = 4, pval.coord = c(2000, .25),
                     risk.table = TRUE, 
                     risk.table.col = "strata", 
                     linetype = "solid", 
                     censor=FALSE,
                     palette = c("#72B000","#009999"),
                     xlab="Temps (en jour)",
                     ylab="Probabilité de survie",
                     legend.title = "Variation taille de la tumeur",
                     legend.labs = c("diminution", "augmentation"),
                     legend = c(.3,.3),
                     ggtheme = theme_classic2(base_size=10, base_family = "Arial"),
                     risk.table.title = "Nombre de personnes à risque",
                     risk.table.fontsize = 4,
                     risk.table.y.text = FALSE
)
fits2[[3]]$table <- fits2[[3]]$table + theme(plot.title = element_text(size = 10, hjust = 0, face = "bold")) + labs(y = NULL)
fits2[[3]]$plot <- fits2[[3]]$plot + theme(plot.title = element_text(hjust = .5))+ labs(title = "C")
fits2[[3]]

# Kaplan- Meier de la variation de la taille de la tumeur durant première phase de traitement (entre baseline et entre deux cycle de traitements) : significatif (à mettre dans résultats)
fit4 <- survfit(Surv(delai,statut) ~ cm12, data = projet)
fits2[[2]]<- ggsurvplot(fit4, 
                     pval = TRUE, conf.int = FALSE,
                     pval.size = 4, pval.coord = c(1900, .25),
                     risk.table = TRUE, 
                     risk.table.col = "strata", 
                     linetype = "solid", 
                     censor=FALSE,
                     palette = c("#72B000","#009999"),
                     xlab="Temps (en jour)",
                     ylab="Probabilité de survie",
                     legend.title = "Variation taille de la tumeur",
                     legend.labs = c("diminution", "augmentation"),
                     legend = c(.3,.3),
                     ggtheme = theme_classic2(base_size=10, base_family = "Arial"),
                     risk.table.title = "Nombre de personnes à risque",
                     risk.table.fontsize = 4,
                     risk.table.y.text = FALSE
)
fits2[[2]]$table <- fits2[[2]]$table + theme(plot.title = element_text(size = 10, hjust = 0, face = "bold")) + labs(y = NULL)
fits2[[2]]$plot <-fits2[[2]]$plot + theme(plot.title = element_text(hjust = .5))+ labs(title = "B")
fits2[[2]]

arrange_ggsurvplots(fits2, print=TRUE, ncol=3, nrow=1, risk.table.height = 0.25, surv.plot.height = 0.1)


# Kaplan Meier avec la variable csucess : non significatif (à mettre en annexe)
fit6 <- survfit(Surv(delai,statut) ~ succes, data = projet)
ggsurv <- ggsurvplot(fit6, 
           pval = TRUE, conf.int = FALSE,
           pval.size = 4, pval.coord = c(2000, .25),
           risk.table = TRUE, 
           risk.table.col = "strata", 
           linetype = "solid", 
           censor=FALSE,
           palette = c("#009999","#72B000"),
           xlab="Temps (en jour)",
           ylab="Probabilité de survie",
           legend.title = "Tumeur avant operation (en mm)",
           legend.labs = c("égale à 0", "supérieur à 0"),
           legend = c(.3,.3),
           ggtheme = theme_classic2(base_size=12, base_family = "Arial"),
           risk.table.title = "Nombre de personnes à risque",
           risk.table.fontsize = 4,
           risk.table.y.text = FALSE
        )
ggsurv$table <- ggsurv$table + theme(plot.title = element_text(size = 12, hjust = 0, face = "bold")) + labs(y = NULL)
ggsurv$plot <- ggsurv$plot + theme(plot.title = element_text(hjust = .5))
ggsurv

# Kaplan Meier avec la variable chormones : significatif (à mettre dans résultats)
fit7 <- survfit(Surv(delai,statut) ~ chormones, data = projet)
ggsurv <- ggsurvplot(fit7, 
                     pval = TRUE, conf.int = FALSE,
                     pval.size = 4, pval.coord = c(2000, .25),
                     risk.table = TRUE, 
                     risk.table.col = "strata", 
                     linetype = "solid", 
                     censor=FALSE,
                     palette = c("#B8255F","#009999","#FF9933","#72B000"),
                     xlab="Temps (en jour)",
                     ylab="Probabilité de survie",
                     legend.title = "Récepteurs hormonaux",
                     legend.labs = c("ERP positif", "ERP et PRP positifs", "PRP positif", "ERP et PRP négatifs"),
                     legend = c(.2,.2),
                     ggtheme = theme_classic2(base_size=12, base_family = "Arial"),
                     risk.table.title = "Nombre de personnes à risque",
                     risk.table.fontsize = 4,
                     risk.table.y.text = FALSE
                     )
ggsurv$table <- ggsurv$table + theme(plot.title = element_text(size = 12, hjust = 0, face = "bold")) + labs(y = NULL)
ggsurv$plot <- ggsurv$plot + theme(plot.title = element_text(hjust = .5))
ggsurv
