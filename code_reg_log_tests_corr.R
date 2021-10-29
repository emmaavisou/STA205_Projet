# =======Régression Logistique ======= #

projet$succes <- ifelse(projet$MRI_LD_PreSurg == 0, 1, 0) # 1 si taille tumeur = 0 avant operation, 0 sinon

m1 <- glm(succes ~ age + prp + erp + MRI_LD_Baseline, data = projet, family = binomial)
logistic.display(m1)


# =======Tests bruts======= #

###Asso entre mesures (ou variations mesures) et erp###

var.test(projet$m12[projet$erp==0], projet$m12[projet$erp==1])
#variances égales

t.test(projet$m12[projet$erp==0], projet$m12[projet$erp==1], var.equal = T)
#pas d'asso stat signi



var.test(projet$m13[projet$erp==0], projet$m13[projet$erp==1])
#variances égales

t.test(projet$m13[projet$erp==0], projet$m13[projet$erp==1], var.equal = T)
#pas d'asso stat signi (p-value=0,07)



var.test(projet$m23[projet$erp==0], projet$m23[projet$erp==1])
#variances non égales (p-valeur=0,03)

t.test(projet$m23[projet$erp==0], projet$m23[projet$erp==1], var.equal = F)
#pas d'asso stat signi 





var.test(projet$MRI_LD_Baseline[projet$erp==0], projet$MRI_LD_Baseline[projet$erp==1])
#égalité

t.test(projet$MRI_LD_Baseline[projet$erp==0], projet$MRI_LD_Baseline[projet$erp==1], var.equal = T)
#pas d'asso




var.test(projet$MRI_LD_Int_Reg[projet$erp==0], projet$MRI_LD_Int_Reg[projet$erp==1])
#égalité

t.test(projet$MRI_LD_Int_Reg[projet$erp==0], projet$MRI_LD_Int_Reg[projet$erp==1], var.equal = T)
#pas d'asso




var.test(projet$MRI_LD_PreSurg[projet$erp==0], projet$MRI_LD_PreSurg[projet$erp==1])
#égalité

t.test(projet$MRI_LD_PreSurg[projet$erp==0], projet$MRI_LD_PreSurg[projet$erp==1], var.equal = T)
#pas d'asso


### Matrice de corrélation ###
library(corrplot)
Corr <- cor(sapply(projet,as.numeric), use="pairwise.complete.obs", method="spearman")
corrplot(Corr, method="square", type="upper", tl.col="black")

