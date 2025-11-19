### TESTE SEM CSA ###

library(writexl)
library(readxl)
library(lmerTest) #MODELO MISTO
library(MuMIn) #AICc

dadosmisto1 <- dadosmisto[-c(1:7), ]
modelo_nulo <- lmer(log_produt ~ (1 | site), data = dadosmisto1, REML = FALSE)

m1 <- lmer(log_produt ~ SR + silte + (1 | site), data = dadosmisto1, REML = FALSE)
summary(m1) 
anova(modelo_nulo, m1) # p-value = 0.002509 **
r.squaredGLMM(m1) # 0.4281354
AICc(m1) # 79.53395

m2 <- lmer(log_produt ~ SR + silte + pcps1ab + (1 | site), data = dadosmisto1, REML = FALSE)
summary(m2) 
anova(modelo_nulo, m2) # p-value = 0.0007806 ***
r.squaredGLMM(m2) # 0.4347053
AICc(m2) # 77.15064

m3 <- lmer(log_produt ~ SR + pcps1ab + (1 | site), data = dadosmisto1, REML = FALSE)
summary(m3) 
anova(modelo_nulo, m3) # p-value = 0.0007313 ***
r.squaredGLMM(m3) # 0.4698692
AICc(m3) # 77.06796


m4 <- lmer(log_produt ~ SESPDab + (1 | site), data = dadosmisto1, REML = FALSE)
summary(m4) 
anova(modelo_nulo, m4) # p-value = 0.621    
r.squaredGLMM(m4) # 0.40311
AICc(m4) # 88.91338

m5 <- lmer(log_produt ~ SR * pcps1ab + (1 | site), 
           data = dadosmisto1, REML = FALSE)
summary(m5)
AICc(m5) # 79.36204

m6 <- lmer(log_produt ~ riq_inicial + (1 | site), data = dadosmisto1, REML = FALSE)
summary(m6) 
anova(modelo_nulo, m6) # p-value = 0.88    
r.squaredGLMM(m6) # 0.4070914
AICc(m6) # 89.1361

m7 <- lmer(log_produt ~ ph + (1 | site), data = dadosmisto1, REML = FALSE)
summary(m7) # p-value = 0.0422 * 
r.squaredGLMM(m7) # 0.6564318
AICc(m7) # 85.92103

m8 <- lmer(log_produt ~ n_solo + (1 | site), data = dadosmisto1, REML = FALSE)
summary(m8) # p-value = 0.0449 * 
r.squaredGLMM(m8) # 0.3898037
AICc(m8) # 85.24197

m9 <- lmer(log_produt ~ c.n_soloid + (1 | site), data = dadosmisto1, REML = FALSE)
summary(m9) # p-value = 0.0064 ** 
r.squaredGLMM(m9) # 0.4032206
AICc(m9) # 81.93834


m10 <- lmer(log_produt ~ c.n_soloid + SR + ph + silte + pcps1ab + (1 | site), data = dadosmisto1, REML = FALSE)
summary(m10) 
anova(modelo_nulo, m10) # p-value = 0.0003433 ***    
r.squaredGLMM(m10) # 0.4529068
AICc(m10) # 76.10222


m11 <- lmer(log_produt ~ c.n_soloid + SR + silte + pcps1ab + (1 | site), data = dadosmisto1, REML = FALSE)
summary(m11) 
anova(modelo_nulo, m11) # p-value = 0.0001304 ***   
r.squaredGLMM(m11) # 0.4526873
AICc(m11) # 73.52156

m12 <- lmer(log_produt ~ c.n_soloid + SR + pcps1ab + (1 | site), data = dadosmisto1, REML = FALSE)
summary(m12) 
anova(modelo_nulo, m12) # p-value = 0.000116 ***  
r.squaredGLMM(m12) # 0.4417638
AICc(m12) # 73.14348
