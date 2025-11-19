

## Only diversity

m1 <- lmer(log_produt ~ SR + (1 | site), data = dadosmisto1, REML = FALSE)
summary(m1) # p-value = 0.00766 **
r.squaredGLMM(m1) # 0.5957724
AICc(m1) # 77.13385

m2 <- lmer(log_produt ~ SESPDab + (1 | site), data = dadosmisto1, REML = FALSE)
summary(m2) # p-value = 0.953
r.squaredGLMM(m2) # 0.5117967
AICc(m2) # 81.83127

m3 <- lmer(log_produt ~ SESMPDab + (1 | site), data = dadosmisto1, REML = FALSE)
summary(m3) # p-value = 0.514
r.squaredGLMM(m3) # 0.51117
AICc(m3) # 81.4064

m4 <- lmer(log_produt ~ SESMNTDab + (1 | site), data = dadosmisto1, REML = FALSE)
summary(m4) # p-value = 0.622    
r.squaredGLMM(m4) # 0.5288856
AICc(m4) # 81.60251

mpcps <- lmer(log_produt ~ pcps1ab + (1 | site), data = dadosmisto1, REML = FALSE)
summary(mpcps) # p-value = 0.0852 .    
r.squaredGLMM(mpcps) # 0.4888937
AICc(mpcps) # 78.97703

## Only environmental (climate / soil / slope)

m5 <- lmer(log_produt ~ silte + (1 | site), data = dadosmisto1, REML = FALSE)
summary(m5) 
anova(modelo_nulo, m5) # p-value = 0.08606 .
r.squaredGLMM(m5) # 0.4771867
AICc(m5) # 78.88817

m6 <- lmer(log_produt ~ c.n_soloid + (1 | site), data = dadosmisto1, REML = FALSE)
summary(m6) 
anova(modelo_nulo, m6) # p-value = 0.02618 *
r.squaredGLMM(m6) # 0.5271506
AICc(m6) # 76.89054

m7 <- lmer(log_produt ~ altitude + (1 | site), data = dadosmisto1, REML = FALSE)
summary(m7) 
anova(modelo_nulo, m7) # p-value = 0.0828 .
r.squaredGLMM(m7) # 0.5066781
AICc(m7) # 78.8257

m8 <- lmer(log_produt ~ c.n_soloid + silte + (1 | site), data = dadosmisto1, REML = FALSE)
summary(m8) 
anova(modelo_nulo, m8) # p-value = 0.02314 *
r.squaredGLMM(m8) # 0.4828996
AICc(m8) # 76.6528

m9 <- lmer(log_produt ~ PC1_clima + PC1nutri + (1 | site), data = dadosmisto1, REML = FALSE)
summary(m9) 
anova(modelo_nulo, m9) # p-value = 0.3342 
r.squaredGLMM(m9) # 0.5085717
AICc(m9) # 81.99316

## Diversity + Environmental

m10 <- lmer(log_produt ~ altitude + SR + pcps1ab + (1 | site), data = dadosmisto1, REML = FALSE)
summary(m10) 
anova(modelo_nulo, m10) # p-value = 0.0008819 ***
r.squaredGLMM(m10) # 0.5697163
AICc(m10) # 70.08434

m11 <- lmer(log_produt ~ c.n_soloid + SR + pcps1ab + (1 | site), data = dadosmisto1, REML = FALSE)
summary(m11) 
anova(modelo_nulo, m11) # p-value = 0.0003766 *** 
r.squaredGLMM(m11) # 0.5667528
AICc(m11) # 68.28968
