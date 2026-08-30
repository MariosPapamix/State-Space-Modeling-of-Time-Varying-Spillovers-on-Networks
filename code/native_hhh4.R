suppressMessages(library(surveillance))
# ---------- measles: export data + native hhh4 with rolling one-step ----------
data(measlesWeserEms)
Y <- observed(measlesWeserEms); NB <- neighbourhood(measlesWeserEms)
pop <- population(measlesWeserEms)[1,]
write.csv(Y, "../data/measles_counts.csv", row.names=FALSE)
write.csv(NB, "../data/measles_nborder.csv", row.names=FALSE)
write.csv(data.frame(pop=pop), "../data/measles_pop.csv", row.names=FALSE)
TRAIN <- 78
ctrl <- list(end = list(f = addSeason2formula(~ -1 + fe(1, unitSpecific=TRUE), S=1, period=52)),
             ar  = list(f = ~1),
             ne  = list(f = ~1, weights = NB == 1),
             family = "NegBin1", subset = 2:TRAIN)
fitM <- hhh4(measlesWeserEms, ctrl)
cat("[measles fit] logLik", logLik(fitM), " ar", coef(fitM)["ar.1"], " ne", coef(fitM)["ne.1"],
    " psi", coef(fitM)["overdisp"], "\n")
osaM <- oneStepAhead(fitM, tp=c(TRAIN, nrow(Y)-1), type="rolling", verbose=FALSE)
logsM <- scores(osaM, which="logs", individual=TRUE)
write.csv(logsM, "../data/measles_hhh4_logs.csv", row.names=FALSE)
write.csv(osaM$pred, "../data/measles_hhh4_pred.csv", row.names=FALSE)
write.csv(data.frame(psi=osaM$psi), "../data/measles_hhh4_psi.csv", row.names=FALSE)
co <- coef(fitM)
write.csv(data.frame(name=names(co), value=as.numeric(co)), "../data/measles_hhh4_coef.csv", row.names=FALSE)
cat("[measles osa] mean composite logs/week", mean(rowSums(logsM)), "\n")
# ---------- texas: same spec, native fit + rolling one-step ----------
Yt <- as.matrix(read.csv("../data/tx_counts.csv")); W <- as.matrix(read.csv("../data/tx_W.csv"))
stsT <- sts(observed=Yt, neighbourhood=W)
TR <- 72
ctrlT <- list(end = list(f = addSeason2formula(~ -1 + fe(1, unitSpecific=TRUE), S=1, period=52)),
              ar  = list(f = ~1),
              ne  = list(f = ~1, weights = W),
              family = "NegBin1", subset = 2:TR)
fitT <- hhh4(stsT, ctrlT)
cat("[texas fit] logLik", logLik(fitT), " ar", coef(fitT)["ar.1"], " ne", coef(fitT)["ne.1"],
    " psi", coef(fitT)["overdisp"], "\n")
co <- coef(fitT)
write.csv(data.frame(name=names(co), value=as.numeric(co)), "../data/tx_hhh4_coef.csv", row.names=FALSE)
osaT <- oneStepAhead(fitT, tp=c(TR, nrow(Yt)-1), type="rolling", verbose=FALSE)
logsT <- scores(osaT, which="logs", individual=TRUE)
write.csv(logsT, "../data/tx_hhh4_logs.csv", row.names=FALSE)
write.csv(osaT$pred, "../data/tx_hhh4_pred.csv", row.names=FALSE)
write.csv(data.frame(psi=osaT$psi), "../data/tx_hhh4_psi.csv", row.names=FALSE)
cat("[texas osa] mean composite logs/week", mean(rowSums(logsT)), "\n")
cat("DONE-R\n")
