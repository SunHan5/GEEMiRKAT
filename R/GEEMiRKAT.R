#' @export
GEEMiRKAT <- function (y, id, covs = NULL, Ks, model, n.perm = 5000) {
  options(warn = -1)
  if (is.null(covs)) {
    fit.GEE.AR <- MGEE(y ~ 1, id = id, family = model, corstr = "AR-1")
    f.y.GEE.AR <- fitted.values(fit.GEE.AR)
    r.GEE.AR <- y - f.y.GEE.AR
    fit.GEE.EX <- MGEE(y ~ 1, id = id, family = model, corstr = "exchangeable")
    f.y.GEE.EX <- fitted.values(fit.GEE.EX)
    r.GEE.EX <- y - f.y.GEE.EX
    fit.GEE.IN <- MGEE(y ~ 1, id = id, family = model, corstr = "independence")
    f.y.GEE.IN <- fitted.values(fit.GEE.IN)
    r.GEE.IN <- y - f.y.GEE.IN
  }
  else {
    fit.GEE.AR <- MGEE(y ~ ., id = id, family = model, corstr = "AR-1", as.data.frame(covs))
    f.y.GEE.AR <- fitted.values(fit.GEE.AR)
    r.GEE.AR <- y - f.y.GEE.AR
    fit.GEE.EX <- MGEE(y ~ ., id = id, family = model, corstr = "exchangeable", as.data.frame(covs))
    f.y.GEE.EX <- fitted.values(fit.GEE.EX)
    r.GEE.EX <- y - f.y.GEE.EX
    fit.GEE.IN <- MGEE(y ~ ., id = id, family = model, corstr = "independence", as.data.frame(covs))
    f.y.GEE.IN <- fitted.values(fit.GEE.IN)
    r.GEE.IN <- y - f.y.GEE.IN
  }
  
  if (model == "gaussian") {
    mu.GEE.AR <- fit.GEE.AR$linear.predictors
    mu.GEE.EX <- fit.GEE.EX$linear.predictors
    mu.GEE.IN <- fit.GEE.IN$linear.predictors
    inv.de.mu.GEE.AR <- inv.de.mu.GEE.EX <- inv.de.mu.GEE.IN <- diag(length(y))
    e.var.GEE.AR <- diag(rep(var(mu.GEE.AR), length(y)))
    e.var.GEE.EX <- diag(rep(var(mu.GEE.EX), length(y)))
    e.var.GEE.IN <- diag(rep(var(mu.GEE.IN), length(y)))
  }
  if (model == "binomial") {
    mu.GEE.AR <- exp(fit.GEE.AR$linear.predictors)/(1 + exp(fit.GEE.AR$linear.predictors))
    mu.GEE.EX <- exp(fit.GEE.EX$linear.predictors)/(1 + exp(fit.GEE.EX$linear.predictors))
    mu.GEE.IN <- exp(fit.GEE.IN$linear.predictors)/(1 + exp(fit.GEE.IN$linear.predictors))
    inv.de.mu.GEE.AR <- diag(as.numeric(mu.GEE.AR * (1 - mu.GEE.AR)))
    inv.de.mu.GEE.EX <- diag(as.numeric(mu.GEE.EX * (1 - mu.GEE.EX)))
    inv.de.mu.GEE.IN <- diag(as.numeric(mu.GEE.IN * (1 - mu.GEE.IN)))
    e.var.GEE.AR <- e.var.GEE.EX <- e.var.GEE.IN <- diag(rep(1, length(y)))
  }
  if (model == "poisson") {
    mu.GEE.AR <- exp(fit.GEE.AR$linear.predictors)
    mu.GEE.EX <- exp(fit.GEE.EX$linear.predictors)
    mu.GEE.IN <- exp(fit.GEE.IN$linear.predictors)
    inv.de.mu.GEE.AR <- diag(as.numeric(mu.GEE.AR))
    inv.de.mu.GEE.EX <- diag(as.numeric(mu.GEE.EX))
    inv.de.mu.GEE.IN <- diag(as.numeric(mu.GEE.IN))
    e.var.GEE.AR <- e.var.GEE.EX <- e.var.GEE.IN <- diag(rep(1, length(y)))
  }
  
  y.star.GEE.AR <- as.numeric(fit.GEE.AR$linear.predictors + ginv(inv.de.mu.GEE.AR) %*% (y - mu.GEE.AR))
  y.star.GEE.EX <- as.numeric(fit.GEE.EX$linear.predictors + ginv(inv.de.mu.GEE.EX) %*% (y - mu.GEE.EX))
  y.star.GEE.IN <- as.numeric(fit.GEE.IN$linear.predictors + ginv(inv.de.mu.GEE.IN) %*% (y - mu.GEE.IN))
  r.GEE.AR <- y.star.GEE.AR - fit.GEE.AR$linear.predictors
  r.GEE.EX <- y.star.GEE.EX - fit.GEE.EX$linear.predictors
  r.GEE.IN <- y.star.GEE.IN - fit.GEE.IN$linear.predictors
  
  v.cov.GEE.AR <- e.var.GEE.AR
  v.cov.GEE.EX <- e.var.GEE.EX
  v.cov.GEE.IN <- e.var.GEE.IN
  inv.vcov.GEE.AR <- ginv(v.cov.GEE.AR)
  inv.vcov.GEE.EX <- ginv(v.cov.GEE.EX)
  inv.vcov.GEE.IN <- ginv(v.cov.GEE.IN)
  
  r.s.GEE.AR <- r.s.GEE.EX <- r.s.GEE.IN <- list()
  for (j in 1:n.perm) {
    ind_r.s <- shuffle(length(y))
    r.s.GEE.AR[[j]] <- r.GEE.AR[ind_r.s]
    r.s.GEE.EX[[j]] <- r.GEE.EX[ind_r.s]
    r.s.GEE.IN[[j]] <- r.GEE.IN[ind_r.s]
  }
  
  Qs.GEE.AR <- Qs.GEE.EX <- Qs.GEE.IN <- rep(NA, length(Ks))
  for (j in 1:length(Ks)) {
    Qs.GEE.AR[j] <- (t(r.GEE.AR) %*% inv.vcov.GEE.AR %*% Ks[[j]] %*% inv.vcov.GEE.AR %*% r.GEE.AR)/(t(r.GEE.AR) %*% inv.vcov.GEE.AR %*% r.GEE.AR)
    Qs.GEE.EX[j] <- (t(r.GEE.EX) %*% inv.vcov.GEE.EX %*% Ks[[j]] %*% inv.vcov.GEE.EX %*% r.GEE.EX)/(t(r.GEE.EX) %*% inv.vcov.GEE.EX %*% r.GEE.EX)
    Qs.GEE.IN[j] <- (t(r.GEE.IN) %*% inv.vcov.GEE.IN %*% Ks[[j]] %*% inv.vcov.GEE.IN %*% r.GEE.IN)/(t(r.GEE.IN) %*% inv.vcov.GEE.IN %*% r.GEE.IN)
  }
  
  Q0s.GEE.AR <- Q0s.GEE.EX <- Q0s.GEE.IN <- list()
  for (j in 1:length(Ks)) {
    Q0s.inv.GEE.AR <- Q0s.inv.GEE.EX <- Q0s.inv.GEE.IN <- rep(NA, n.perm)
    for (k in 1:n.perm) {
      Q0s.inv.GEE.AR[k] <- (t(r.s.GEE.AR[[k]]) %*% inv.vcov.GEE.AR %*% Ks[[j]] %*% inv.vcov.GEE.AR %*% r.s.GEE.AR[[k]])/(t(r.s.GEE.AR[[k]]) %*% inv.vcov.GEE.AR %*% r.s.GEE.AR[[k]])
      Q0s.inv.GEE.EX[k] <- (t(r.s.GEE.EX[[k]]) %*% inv.vcov.GEE.EX %*% Ks[[j]] %*% inv.vcov.GEE.EX %*% r.s.GEE.EX[[k]])/(t(r.s.GEE.EX[[k]]) %*% inv.vcov.GEE.EX %*% r.s.GEE.EX[[k]])
      Q0s.inv.GEE.IN[k] <- (t(r.s.GEE.IN[[k]]) %*% inv.vcov.GEE.IN %*% Ks[[j]] %*% inv.vcov.GEE.IN %*% r.s.GEE.IN[[k]])/(t(r.s.GEE.IN[[k]]) %*% inv.vcov.GEE.IN %*% r.s.GEE.IN[[k]])
    }
    Q0s.GEE.AR[[j]] <- Q0s.inv.GEE.AR
    Q0s.GEE.EX[[j]] <- Q0s.inv.GEE.EX
    Q0s.GEE.IN[[j]] <- Q0s.inv.GEE.IN
  }
  
  pvs.GEE.AR <- pvs.GEE.EX <- pvs.GEE.IN <- rep(NA, length(Ks))
  for (j in 1:length(Ks)) {
    pvs.GEE.AR[j] <- (length(which(Q0s.GEE.AR[[j]] > Qs.GEE.AR[[j]])) + 1)/(n.perm + 1)
    pvs.GEE.EX[j] <- (length(which(Q0s.GEE.EX[[j]] > Qs.GEE.EX[[j]])) + 1)/(n.perm + 1)
    pvs.GEE.IN[j] <- (length(which(Q0s.GEE.IN[[j]] > Qs.GEE.IN[[j]])) + 1)/(n.perm + 1)
  }

  GEEMiRKAT.pvs <- rbind(pvs.GEE.AR, pvs.GEE.EX, pvs.GEE.IN)
  pv.opt.GEE <- apply(GEEMiRKAT.pvs, 1, ACAT)
  pv.opt <- ACAT(c(pvs.GEE.AR, pvs.GEE.EX, pvs.GEE.IN))
  p.aGEEMiRKAT <- c(pv.opt.GEE, pv.opt)
  colnames(GEEMiRKAT.pvs) <- names(Ks)
  rownames(GEEMiRKAT.pvs) <- c("GEE.AR", "GEE.EX", "GEE.IN")
  names(p.aGEEMiRKAT) <- c("GEEMiRKAT.AR", "GEEMiRKAT.EX", "GEEMiRKAT.IN", "aGEEMiRKAT")
  aGEEMiRKAT.output <- list(Ind.pvs = GEEMiRKAT.pvs, GEEMiRKAT.pvs = p.aGEEMiRKAT)
  return(aGEEMiRKAT.output = aGEEMiRKAT.output)
}
