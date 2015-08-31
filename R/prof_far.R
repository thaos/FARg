#' @importFrom parallel mclapply

#' @export
prof_ic=function(xp, t0, t1, y_fit, ci_p=0.95 ,to_plot=FALSE, ...){
  # 1) fit_gauss/gev/gpd -> init + overallmax likelihood
  # 2) get_far -> p0/p1= 1-FAR dans la fonction contraindre
  # 3) creer une fonction contraindre par type d'ajustement (gev, gpd, gaussien)
  # 4) pareil pour la fonction oprmiser profil
  # 5) on vérifie que l'optimisation sous contrainte du bon far donne les même résultats que l'optimisation sans contraine
  # 6) construire le profil comme précédemment
  # 7) normalement on ne pourrait faire qu'une seule fonction pour construire le profile des 3 types d'ajusteùents
  pnt0 <- set_pnt(t0, xp, y_fit$data)
  pnt1 <- set_pnt(t1, xp, y_fit$data)
  init <- y_fit$par
  overall_max <- y_fit$objective
  far <- get_far(y_fit, pnt0, pnt1, ...)
  p0p1 <- 1 - far[1]
  print("--- FAR -----")
  print(far)
  fit_check <- optimize_constr(y_fit, pnt0, pnt1, p0p1)
  if(!isTRUE(all.equal(fit_check$value, overall_max))){
    message("Slightly different results with constraint optimization")
    print(fit_check$value)
    print(overall_max)
    #             stop()
  }
  if(!isTRUE(all.equal(fit_check$par, init))){
    message("Slightly different results with constraint optimization")
    print(fit_check$par)
    print(init)
    #             stop()
  }
  aalpha <- qchisq(ci_p, 1)
  f_roots <- function(ratio){
    print(ratio)
    parmax  <-  -profil_optim(y_fit, pnt0, pnt1, R=ratio)
    parmax + aalpha/2 + overall_max 
  }
  print("BORNE INF-------------------------------------------")
  ic_inf <- binf_time(c(-0.1, p0p1), fun=f_roots, nbdiv=2 ,xmax=p0p1,fmax=aalpha/2)
  print("BORNE SUP-------------------------------------------")
  if(p0p1 <= 1)
    ic_sup <- bsup_time(c(p0p1,1.1),fun=f_roots,nbdiv=2,xmax=p0p1,fmax=aalpha/2)
  else
    ic_sup <- bsup_time(c(p0p1 ,p0p1*2), fun=f_roots, nbdiv=2, xmax=p0p1, fmax=aalpha/2)
  ratio_l <- c(ic_inf$x, ic_sup$x)
  parmax <- c(ic_inf$likel ,ic_sup$likel)
  if(abs(max(parmax)-aalpha/2)>0.0001){
    print("non equal mle")
    print(max(parmax))
    print(aalpha/2)
  }
  crit <- aalpha/2 - qchisq(0.999, 1)/2
  cond <- parmax > crit
  ratio_l <- ratio_l[cond]
  parmax <- parmax[cond]
  if(to_plot){
    smth <- spline(ratio_l[cond], parmax[cond], n = 500,method="natural")
    plot(ratio_l, parmax, type = "l", xlab = "", ylab = "")
    abline(h = 0, lty = 2, col = 2)
    lines(smth$x,smth$y, lty = 2, col = 2)
  }
  ci <- ratio_l[parmax > 0]
  stopifnot(length(ci)>1)
  out  <- c(1-max(ci), far[1], 1-min(ci),far[-1])
  names(out)[c(1,3)] <- c("IC_inf", "IC_sup")
  # out  <- c(min(ci),get_far()[1],max(ci),get_far()[-1])
  out
}

