   tempdir <- function() '/temp'
   unlockBinding("tempdir", baseenv())
   utils::assignInNamespace("tempdir", tempdir, ns="base", envir=baseenv())
   assign("tempdir", tempdir, baseenv())
   lockBinding("tempdir", baseenv())
   marcox(formula = Surv(time, cens) ~ sex, dat =init('D:/kidney.txt', div=2, dummy=dum(c('type'),list(c(0,1,2)))))
   Sys.setenv(CXXFLAGS = "-Wno-unused -Wno-templates")
   options(warn = -1)
   Rcpp::sourceCpp("R/kernal.cpp")



   kernal(dd_c=dd, g11_c=g11, t11_c=t11, tt1_c=tt1,
          t2_c=t2, betainit_h=betainit, x111_c=x111, xxx_c=xxx,
          c1_c=c1, n_c=n, id_c=id,
          Kn_ls_c=Kn_ls, kk_c=kk, Kn_c=Kn, K_c=K)

