rm(list = ls())

##################################################################
#Load library and auxiliar functions
library(CDVine)
library(energy)
library(mvtnorm)
#library(VineCopula)

##################################################################
source('simulate_multivariate.r')
#source('simulate_multivariate-CDVINEold.r')
source('runif_corel.R')


##################################################################
#SIMULATION
##################################################################
# parametres a faire evoluer si besoin
vect_var     <- 4 # nb_var<-seq(3,6, by=1) # nombre de variables
vect_corel   <- seq(0.5, 0.95, by = 0.05) #correlation
vect_nb_data <- c(seq(20, 50, by = 10), 100) #nombre de donnes
nb.iter <- 30 
nb.test <- 30 
vect_method <- c("indep", "indepPCA", "copula")
vect_dataset_form <- c("banana", "cigar")


#fonction de simulation

simulating_random_generation <- function(vect_var,
                                         vect_corel,
                                         vect_nb_data,
                                         nb.iter,
                                         nb.test,
                                         vect_method,
                                         vect_dataset_form)
{
 
  nb_block <- (length(vect_var) * length(vect_corel) * length(vect_nb_data) * 
              nb.iter * length(vect_dataset_form))
  l_block <- nb.test * length(vect_method) # block length
  
  nb_total_iteration <- nb_block * l_block

  multitest <- data.frame(var     = numeric(nb_total_iteration),
                          corel   = numeric(nb_total_iteration),
                          nb_data = numeric(nb_total_iteration),
                          form    = factor(character(nb_total_iteration),
                                           levels = vect_dataset_form),
                          method  = factor(character(nb_total_iteration),
                                           levels = vect_method),
                          p_value = numeric(nb_total_iteration))
  
  index_iteration <- 0
  for (nb_var in vect_var)
  { 
    for (corel in vect_corel)
    { 
      for (nb_data in vect_nb_data)
      { 
        for (form in vect_dataset_form)
          for (i in 1:nb.iter)
          {
            df <- runif_corel(nb_data, nb_var, corel ^ 2)
            #a changer si nb_var!=4
            if (form == "banana") {
              df[, 2] <-
                (2 * (df[, 2] - mean(df, 2))) ^ 2
            } #la variable 2 est mise en "banane"
            df[, 3] <- log(df[, 3])                  # var3 devient log
            df[, 4] <- exp((df[, 4] - min(df[, 4]))) # var4 devient exp
            df <- scale(df)

            multitest[index_iteration + 1:l_block, ] <-
            foreach(met = vect_method, .combine = rbind, .inorder = FALSE) %:%
              foreach (j = 1:nb.test, .combine = rbind, .inorder = FALSE) %dopar%
              {
                res <- simulatemulti(data = df, method = met)
                p_value <-
                  eqdist.etest(rbind(df, res), c(nrow(df), nrow(res)), R = 999)$p.value
                data.frame(var     = nb_var, 
                           corel   = corel, 
                           nb_data = nb_data, 
                           form    = form,
                           method  = met,
                           p_value = p_value)
              }
            index_iteration <- index_iteration + l_block
          }
      }
    }
  }
  simulating_random_generation <- multitest
  
}

##########################
#sauvegarde des simu!
##########################
library(doParallel)
library(foreach)

registerDoParallel(detectCores())
#registerDoParallel(8)

#T1 <- Sys.time()
system.time(resultat <-
  simulating_random_generation(vect_var,
                               vect_corel,
                               vect_nb_data,
                               nb.iter,
                               nb.test,
                               vect_method,
                               vect_dataset_form)
)
stopImplicitCluster()
#T2 <- Sys.time()
#Tdiff = difftime(T2, T1)

write.csv2(resultat , file = "~/Dropbox/Jairo Antoine/resultat_simu_complete-par.csv")
