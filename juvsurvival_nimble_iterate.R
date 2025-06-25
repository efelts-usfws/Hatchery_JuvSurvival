## Work out code to run a year of our data for juvenile
# survival through the hydrosystem

library(tidyverse)
library(nimble)
library(nimbleEcology)
library(MCMCvis)
library(arrow)
library(tictoc)
library(furrr)

# read in data from 2025 releases

ch.dat25 <- read_rds("data/travel")|> 
  mutate(capture_history=str_c("1",LGR_ch,LGS_ch,LOMO_ch,ICH_ch,
                               MCN_ch,JD_ch,BONN_ch,
                               TWX_ch,sep=""))


# use NIMBLE code to build the model

hmm.phitpt <- nimbleCode({
  # parameters
  delta[1] <- 1          # Pr(alive t = 1) = 1
  delta[2] <- 0          # Pr(dead t = 1) = 0
  for (t in 1:(T-1)){
    phi[t] ~ dunif(0, 1) # prior survival
    gamma[1,1,t] <- phi[t]      # Pr(alive t -> alive t+1)
    gamma[1,2,t] <- 1 - phi[t]  # Pr(alive t -> dead t+1)
    gamma[2,1,t] <- 0        # Pr(dead t -> alive t+1)
    gamma[2,2,t] <- 1        # Pr(dead t -> dead t+1)
    p[t] ~ dunif(0, 1) # prior detection
    omega[1,1,t] <- 1 - p[t]    # Pr(alive t -> non-detected t)
    omega[1,2,t] <- p[t]        # Pr(alive t -> detected t)
    omega[2,1,t] <- 1        # Pr(dead t -> non-detected t)
    omega[2,2,t] <- 0        # Pr(dead t -> detected t)
  }
  # likelihood
  for (i in 1:N){
    z[i,first[i]] ~ dcat(delta[1:2])
    for (j in (first[i]+1):T){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:2, j-1])
      y[i,j] ~ dcat(omega[z[i,j], 1:2, j-1])
    }
  }
})

# specify the parameters to monitor (phi and p)

parameters.to.save <- c("phi", "p")

# do 10,000 iterations with a 1000 burn-in period
# across 3 chains

n.iter <- 15000
n.burnin <- 5000
n.chains <- 3

# write a function to filter capture history info
# by release year and release site and then put those
# into the format that starts the process to go into nimble

nimble_dat.f <- function(hatchery.name,
                         release.year,release.group,
                         spp,
                         input.dat) {
  
  dat1 <- input.dat %>% 
    filter(hatchery==hatchery.name,
           release_year==release.year,
           release_group==release.group,
           species==spp,)
  
  dat2 <- dat1 %>% 
    ungroup() %>% 
    mutate(release_ch=1,LGR_ch=as.numeric(LGR_ch),
           LGS_ch=as.numeric(LGS_ch),LOMO_ch=as.numeric(LOMO_ch),
           ICH_ch=as.numeric(ICH_ch),MCN_ch=as.numeric(MCN_ch),
           JD_ch=as.numeric(JD_ch),BONN_ch=as.numeric(BONN_ch),
           TWX_ch=as.numeric(TWX_ch)) %>% 
    select(release_ch,LGR_ch,LGS_ch,LOMO_ch,ICH_ch,
           MCN_ch,JD_ch,BONN_ch,TWX_ch) %>% 
    as.matrix()
  
}



# define a function to be applied across groups

run_nimble_group.f <- function(hatchery.name,
                               release.year,release.group,
                               spp,
                               input.dat){
  
  
  dat <- nimble_dat.f(hatchery.name,release.year,release.group,
                      spp,input.dat)
  
  first <- apply(dat,1,function(x) min(which(x != 0)))
  
  constants <- list(N = nrow(dat),T=ncol(dat),first=first)
  
  data_list <- list(y=dat+1)
  
  zinits <- dat+1
  
  zinits[zinits==2] <- 1
  
  inits <- list(phi=runif(constants$T-1),
                p=runif(constants$T-1),
                z=zinits)
  
  model <- nimbleModel(code=hmm.phitpt,
                       constants = constants,
                       data=data_list,
                       inits=inits)
  
  c_model <- compileNimble(model)
  
  conf <- configureMCMC(model,monitors=parameters.to.save)
  
  mcmc <- buildMCMC(conf)
  
  c_mcmc <- compileNimble(mcmc,project=c_model)
  
  samples <- runMCMC(c_mcmc,
                     niter=n.iter,
                     nburnin=n.burnin,
                     nchains=n.chains,
                     samplesAsCodaMCMC = TRUE,
                     summary=FALSE)
  
  return(samples)
  
}

# for now try that function just on the clear creek chinook 

test_group <- ch.dat25 |> 
  distinct(hatchery,species,release_group,release_year) |> 
  slice(1:2)



plan(multisession,workers=2)


test_run <- future_pmap(test_group,
                        ~run_nimble_group.f(..1,..2,..3,..4,ch.dat25),
                        .progress = TRUE)

# write a function to extract paremeters of interest to export

chains_extract.f <- function(mcmc.object){
  
  phi1.1 <- mcmc.object$chain1[,"phi[1]"] %>% 
    as_tibble() %>% 
    mutate(chain="1",
           descriptor="release_lgr_phi",
           iteration=row_number())
  
  phi1.2 <-  mcmc.object$chain2[,"phi[1]"] %>% as_tibble() %>% 
    mutate(chain="2",
           descriptor="release_lgr_phi",
           iteration=row_number())
  
  phi2.1 <- mcmc.object$chain1[,"phi[2]"] %>% 
    as_tibble() %>% 
    mutate(chain="1",
           descriptor="lgr_lgs_phi",
           iteration=row_number())
  
  phi2.2 <- mcmc.object$chain2[,"phi[2]"] %>% 
    as_tibble() %>% 
    mutate(chain="2",
           descriptor="lgr_lgs_phi",
           iteration=row_number())
  
  phi3.1 <- mcmc.object$chain1[,"phi[3]"] %>% 
    as_tibble() %>% 
    mutate(chain="1",
           descriptor="lgs_lomo_phi",
           iteration=row_number())
  
  phi3.2 <- mcmc.object$chain2[,"phi[3]"] %>% 
    as_tibble() %>% 
    mutate(chain="2",
           descriptor="lgr_lomo_phi",
           iteration=row_number())
  
  phi4.1 <- mcmc.object$chain1[,"phi[4]"] %>% 
    as_tibble() %>% 
    mutate(chain="1",
           descriptor="lomo_ich_phi",
           iteration=row_number())
  
  phi4.2 <- mcmc.object$chain2[,"phi[4]"] %>% 
    as_tibble() %>% 
    mutate(chain="2",
           descriptor="lomo_ich_phi",
           iteration=row_number())
  
  phi5.1 <- mcmc.object$chain1[,"phi[5]"] %>% 
    as_tibble() %>% 
    mutate(chain="1",
           descriptor="ich_mcn_phi",
           iteration=row_number())
  
  phi5.2 <- mcmc.object$chain2[,"phi[5]"] %>% 
    as_tibble() %>% 
    mutate(chain="2",
           descriptor="ich_mcn_phi",
           iteration=row_number())
  
  phi6.1 <- mcmc.object$chain1[,"phi[6]"] %>% 
    as_tibble() %>% 
    mutate(chain="1",
           descriptor="mcn_jda_phi",
           iteration=row_number())
  
  phi6.2 <- mcmc.object$chain2[,"phi[6]"] %>% 
    as_tibble() %>% 
    mutate(chain="2",
           descriptor="mcn_jda_phi",
           iteration=row_number())
  
  phi7.1 <- mcmc.object$chain1[,"phi[7]"] %>% 
    as_tibble() %>% 
    mutate(chain="1",
           descriptor="jda_bon_phi",
           iteration=row_number())
  
  phi7.2 <- mcmc.object$chain2[,"phi[7]"] %>% 
    as_tibble() %>% 
    mutate(chain="2",
           descriptor="jda_bon_phi",
           iteration=row_number())
  
  phi8.1 <- mcmc.object$chain1[,"phi[8]"] %>% 
    as_tibble() %>% 
    mutate(chain="1",
           descriptor="bon_twx_phi",
           iteration=row_number())
  
  phi8.2 <- mcmc.object$chain2[,"phi[8]"] %>% 
    as_tibble() %>% 
    mutate(chain="2",
           descriptor="bon_twx_phi",
           iteration=row_number())
  
  
  p1.1 <- mcmc.object$chain1[,"p[1]"] %>% 
    as_tibble() %>% 
    mutate(chain="1",
           descriptor="lgr_p",
           iteration=row_number())
  
  p1.2 <- mcmc.object$chain2[,"p[1]"] %>% 
    as_tibble() %>% 
    mutate(chain="2",
           descriptor="lgr_p",
           iteration=row_number())
  
  p2.1 <- mcmc.object$chain1[,"p[2]"] %>% 
    as_tibble() %>% 
    mutate(chain="1",
           descriptor="lgs_p",
           iteration=row_number())
  
  p2.2 <- mcmc.object$chain2[,"p[2]"] %>% 
    as_tibble() %>% 
    mutate(chain="2",
           descriptor="lgs_p",
           iteration=row_number())
  
  p3.1 <- mcmc.object$chain1[,"p[3]"] %>% 
    as_tibble() %>% 
    mutate(chain="1",
           descriptor="lomo_p",
           iteration=row_number())
  
  p3.2 <- mcmc.object$chain2[,"p[3]"] %>% 
    as_tibble() %>% 
    mutate(chain="2",
           descriptor="lomo_p",
           iteration=row_number())
  
  p4.1 <- mcmc.object$chain1[,"p[4]"] %>% 
    as_tibble() %>% 
    mutate(chain="1",
           descriptor="ich_p",
           iteration=row_number())
  
  p4.2 <- mcmc.object$chain2[,"p[4]"] %>% 
    as_tibble() %>% 
    mutate(chain="2",
           descriptor="ich_p",
           iteration=row_number())
  
  p5.1 <- mcmc.object$chain1[,"p[5]"] %>% 
    as_tibble() %>% 
    mutate(chain="1",
           descriptor="mcn_p",
           iteration=row_number())
  
  p5.2 <- mcmc.object$chain2[,"p[5]"] %>% 
    as_tibble() %>% 
    mutate(chain="2",
           descriptor="mcn_p",
           iteration=row_number())
  
  p6.1 <- mcmc.object$chain1[,"p[6]"] %>% 
    as_tibble() %>% 
    mutate(chain="1",
           descriptor="jda_p",
           iteration=row_number())
  
  p6.2 <- mcmc.object$chain2[,"p[6]"] %>% 
    as_tibble() %>% 
    mutate(chain="2",
           descriptor="jda_p",
           iteration=row_number())
  
  p7.1 <- mcmc.object$chain1[,"p[7]"] %>% 
    as_tibble() %>% 
    mutate(chain="1",
           descriptor="bon_p",
           iteration=row_number())
  
  p7.2 <- mcmc.object$chain2[,"p[7]"] %>% 
    as_tibble() %>% 
    mutate(chain="2",
           descriptor="bon_p",
           iteration=row_number())
  
  p8.1 <- mcmc.object$chain1[,"p[8]"] %>% 
    as_tibble() %>% 
    mutate(chain="1",
           descriptor="twx_p",
           iteration=row_number())
  
  p8.2 <- mcmc.object$chain2[,"p[8]"] %>% 
    as_tibble() %>% 
    mutate(chain="2",
           descriptor="twx_p",
           iteration=row_number())
  
  bind_rows(phi1.1,phi1.2,phi2.1,phi2.2,
            phi3.1,phi3.2,phi4.1,phi4.2,
            phi5.1,phi5.2,phi6.1,phi6.2,
            phi7.1,phi7.2,phi8.1,phi8.2,
            p1.1,p1.2,
            p2.1,p2.2,p3.1,p3.2,p4.1,p4.2,
            p5.1,p5.2,p6.1,p6.2,p7.1,p7.2,
            p8.1,p8.2)
  
  
}




# Precompile NIMBLE model 

compiled_model <- nimbleModel(code=hmm.phitpt,
                              constant=constants.example,
                              data=example.list,
                              inits=initial.values())

compiled_model_c <- compileNimble(compiled_model)

compiled_mcmc <- buildMCMC(compiled_model)

compiled_mcmc_c <- compileNimble(compiled_mcmc,project=compiled_model_c)

# get all distinct groups to generate estimates for

group_grid <- ch.dat25 |> 
  distinct(hatchery,species,release_group,release_year)

# wrapper to run one model

run_nimble_group <- function(hatchery.name,spp,release.group,release.year,
                             input.dat){
  
  dat <- nimble_dat.f(hatchery.name,release.year,release.group,spp,
                      input.dat=ch.dat25)
  
  first <- apply(dat,1,function(x) min(which(x != 0)))
  
  constants <- list(N=nrow(dat),T=ncol(dat),first=first)
  
  data_list <- list(y=dat+1)
  
 zinits <- dat +1
 
 zinits[zinits==2] <- 1
 
 inits <- list(phi=runif(constants$T-1),
               p=runif(constants$T-1),
               z=zinits)
 
 c_model$setConstants(constants)
 c_model$setData(data_list)
 c_model$setInits(inits)
 
 samples <- runMCMC(c_mcmc,
                    niter=n.iter,
                    nburnin=n.burnin,
                    nchains=n.chains,
                    samplesAsCodaMCMC=TRUE,
                    summary=FALSE)
 
 return(samples)
  
  
}



# run the nimble code (takes a while)

tic()
example.phitpt <- nimbleMCMC(code=hmm.phitpt,
                             constants=constants.example,
                             data=example.list,
                             inits=initial.values,
                             monitors=parameters.to.save,
                             niter=n.iter,
                             nburnin=n.burnin,
                             nchains=n.chains)
toc()



# look through outputs/diagniostics

MCMCsummary(example.phitpt,params=c("phi","p"),round=2)

MCMCplot(object=example.phitpt,params="phi")


MCMCtrace(object=example.phitpt,
          pdf=F,ind=TRUE,
          params="phi")

MCMCtrace(object=example.phitpt,
          pdf=F,ind=TRUE,
          params="p")



example.chains <- chains_extract.f(mcmc.object = example.phitpt)
