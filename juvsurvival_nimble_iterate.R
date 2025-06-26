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
  distinct(hatchery.name=hatchery,
           spp=species,
           release.group=release_group,
           release.year=release_year) 


plan(multisession,workers=15)


tic()
test_run <- future_pmap(test_group,
                        run_nimble_group.f,
                        input.dat=ch.dat25,
                        .progress = TRUE)
toc()


MCMCsummary(test_run[[2]],round=2)

# get the posterior summaries from that run

test_group_ref <- test_group |> 
  mutate(group_index=row_number())

test_summary.tbl <- map2_dfr(
  .x=test_run,
  1:length(test_run),
  ~ MCMCsummary(.x, params=c("phi","p"),round=2) |> 
    as_tibble(rownames="paramaeter") |> 
    mutate(group_index=.y)
)#|> 
  left_join(test_group_ref,by="group_index")




chains_extract.f <- function(mcmc_object, group_index) {
  mcmc_object %>%
    map_dfr(as.data.frame, .id = "chain") %>%
    pivot_longer(-chain, names_to = "param", values_to = "value") %>%
    group_by(chain, param) %>%
    mutate(iteration = row_number()) %>%
    ungroup() %>%
    mutate(group_index = group_index)

}


# Join metadata
all_chains_tbl <-map2_dfr(test_run,seq_along(test_run), chains_extract.f) |> 
  left_join(test_group_ref,by="group_index")

saveRDS(test_summary.tbl,
        "outputs/parameter_estimates_2025")

saveRDS(all_chains_tbl,
        "outputs/chains_2025")

# plot chain diagnostics, the main ones 
# to look at are at LGR; chinook only here,
# but steelhead are in the data set



lgr_survival.chains <- all_chains_tbl |> 
  filter(param=="phi[1]",
         spp=="Chinook") |> 
  ggplot(aes(x=iteration,y=value,color=chain))+
  geom_line(alpha=0.8)+
  facet_grid(hatchery.name~release.group)+
  theme_bw()
lgr_survival.chains


lgr_detection.chains <- all_chains_tbl |> 
  filter(param=="p[1]",
         spp=="Chinook") |> 
  ggplot(aes(x=iteration,y=value,color=chain))+
  geom_line(alpha=0.8)+
  facet_grid(hatchery.name~release.group)+
  theme_bw()
lgr_detection.chains



