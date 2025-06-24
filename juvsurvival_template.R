library(tidyverse)
library(RMark)

# the inputs here are capture histories which are compiled using
# code that's in my DNFH_Dashboard repository; basically for each 
# fish we look in PTAGIS for detections at each of the mainstem
# dams and the towed estuary array; also, release groups are
# distinguished (hatchery, release site, release date) so
# they can be treated separately for survival estimates

# using 2025 data as an example of how to run these estimates, 
# reading those in here

ch.dat25 <- read_rds("data/travel")

# for now, just demonstrating one group, so filter steelhead
# released from DNFH in the main clearwater

dnfh_sthd25 <- ch.dat25 |> 
  filter(species=="Steelhead",
         hatchery=="DWOR",
         release_sitecode=="DWORMS")

# the capture history is a series of 1's and 0's 
# representing either observed or not at each dam;
# there are different approaches to what to include but
# here I'm mirroring what the typical DART process 
# uses which is 7 of 8 dams (no Dalles bc low detections)
# and the towed array; the columns
# of 1's and 0's for observation sites
# are already present in this data frame
# from the PIT tag processing code so they just need to 
# be pulled together in a single column; note that i'm 
# adding a leading 1 which represents the initial releast
# occasion

dnfh_sthd25 <- dnfh_sthd25 |> 
  mutate(capture_history=str_c("1",LGR_ch,LGS_ch,LOMO_ch,ICH_ch,
                               MCN_ch,JD_ch,BONN_ch,
                               TWX_ch,sep=""))

# CJS models will be run in RMark which uses 
# Maximum Likelihood Estimation to estimate parameters;
# the first step needed is to define models, which 
# RMark handles in a list format. For this analysis I defined
# fully time-dependent survival and detection probabilities; the
# following code makes list objects for those parameters using
# the naming convention in MARK of survival as "phi" and 
# detection probability as "p"

phi.time <- list(formula=~time)
p.time <- list(formula=~time)

# RMark also expects a couple of other objects. First there is the input
# data which will just be a column for the capture histories and another
# for group, although the way we're applying these models won't really use 
# the group

inp.dat <- dnfh_sthd25 |> 
  ungroup() |> 
  select(ch=capture_history,group=release_group)

# The input data also needs to be processed by the RMArk function
# process.data() which essentially translates the input data 
# into a list with formatting and inputs (such as the number of occasions)
# expected by MARK; that piece is also used in a different function
# that feeds into the MARK framework by using the make.design.data
# function

inp.processed <- process.data(inp.dat)

inp.ddl <- make.design.data(inp.processed)

# Now the model can be run using the mark() function. The 
# output of this will be a list which is saved here
# in the mark.output object

mark.output <- mark(inp.processed,inp.ddl,
                    model.parameters = list(Phi=phi.time,
                                            p=p.time),
                    delete = TRUE)

# The output list has lots of information that can be used in things
# like model selection; the main estimates I tend to use come from the
# results$real table which is the parameter estimates and confidence
# intervals for each time step

mark.output$results$real

# I also like to put these outputs into a format that
# identifies what the parameters are. For this model configuration,
# the outputs are going to follow the same format as subsequent models
# area run so here I'll make a key to assign names to the parameters
# for given rows in the results output

row_key <- tibble(rowname=seq(1,16,1),
                  descriptor=c("release_lgr_phi",
                               "lgr_lgs_phi",
                               "lgs_lomo_phi",
                               "lomo_iha_phi",
                               "iha_mcn_phi",
                               "mcn_jda_phi",
                               "jda_bon_phi",
                               "bon_twx_phi",
                               "lgr_p","lgs_p",
                               "lomo_p","iha_p",
                               "mcn_p","jda_p",
                               "bon_p","twx_p"),
                  category=c("phi","phi","phi","phi",
                             "phi","phi","phi","phi",
                             "p","p","p","p","p","p","p","p")) %>% 
  mutate(rowname=as.character(rowname))

# now the model results can be joined to that key so names 
# are associated with parameters

named.results <- mark.output$results$real |> 
  remove_rownames() |> 
  rownames_to_column() |> 
  left_join(row_key,by="rowname") |> 
  select(descriptor,category,estimate,se,lcl,ucl)

