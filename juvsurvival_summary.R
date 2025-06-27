

library(tidyverse)
library(arrow)

# read in nimble parameter estimates

nimble25.parms <- read_rds("outputs/parameter_estimates_2025")

nimble25.chains <- read_rds("outputs/chains_2025")

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

# need to get release group, hatchery, species, as well as
# parm desriptors into nimble parameter estimates

nimble25.id <- nimble25.chains |> 
  distinct(hatchery.name,
           spp,
           release.group,
           release.year,group_index) 


nimble25.parms1 <- nimble25.parms |> 
  group_by(group_index) |> 
  mutate(rowname=as.character(row_number())) |> 
  left_join(row_key,by="rowname") |> 
  left_join(nimble25.id,by="group_index") |> 
  ungroup() |> 
  select(hatchery=hatchery.name,species=spp,
         release_group=release.group,
         release_year=release.year,
         descriptor,category,estimate=`50%`,
         lcl=`2.5%`,ucl=`97.5%`) |> 
  mutate(estimate_source="NIMBLE") 


nimble_totalphi.df <-nimble25.parms1 %>% 
  filter(category=="phi"&
           !descriptor=="bon_twx_phi") %>% 
  group_by(hatchery,species,release_group,
           release_year,estimate_source) |> 
  summarize(estimate=prod(estimate)) %>% 
  mutate(descriptor="overall_phi",
         category="phi",
         lcl=as.numeric(NA),
         ucl=as.numeric(NA)) %>% 
  select(hatchery,species,release_group,release_year,
         descriptor,category,estimate,lcl,ucl,estimate_source)



nimble_finaljoint.df <- nimble25.parms1 %>% 
  filter(descriptor %in% c("bon_twx_phi","twx_p")) %>% 
  group_by(hatchery,species,release_group,
           release_year,estimate_source) |> 
  summarize(estimate=prod(estimate)) %>% 
  mutate(descriptor="final_joint",
         category="final",
         lcl=as.numeric(NA),
         ucl=as.numeric(NA)) %>% 
  select(hatchery,species,release_group,release_year,
         descriptor,category,estimate,lcl,ucl,estimate_source)


nimble25.parms2 <-nimble25.parms1 %>% 
  filter(!descriptor %in% c("bon_twx_phi","twx_p")) %>% 
  bind_rows(nimble_totalphi.df,nimble_finaljoint.df) %>% 
  select(hatchery,species,release_group,release_year,
         descriptor,category,estimate,lcl,ucl,estimate_source)

# read in rmark parameter estimates

rmark25.parms <- read_rds("outputs/rmark_parameter_estimates_2025") |> 
  select(-se)


# bind together estimate sources for comparisons

source_comp25 <- rmark25.parms |> 
  bind_rows(nimble25.parms2) |> 
  mutate(comp_group=str_c(hatchery,species,release_group,sep="_"))

# plot lgr survival and detection between estimate sources

lgr_phi.plot <- source_comp25 |> 
  filter(descriptor=="release_lgr_phi") |> 
  ggplot(aes(x=comp_group,y=estimate,
             color=estimate_source))+
  geom_point(position=position_dodge(width=0.4),
             size=3)+
  geom_errorbar(aes(ymin=lcl,ymax=ucl),
                position=position_dodge(width=0.4))+
  theme_bw()+
  scale_color_manual(values=c("red","blue"))+
  labs(x="Release Group",
       y="Survival Estimate",
       color="Method")+
  theme(axis.text.x = element_text(angle=45,hjust=1))
lgr_phi.plot
