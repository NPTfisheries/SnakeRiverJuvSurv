
#### Note that Program MARK needs to be 
#### installed for this code to run as
#### it will call program MARK from 
#### within R


library(tidyverse)
library(RMark)

capture_history <- read_rds("juv_hydro_25") %>% 
  mutate(capture_history=str_c("1",LGR_ch,LGS_ch,LOMO_ch,ICH_ch,
                               MCN_ch,JD_ch,BONN_ch,
                               TWX_ch,sep=""))

# compute summaries of detections by dam within 
# release groups...this isn't necessary to run
# the survival estimates


detection_summary <- capture_history %>% 
  group_by(hatchery,species,release_group,release_year) %>% 
  summarize(total_marked=n(),
            LGR_count=sum(LGR_ch=="1"),
            LGS_count=sum(LGS_ch=="1"),
            LOMO_count=sum(LOMO_ch=="1"),
            ICH_count=sum(ICH_ch=="1"),
            MCN_count=sum(MCN_ch=="1"),
            JD_count=sum(JD_ch=="1"),
            TDA_count=sum(TDA_ch=="1"),
            BON_count=sum(BONN_ch=="1"),
            TWX_count=sum(TWX_ch=="1")) %>% 
  mutate(across(
    -total_marked,
    ~ .x/total_marked*100,
    .names="{.col}_percent"
  ))

# define CJS models options to be run in RMark;
# the first step needed is to define models, which 
# RMark handles in a list format. For this analysis I defined
# fully time-dependent survival and detection probabilities; the
# following code makes list objects for those parameters using
# the naming convention in MARK of survival as "phi" and 
# detection probability as "p"

phi.time <- list(formula=~time)
p.time <- list(formula=~time)

# I also like to put these outputs into a format that
# identifies what the parameters are. For this model configuration,
# the outputs are going to follow the same format as subsequent models
# are run so here I'll make a key to assign names to the parameters
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

# with multiple groups and years it's worth defining
# a function to iterate over; these are the steps
# to actually run the fully time-dependent CJS models
# in MARK for each release group/year; the function will
# also match those outputs to my key so the resulting 
# paramater estimates come through in a format that's
# easy to deal with for plotting, etc. It also calculates
# total phi as the product of all the interval survival
# estimates excluding the final one (bon_txw) because
# the final location can't separately estimate 
# phi and p

cjs.f <- function(hatchery.name,
                  releaseyear,spp,
                  releasegroup,input.data){
  
  inp1 <-input.data %>% 
    ungroup() %>% 
    filter(hatchery==hatchery.name,
           release_year==releaseyear&
             release_group==releasegroup&
             species==spp) %>% 
    select(ch=capture_history,group=release_group)
  
  inp.processed <- process.data(inp1, model="CJS")
  
  inp.ddl <- make.design.data(inp.processed)
  
  time.cjs <- mark(inp.processed,
                   inp.ddl,
                   model.parameters = list(Phi=phi.time,
                                           p=p.time),
                   delete=TRUE)
  
  results1 <- time.cjs$results$real %>% 
    remove_rownames() %>% 
    rownames_to_column() %>% 
    left_join(row_key,by="rowname") %>% 
    select(descriptor,category,estimate,se,lcl,ucl)
  
  totalphi.df <- results1 %>% 
    filter(category=="phi"&
             !descriptor=="bon_twx_phi") %>% 
    summarize(estimate=prod(estimate)) %>% 
    mutate(descriptor="overall_phi",
           category="phi",
           se=as.numeric(NA),
           lcl=as.numeric(NA),
           ucl=as.numeric(NA)) %>% 
    select(descriptor,category,estimate,se,lcl,ucl)
  
  finaljoin.df <- results1 %>% 
    filter(descriptor %in% c("bon_twx_phi","twx_p")) %>% 
    summarize(estimate=prod(estimate)) %>% 
    mutate(descriptor="final_joint",
           category="final",
           se=as.numeric(NA),
           lcl=as.numeric(NA),
           ucl=as.numeric(NA)) %>% 
    select(descriptor,category,estimate,se,lcl,ucl)
  
  results2 <- results1 %>% 
    filter(!descriptor %in% c("bon_twx_phi","twx_p")) %>% 
    bind_rows(totalphi.df,finaljoin.df) %>% 
    mutate(hatchery=hatchery.name,species=spp,
           release_year=releaseyear,release_group=releasegroup,
           estimate_source="RMark") %>% 
    select(hatchery,species,release_year,release_group,
           descriptor,category,estimate,se,lcl,ucl,
           estimate_source)
  
  
}

# to iterate over all inputs set up a grid of the possible
# combinations

iteration.vars <- capture_history |> 
  group_by(hatchery,release_year,
           species,release_group) |> 
  tally() |> 
  ungroup() |> 
  select(hatchery.name=hatchery,
         releaseyear=release_year,
         spp=species,
         releasegroup=release_group)

# now run over all of those groups

parameter_estimates <- pmap(iteration.vars,
                             cjs.f,
                             input.data=capture_history) %>% 
  bind_rows()
