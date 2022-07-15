library(tidyverse)
library(here)

## This input file is a tabulized version of the blast output. It has all the fields from 
## Blast -outfmt 6 plus the following taxonomy fields: phylum, class, order, family, genus, species
## Each row is a match


cleaned.input <- read_csv(here("test_files", "cleaned_input.test.csv"))



### Then for each query we will try to detect outliers in their matches
### We turn the taxonomy into numbers -> so crustacea:hexanauplia:Metrinidae is 1:2:2 and then 111 and 
### Crustacea:decapoda:cancridae 1:1:1 and then 111
### It is so blunt I am ashamed of sharing it - the idea is that for a given query, matches that have been
### missasigned will be 2:3:3, and then 233 - and in a distribution they will look like outliers


cleaned.input %>%
  group_by(qseqid) %>% 
  nest() %>% 
  mutate(data = map(data, function(.x){
                  .x %>% 
                      mutate(across(phylum:genus,.fns = ~as.numeric(as.factor(.x)))) %>% 
                      unite(phylum:genus, col = "num", sep = "") %>% 
                      mutate(num = as.numeric(num)) -> temp 
                    
   if(n_distinct(temp$num) == 1 ){ # If there is only one match, keep it
     .x %>% 
       mutate(probs = 1,
              keepers = "Keep")
   }else{  # If there are more than one, then fit the values to a normal distribution and label outliers
    
    
    pnorm(q = temp$num,
          mean = mean(temp$num, na.rm = T), 
          sd = 10) -> probas 

      .x %>%
        mutate(probs  =  probas,
               keepers = case_when(probs < 0.95 ~ "Keep",
                                   TRUE         ~ "No Keep"))
     } 
    })) -> output_with_decision

### From here I had a df with the discarded and had a look at it.
### And another one with all the ones I should keep

output_with_decision %>% 
  unnest(data) %>% 
  filter(keepers == "Keep") -> good.matches

output_with_decision %>% 
  unnest(data) %>% 
  filter(keepers != "Keep") -> bad.matches



## Let's see which ones are we keeping
output_with_decision %>% unnest() %>% filter (keepers == "No Keep") %>% distinct(qseqid) -> interesting


output_with_decision %>% 
  ungroup() %>% 
  semi_join(interesting) %>% 
  mutate(overlook = map(data, ~.x %>% group_by(keepers) %>% nest())) %>%
  select(-data) %>% 
  unnest(overlook) -> check 
  
check %>% 
  mutate(data = case_when(keepers == "Keep" ~ map(data, ~.x %>%
                                                    select(phylum:genus) %>%
                                                    taxonomizr::condenseTaxa()%>% 
                                                    set_names(nm = c("phylum", "class", "order", "family", "genus")) %>%
                                                    bind_cols()),
                          TRUE              ~ map(data, ~.x %>% 
                                                    select(phylum:genus) %>%
                                                    taxonomizr::condenseTaxa()%>% 
                                                    set_names(nm = c("phylum", "class", "order", "family", "genus")) %>%
                                                    bind_cols()))) %>% 
  unnest(data) -> to.view

## are the consensus different

to.view %>% 
  unite(phylum:genus, col = "taxa", sep = "|") %>% 
  pivot_wider(names_from = keepers, values_from = taxa) %>% 
  group_by(Keep == `No Keep`) %>% 
  tally
 
## It is not really working

