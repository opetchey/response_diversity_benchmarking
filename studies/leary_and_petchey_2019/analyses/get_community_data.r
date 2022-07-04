rm(list=ls())

library(tidyverse)
library(here)

pop.data <- read.csv(here("data/raw/raw_community_data.csv"))
biomass <- read.csv("data/raw/cell_volumes.csv")
tt.rr <- read.csv(here("data/raw/expt micro temp ranges.csv"))

##pop.data$count = as.numeric(pop.data$count)

pop.data = transform(pop.data, per.ml=count/weight.1*weight.2/weight.3*weight.4/weight.5)

summary(pop.data)

####first edit of pop.data (jars that I forgot to count)
m.42 = subset(pop.data, microcosm==42)####paramecium not counted on day 27
m.24 = subset(pop.data, microcosm==24)####no count for tetrahymena day 43
m.36 = subset(pop.data, microcosm==36)
unique(m.42$species)
unique(m.24$species)
unique(m.36$species)
edit.1 = 	pop.data$microcosm==24 & pop.data$day==43 |
		pop.data$microcosm==42 & pop.data$day==27
pop.data = pop.data[!edit.1,]

########including biomass data
mean.biomass <- tapply(biomass$volume.mg.,
			list(biomass$species),
			mean)
mean.biomass <- data.frame(mean.biomass)
species <- cbind(unique(biomass$species))
species <- data.frame(species)
biomass <- data.frame(mean.biomass, species)
pop.data <- data.frame(pop.data, biomass=biomass
		[match(pop.data[,3], biomass[,2]),1])

####calculating biomass of each pop count
pop.data <- transform(pop.data, pop.biomass=biomass*per.ml)
pop.data <- transform(pop.data, log.biomass=log10(pop.biomass))

####
pop.data.edit.1 = 	pop.data$microcosm==36 & pop.data$day==51 |
			pop.data$microcosm==24 & pop.data$day==43 |
			pop.data$microcosm==42 & pop.data$day==27
pop.data = pop.data[!pop.data.edit.1,]
pop.data.edit.2 = 	pop.data$microcosm==2 |
			pop.data$microcosm==8 |
			pop.data$microcosm==14 |
			pop.data$microcosm==20 |
			pop.data$microcosm==26 |
			pop.data$microcosm==32 |
			pop.data$microcosm==38 |
			pop.data$microcosm==44 |
			pop.data$microcosm==50 
pop.data = pop.data[!pop.data.edit.2,]


#### First match up the microcosm number to correlation value
spp <- unique(pop.data$species)
species.codes <- data.frame(spp=spp,
                            code=c("C", "L", "T", "P"))

all.jars <- unique(pop.data$microcosm)

counter = 1
for(i in 1:length(all.jars)) {
    spps <- unique(pop.data$species[pop.data$microcosm==all.jars[i]])
    if(counter==1)
        jar.spp <- data.frame(microcosm=all.jars[i],
                              spp1=spps[1],
                              spp2=spps[2])
    if(counter>1)
        jar.spp <- rbind(jar.spp, data.frame(microcosm=all.jars[i],
                                                  spp1=spps[1],
                                                  spp2=spps[2]))
    counter = counter + 1
}

## stick the temp ranges on
temp.ranges <- tt.rr[match(jar.spp[,1], tt.rr[,1]),2]
jar.spp <- cbind(jar.spp, temp.ranges)

## and the compositions
m.code <- unique(paste(species.codes[match(jar.spp[,2], species.codes[,1]),2],
                       species.codes[match(jar.spp[,3], species.codes[,1]),2], sep=""))
jar.spp <- cbind(jar.spp, m.code)

## code.match <- data.frame(c1=m.code,
##                          c2=c("LC", "TC", "LP", "TL", "TP"))




## now get the cov of total biomasses
all.jars
covs <- data.frame(microcosm=numeric(length=length(all.jars)),
                   cov=numeric(length=length(all.jars)))
vars <- data.frame(microcosm=numeric(length=length(all.jars)),
                   cov=numeric(length=length(all.jars)))
cors <- data.frame(microcosm=numeric(length=length(all.jars)),
                   cor=numeric(length=length(all.jars)))


##par(ask=T)
for(i in 1:length(all.jars)) {
    sp.in.jar <- as.character(unique(pop.data$species[pop.data$microcosm==all.jars[i]]))
    these.x <- pop.data$species==sp.in.jar[1] & pop.data$microcosm==all.jars[i]
    these.y <- pop.data$species==sp.in.jar[2] & pop.data$microcosm==all.jars[i]
    ##plot(pop.data$pop.biomass[these.x], pop.data$pop.biomass[these.y]) 
    covs[i,1] <- all.jars[i]
    covs[i,2] <- cov(pop.data$pop.biomass[these.x], pop.data$pop.biomass[these.y]) 
    vars[i,1] <- all.jars[i]
    vars[i,2] <- var(pop.data$pop.biomass[these.x]) + var(pop.data$pop.biomass[these.y])
    cors[i,1] <- all.jars[i]
    cors[i,2] <- cor(pop.data$pop.biomass[these.x], pop.data$pop.biomass[these.y])

}
    
total.biomass = aggregate(pop.data$pop.biomass,
    list(day=pop.data$day,
         microcosm=pop.data$microcosm),
    sum)

cv.biomass <- aggregate(total.biomass$x,
                        list(microcosm=total.biomass$microcosm),
                        function(x) sd(x)/mean(x))

total.biomass = aggregate(pop.data$pop.biomass,
    list(microcosm=pop.data$microcosm),
    mean)


## also check out evenness
evenness = aggregate(pop.data$pop.biomass,
    list(day=pop.data$day,
         microcosm=pop.data$microcosm),
    function(x) 1/sum((x/sum(x))^2)/2 )

mean.eve <- aggregate(evenness$x,
                      list(microcosm=evenness$microcosm),
                      mean)



## and add these to the jar.spp data

jar.spp <- cbind(jar.spp, cv=cv.biomass[match(jar.spp$microcosm, cv.biomass$microcosm),2])

jar.spp <- cbind(jar.spp, cov=covs[match(jar.spp$microcosm, covs$microcosm),2])

jar.spp <- cbind(jar.spp, cor=cors[match(jar.spp$microcosm, cors$microcosm),2])

jar.spp <- cbind(jar.spp, eve=mean.eve[match(jar.spp$microcosm, mean.eve$microcosm),2])

jar.spp <- cbind(jar.spp, var=vars[match(jar.spp$microcosm, vars$microcosm),2])

jar.spp <- cbind(jar.spp, tot=total.biomass[match(jar.spp$microcosm, total.biomass$microcosm),2])

str(jar.spp)

community_data <- jar.spp %>%
  mutate(study_name = "Leary and Petchey 2009",
         community_name = paste0("community-", microcosm)) %>%
  rename(total_biomass = tot) %>%
  pivot_longer(names_to = "community_metric_name",
               values_to = "community_metric_value",
               cols = 6:11) %>%
  select(study_name, community_name, community_metric_name, community_metric_value)
write_csv(community_data, here("data/derived1/community_data.csv"))
         
community_composition_data <- jar.spp %>%
  mutate(study_name = "Leary and Petchey 2009",
         community_name = paste0("community-", microcosm)) %>%
  select(study_name, community_name, spp1, spp2) %>%
  pivot_longer(names_to = "junk", values_to = "species_name", cols = 3:4) %>%
  select(-junk)
write_csv(community_composition_data, here("data/derived1/community_composition_data.csv"))

environment_data <- jar.spp %>%
  mutate(study_name = "Leary and Petchey 2009",
         community_name = paste0("community-", microcosm)) %>%
  select(study_name, community_name, temp.ranges) %>%
  rowwise() %>%
  mutate(temp = case_when(temp.ranges == 1 ~ "18-20-22",
                          temp.ranges == 2 ~ "20-22-24",
                          temp.ranges == 3 ~ "22-24-26")) %>%
  separate(temp, into = c("t1", "t2", "t3"), sep = "-") %>%
  pivot_longer(names_to = "temp", values_to = "environment_value", cols = 4:6) %>%
  mutate(environment_time = parse_number(temp),
         environment_value = parse_number(environment_value),
         environment_name = "temperature") %>%
  select(-temp, -temp.ranges)
write_csv(environment_data, here("data/derived1/environment_data.csv"))




