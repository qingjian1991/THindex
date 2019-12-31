# functions for TH indices estimating.

diversity_partition = function( data  , method = "diversity", region_num = 10, minVaf = 0, maxVaf = 1){
  # devide the VAF into {region_num} regions and calculate the diversity

  if( method !="diversity" &  method != "taxonomic" ){
    stop("The available indicies include diveristy or taxonomic Indicies")
  }

  Summ = data %>%
    mutate(Group =  cut( t_vaf, breaks = seq( minVaf, maxVaf, length.out = region_num+1), labels = seq(1, region_num ) ) ) %>%
    group_by(Group) %>%
    summarise(Num = n() ) %>%
    mutate(prop = Num / sum(Num),
           Group = as.numeric(Group)
    )

  a = Summ$prop

  if(method == "diversity"){
    index1 = -sum(a*log(a)) #Shannon Index.
    index2 = 1/sum(a^2)    # reverse-simpson Index.

    index = c(index1, index2); names(index) = c("shannon","simpson")

  }else if(method == "taxonomic"){

    SummD = Summ %>% dplyr::select(Group, Num) %>%
      spread(key = Group, value = Num )
    x = 1:region_num; names(x) = x; dist = dist(x) #distance.

    taxondis = taxondive(SummD, x )
    index = c(taxondis$D, taxondis$Dstar ); names(index) = c("Delt","Dstar")
  }
  list(index = index, Summary = Summ, method = method)
}
