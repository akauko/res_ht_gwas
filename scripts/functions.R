

create_n_table <- function(my.df, var, covs_v, names_ep_l){
  
  my.df <- my.df %>% 
    filter(!is.na(get(var))) %>%
    group_by_at(var)
  
  ns.n.tmp <- my.df %>% 
    summarise(n=n()) 
  
  ns.all.tmp <- my.df %>%
    summarise_at(c(covs_v, unlist(names_ep_l)), 
                 ~str_glue("{sum(.,na.rm=T)} ({round(mean(., na.rm=T)*100,1)}%)")) %>%
    select(-all_of(var))
  
  ns.p.tmp <- lapply(as.list(c(covs_v, unlist(names_ep_l))), function(x){
    chisq.test(table(select(my.df, all_of(c(var,x)) )))$p.val
  }) %>% 
    unlist %>%
    signif(3)
  
  bind_cols(ns.n.tmp, ns.all.tmp) %>% 
    t() %>%
    as.data.frame() %>%
    setNames(ns.n.tmp[[var]]) %>%
    rownames_to_column(var=var) %>% 
    slice(-1) %>%
    mutate(pval=c(NA_real_, ns.p.tmp))  %>% 
    mutate(pval=if_else(pval < 2e-16, "<2e-16", as.character(pval))) %>%
    rename(`P-value`=pval)
}


my.kable <- function(my.table, font=16) {
  my.table %>%
    knitr::kable() %>%
    kable_classic() %>%
    row_spec(0,bold=TRUE) #%>%
    #kable_styling(font_size = font)
  
}


#Draws histograms of endpoints by age. Cases and controls are coloured separately-

my_endpoint_histograms <-  function(endpoints_all, endpoint_ages_all, df=df, ncol=4, xlab="age", ylim=NA){

  df_long_ep <-   df %>%
    select("FINNGENID", all_of(endpoints_all)) %>%
    pivot_longer(cols=all_of(endpoints_all), names_to = "endpoint",values_to = "status")
  
  df_long <-   df %>%
    select("FINNGENID", all_of(endpoint_ages_all)) %>%
    pivot_longer(cols=all_of(endpoint_ages_all), names_to = "endpoint",values_to = "age") %>%
    mutate(endpoint = str_replace(endpoint,"(\\w+)_.+", "\\1")) %>%
    left_join(df_long_ep, by=c("FINNGENID", "endpoint"))
  
  ggplot(df_long, aes(x=age, color=as.factor(status), fill=as.factor(status)), xlab = xlab) +
    geom_histogram(alpha=0.25, position="identity") +
    theme(legend.title = element_blank()) +
    facet_wrap(~endpoint, ncol=ncol) +
    xlab(xlab) + 
    coord_cartesian(ylim=c(0,ylim))
    #ylim(0,ylim)
  
}



#Function to extract tables
extr_table <- function(fit, var_select, outcome="Variable"){ 
  cbind(summary(fit)$conf.int, pval = summary(fit)$coefficients[,5]) %>%
    as.data.frame %>%
    rownames_to_column("Variable") %>%
    rename(est ="exp(coef)", low = "lower .95", high="upper .95") %>%
    filter(str_detect(Variable, var_select)) %>%
    mutate_at(c("est","low","high"), round,2) %>%
    mutate_at("pval", signif,2) %>%
    mutate(HR = str_glue("{est} ({low}-{high})")) %>% 
    select(Variable, HR, pval) %>%
    rename(!!outcome := Variable)
}  



my.coxsnell.plot <- function(fit, var_event, df, title=""){
  
  #Creates coxsnell plot
  #Based on code at this tutorial: https://rpubs.com/kaz_yos/resid_coxot 
  
  df.tmp <- df %>% filter(!is.na(get(var_event)))
  ##Martingale residuals
  df.tmp$resid_mart <- residuals(fit, type = "martingale")
  
  ## Cox-Snell residuals
  df.tmp$resid_coxsnell <- -(df.tmp$resid_mart - df.tmp[[var_event]])
  
  ## Fit model on Cox-Snell residuals (Approximately Expo(1) distributed under correct model)
  fit_coxsnell <- coxph(as.formula(str_glue("Surv(resid_coxsnell, {var_event}) ~ 1")),
                        data    = df.tmp)
  
  ## Nelson-Aalen estimator for baseline hazard (all covariates zero)
  df_base_haz <- basehaz(fit_coxsnell, centered = FALSE)
  
  ## Plot
  ggplot(data = df_base_haz, mapping = aes(x = time, y = hazard)) +
    geom_point() +
    #labs(x = "Cox-Snell residuals as pseudo observed times",
         #y = "Estimated cumulative hazard at pseudo observed times") +
    labs(x = "Cox-Snell residuals", y = "Estimated cumulative hazard") +
    theme_bw() + theme(legend.key = element_blank()) + 
    geom_abline(intercept = 0, slope=1, size=0.5, linetype="dashed") + 
    ggtitle(title)
  
}

load_fg_variables <- function(variables, table_ep, table_min, projectid ){
  #Reads in FinnGen variables from bigrquery SQL-tables
  #Input: variable names as vector. 
  #Creates following columns:
  #FINNGENID, BL_AGE, BL_YEAR, AGE_AT_DEATH_OR_END_OF_FOLLOWUP, SEX, <variable1>, <variable1>_AGE, <variable1>_YEAR,  <variable2>, ....
  #Finngen control exclusions are set as NA as in official flat files.
  #Finngen minimum
  sql<- str_glue("
    SELECT FINNGENID, BL_AGE, BL_YEAR, AGE_AT_DEATH_OR_END_OF_FOLLOWUP, SEX
    FROM `{table_min}`
  ")
  
  tb <- bq_project_query(projectid, sql)
  df_min <- bq_table_download(tb)  %>%
    mutate_at("SEX", as.factor)
  #df_min %>% summary()
  
  #Phenotype data
  df0 <- lapply(variables, function(my.var){
    
    #The file is in long format 
    sql <- str_glue("
      SELECT FINNGENID,  ENDPOINT, AGE,  APPROX_EVENT_DAY, CONTROL_CASE_EXCL
      FROM `{table_ep}`
      WHERE ENDPOINT = '{my.var}'
    ")
    tb <- bq_project_query(projectid, sql)
    bq_table_download(tb)
  }) %>%
    bind_rows() %>%
    #CONTROL_CASE_EXCL: 1 case, 2: excluded from control.
    #We will later convert 2L as NA and set NA's (values absent from here) as 1L 
    rename(EVENT= CONTROL_CASE_EXCL) %>%    
    mutate(YEAR=lubridate::year(APPROX_EVENT_DAY)) %>%
    select(-APPROX_EVENT_DAY) %>%
    #To wide format
    pivot_wider(names_from=ENDPOINT, values_from = c("EVENT", "AGE", "YEAR"),
                names_glue = "{ENDPOINT}_{.value}")
  
  #Combine and create data in 'endpoint format'
  df.tmp <- df_min %>% 
    left_join(df0, by="FINNGENID") %>%
    mutate_at(vars(ends_with("_EVENT")), ~if_else(is.na(.), 0L, .)) %>%            #These were missing from bigrqurery table, set as 0L 
    mutate_at(vars(ends_with("_EVENT")), ~if_else(.==2, NA, .)) %>%                #These are FG control exclusions, set as NA 
    mutate_at(vars(ends_with("_AGE")), ~if_else(is.na(.), AGE_AT_DEATH_OR_END_OF_FOLLOWUP, .)) %>%      #If age is missing, set as FU_END_AGE
    mutate_at(vars(ends_with("_EVENT")), as.factor) %>% 
    rename_all(~str_remove(.,"_EVENT"))
} 


