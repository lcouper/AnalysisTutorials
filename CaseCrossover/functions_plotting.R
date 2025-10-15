make_cclong <- function(casecontr){
  
  ncol <- ncol(casecontr)
  #colnames(casecontr) <- c(colnames(casecontr)[2:ncol], "X")
  
  ## Now make long, make cases 1 and controls 0 ##
  cclong_di0_1 <- casecontr %>% 
    mutate(id = 1:nrow(casecontr)) %>%
    dplyr::select(c("id", "diagDate", "Age", "Sex", "Race_eth", "Race", "Ethnicity",
                    "case_di0_1", "contr_di0_1")) %>%
    pivot_longer("case_di0_1":"contr_di0_1") %>%
    mutate(case = as.numeric(str_detect(name, "case")),
           year = year(diagDate)) %>% data.frame()
  colnames(cclong_di0_1)[which(colnames(cclong_di0_1) == "value")] <- "exp_di0_1"
  
  cclong_di1_2 <- casecontr %>% 
    dplyr::select(c("case_di1_2", "contr_di1_2")) %>%
    pivot_longer("case_di1_2":"contr_di1_2") %>%
    data.frame()
  colnames(cclong_di1_2)[which(colnames(cclong_di1_2) == "value")] <- "exp_di1_2"
  
  cclong_di2_3 <- casecontr %>% 
    dplyr::select(c("case_di2_3", "contr_di2_3")) %>%
    pivot_longer("case_di2_3":"contr_di2_3") %>%
    data.frame()
  colnames(cclong_di2_3)[which(colnames(cclong_di2_3) == "value")] <- "exp_di2_3"
  
  cclong_di3_4 <- casecontr %>% 
    dplyr::select(c("case_di3_4", "contr_di3_4")) %>%
    pivot_longer("case_di3_4":"contr_di3_4") %>%
    data.frame()
  colnames(cclong_di3_4)[which(colnames(cclong_di3_4) == "value")] <- "exp_di3_4"
  
  cclong_di4_5 <- casecontr %>% 
    dplyr::select(c("case_di4_5", "contr_di4_5")) %>%
    pivot_longer("case_di4_5":"contr_di4_5") %>%
    data.frame()
  colnames(cclong_di4_5)[which(colnames(cclong_di4_5) == "value")] <- "exp_di4_5"
  
  cclong_di0_5 <- casecontr %>% 
    dplyr::select(c("case_di0_5", "contr_di0_5")) %>%
    pivot_longer("case_di0_5":"contr_di0_5") %>%
    data.frame()
  colnames(cclong_di0_5)[which(colnames(cclong_di0_5) == "value")] <- "exp_di0_5"
  
  cclong <- cbind(cclong_di0_1, cclong_di1_2, cclong_di2_3, cclong_di3_4,
                  cclong_di4_5, cclong_di0_5)
  
  # Make binary vars
  cclong$exp_di0_1bin <- as.numeric(cclong$exp_di0_1 > 0)
  cclong$exp_di1_2bin <- as.numeric(cclong$exp_di1_2 > 0)
  cclong$exp_di2_3bin <- as.numeric(cclong$exp_di2_3 > 0)
  cclong$exp_di3_4bin <- as.numeric(cclong$exp_di3_4 > 0)
  cclong$exp_di4_5bin <- as.numeric(cclong$exp_di4_5 > 0)
  cclong$exp_di0_5bin <- as.numeric(cclong$exp_di0_5 > 0)
  
  # Make continuous vars from 0 to 5
  cclong$exp_di0_2 <- cclong$exp_di0_1 + cclong$exp_di1_2
  cclong$exp_di0_3 <- cclong$exp_di0_2 + cclong$exp_di2_3
  cclong$exp_di0_4 <- cclong$exp_di0_3 + cclong$exp_di3_4
  
  cclong$exp_di0_2bin <- as.numeric(cclong$exp_di0_2 > 0)
  cclong$exp_di0_3bin <- as.numeric(cclong$exp_di0_3 > 0)
  cclong$exp_di0_4bin <- as.numeric(cclong$exp_di0_4 > 0)
  
  return(cclong)
  
}

make_cclong2 <- function(casecontr){
  
  ncol <- ncol(casecontr)
  
  ## Now make long, make cases 1 and controls 0 ##
  cclong_di0_1 <- casecontr %>% 
    mutate(id = 1:nrow(casecontr)) %>%
    dplyr::select(c("id", "diagDate", "Age", "Sex", "Race_eth", "Race", "Ethnicity",
                    "case_di0_1", "contr1_di0_1", "contr2_di0_1")) %>%
    pivot_longer("case_di0_1":"contr2_di0_1") %>%
    mutate(case = as.numeric(str_detect(name, "case")),
           year = year(diagDate)) %>% data.frame()
  colnames(cclong_di0_1)[which(colnames(cclong_di0_1) == "value")] <- "exp_di0_1"
  
  cclong_di1_2 <- casecontr %>% 
    dplyr::select(c("case_di1_2", "contr1_di1_2", "contr2_di1_2")) %>%
    pivot_longer("case_di1_2":"contr2_di1_2") %>%
    data.frame()
  colnames(cclong_di1_2)[which(colnames(cclong_di1_2) == "value")] <- "exp_di1_2"
  
  cclong_di2_3 <- casecontr %>% 
    dplyr::select(c("case_di2_3", "contr1_di2_3", "contr2_di2_3")) %>%
    pivot_longer("case_di2_3":"contr2_di2_3") %>%
    data.frame()
  colnames(cclong_di2_3)[which(colnames(cclong_di2_3) == "value")] <- "exp_di2_3"
  
  cclong_di3_4 <- casecontr %>% 
    dplyr::select(c("case_di3_4", "contr1_di3_4", "contr2_di3_4")) %>%
    pivot_longer("case_di3_4":"contr2_di3_4") %>%
    data.frame()
  colnames(cclong_di3_4)[which(colnames(cclong_di3_4) == "value")] <- "exp_di3_4"
  
  cclong_di4_5 <- casecontr %>% 
    dplyr::select(c("case_di4_5", "contr1_di4_5", "contr2_di4_5")) %>%
    pivot_longer("case_di4_5":"contr2_di4_5") %>%
    data.frame()
  colnames(cclong_di4_5)[which(colnames(cclong_di4_5) == "value")] <- "exp_di4_5"
  
  cclong_di0_5 <- casecontr %>% 
    dplyr::select(c("case_di0_5", "contr1_di0_5", "contr2_di0_5")) %>%
    pivot_longer("case_di0_5":"contr2_di0_5") %>%
    data.frame()
  colnames(cclong_di0_5)[which(colnames(cclong_di0_5) == "value")] <- "exp_di0_5"
  
  cclong <- cbind(cclong_di0_1, cclong_di1_2, cclong_di2_3, cclong_di3_4,
                  cclong_di4_5, cclong_di0_5)
  
  # Make binary vars
  cclong$exp_di0_1bin <- as.numeric(cclong$exp_di0_1 > 0)
  cclong$exp_di1_2bin <- as.numeric(cclong$exp_di1_2 > 0)
  cclong$exp_di2_3bin <- as.numeric(cclong$exp_di2_3 > 0)
  cclong$exp_di3_4bin <- as.numeric(cclong$exp_di3_4 > 0)
  cclong$exp_di4_5bin <- as.numeric(cclong$exp_di4_5 > 0)
  cclong$exp_di0_5bin <- as.numeric(cclong$exp_di0_5 > 0)
  
  # Make continuous vars from 0 to 5
  cclong$exp_di0_2 <- cclong$exp_di0_1 + cclong$exp_di1_2
  cclong$exp_di0_3 <- cclong$exp_di0_2 + cclong$exp_di2_3
  cclong$exp_di0_4 <- cclong$exp_di0_3 + cclong$exp_di3_4
  
  cclong$exp_di0_2bin <- as.numeric(cclong$exp_di0_2 > 0)
  cclong$exp_di0_3bin <- as.numeric(cclong$exp_di0_3 > 0)
  cclong$exp_di0_4bin <- as.numeric(cclong$exp_di0_4 > 0)
  
  return(cclong)
  
}

run_models <- function(cclong, subset = F, year_subset = NULL, model5 = F){
  
  if (subset == T){
    cclong <- subset(cclong, year >= year_subset)
  }
  # run the conditional logistic regression 
  mod0_1 <- clogit(case ~ exp_di0_1 + strata(id), cclong)
  mod0_2 <- clogit(case ~ exp_di0_2 + strata(id), cclong)
  mod0_3 <- clogit(case ~ exp_di0_3 + strata(id), cclong)
  mod0_4 <- clogit(case ~ exp_di0_4 + strata(id), cclong)
  mod0_5 <- clogit(case ~ exp_di0_5 + strata(id), cclong)
  
  c1 <- data.frame(rbind(summary(mod0_1)$conf.int,
              summary(mod0_2)$conf.int,
              summary(mod0_3)$conf.int,
              summary(mod0_4)$conf.int,
              summary(mod0_5)$conf.int))
  c1$name <- c("0-1 km", "0-2 km", "0-3 km", "0-4 km", "0-5 km")
  c1$model <- 1

  # continuous models
  mod_cont <- clogit(case ~ exp_di0_1 + exp_di1_2 +
                  exp_di2_3 + exp_di3_4 + exp_di4_5 + strata(id), cclong)
  
  c2 <- data.frame(summary(mod_cont)$conf.int)
  c2$name <- c("0-1 km", "1-2 km", "2-3 km", "3-4 km", "4-5 km")
  c2$model <- 2
  
  # binary models
  # run the conditional logistic regression 
  mod0_1bin <- clogit(case ~ exp_di0_1bin + strata(id), cclong)
  mod0_2bin <- clogit(case ~ exp_di0_2bin + strata(id), cclong)
  mod0_3bin <- clogit(case ~ exp_di0_3bin + strata(id), cclong)
  mod0_4bin <- clogit(case ~ exp_di0_4bin + strata(id), cclong)
  mod0_5bin <- clogit(case ~ exp_di0_5bin + strata(id), cclong)
  
  mod_contbin <- clogit(case ~ exp_di0_1bin + exp_di1_2bin +
                  exp_di2_3bin + exp_di3_4bin + exp_di4_5bin + strata(id), cclong)
  
  c3 <- data.frame(rbind(summary(mod0_1bin)$conf.int,
                         summary(mod0_2bin)$conf.int,
                         summary(mod0_3bin)$conf.int,
                         summary(mod0_4bin)$conf.int,
                         summary(mod0_5bin)$conf.int))
  c3$name <- c("0-1 km", "0-2 km", "0-3 km", "0-4 km", "0-5 km")
  c3$model <- 3
  
  c4 <- data.frame(summary(mod_contbin)$conf.int)
  c4$name <- c("0-1 km", "1-2 km", "2-3 km", "3-4 km", "4-5 km")
  c4$model <- 4
  

  # Quantiles
  if (model5 == T){
    breaks <- quantile(cclong$exp_di0_5[cclong$exp_di0_5 >0])
    breaks <- c(0, breaks)
    
    cclong$diQ <- cut(cclong$exp_di0_5, breaks = breaks, include.lowest = T)
    
    modQ <- clogit(case ~ diQ + strata(id), cclong)
    
    c5 <- data.frame(summary(modQ)$conf.int)
    c5$name <- c("Q1", "Q2", "Q3", "Q4")
    c5$model <- 5
    
    coefs <- rbind(c1,c2,c3,c4,c5)
  } else {
    coefs <- rbind(c1,c2,c3,c4)
  }
  
  return(coefs)

}

run_models_with_neg_control <- function(cclong, subset = F, year_subset = NULL, model5 = T){
  
  if (subset == T){
    cclong <- subset(cclong, year >= year_subset)
  }
  # run the conditional logistic regression 
  mod0_1 <- clogit(case ~ exp_di0_1 + exp_di0_1NC + strata(id), cclong)
  mod0_2 <- clogit(case ~ exp_di0_2 + exp_di0_2NC + strata(id), cclong)
  mod0_3 <- clogit(case ~ exp_di0_3 + exp_di0_3NC + strata(id), cclong)
  mod0_4 <- clogit(case ~ exp_di0_4 + exp_di0_4NC + strata(id), cclong)
  mod0_5 <- clogit(case ~ exp_di0_5 + exp_di0_5NC + strata(id), cclong)
  
  c1 <- data.frame(rbind(summary(mod0_1)$conf.int,
                         summary(mod0_2)$conf.int,
                         summary(mod0_3)$conf.int,
                         summary(mod0_4)$conf.int,
                         summary(mod0_5)$conf.int))
  c1$name <- c("0-1 km", "0-1 km", "0-2 km", "0-2 km", "0-3 km", "0-3 km",
               "0-4 km", "0-4 km", "0-5 km", "0-5 km")
  c1$model <- 1
  c1$Estimate <- rep(c("Observed", "Negative control"), 5)
  
  # continuous models
  mod_cont <- clogit(case ~ exp_di0_1 + exp_di1_2 +
                       exp_di2_3 + exp_di3_4 + exp_di4_5 + 
                       exp_di0_1NC + exp_di1_2NC +
                       exp_di2_3NC + exp_di3_4NC + exp_di4_5NC +
                       strata(id), cclong)

  c2 <- data.frame(summary(mod_cont)$conf.int)
  c2$name <- rep(c("0-1 km", "1-2 km", "2-3 km", "3-4 km", "4-5 km"),2)
  c2$model <- 2
  c2$Estimate <- c(rep("Observed", 5), rep("Negative Control", 5))

  # binary models
  # run the conditional logistic regression
  mod0_1bin <- clogit(case ~ exp_di0_1bin + exp_di0_1binNC + strata(id), cclong)
  mod0_2bin <- clogit(case ~ exp_di0_2bin + exp_di0_2binNC + strata(id), cclong)
  mod0_3bin <- clogit(case ~ exp_di0_3bin + exp_di0_3binNC + strata(id), cclong)
  mod0_4bin <- clogit(case ~ exp_di0_4bin + exp_di0_4binNC + strata(id), cclong)
  mod0_5bin <- clogit(case ~ exp_di0_5bin + exp_di0_5binNC + strata(id), cclong)

  mod_contbin <- clogit(case ~ exp_di0_1bin + exp_di1_2bin + exp_di2_3bin +
                          exp_di3_4bin + exp_di4_5bin +
                          exp_di0_1binNC + exp_di1_2binNC + exp_di2_3binNC +
                          exp_di3_4binNC + exp_di4_5binNC + strata(id), cclong)
  
  c3 <- data.frame(rbind(summary(mod0_1bin)$conf.int,
                         summary(mod0_2bin)$conf.int,
                         summary(mod0_3bin)$conf.int,
                         summary(mod0_4bin)$conf.int,
                         summary(mod0_5bin)$conf.int))
  c3$name <- c("0-1 km", "0-1 km", "0-2 km", "0-2 km", "0-3 km", "0-3 km", 
               "0-4 km", "0-4 km", "0-5 km", "0-5 km")
  c3$model <- 3
  c3$Estimate <- rep(c("Observed", "Negative Control"), 5)

  c4 <- data.frame(summary(mod_contbin)$conf.int)
  c4$name <- rep(c("0-1 km", "1-2 km", "2-3 km", "3-4 km", "4-5 km"),2)
  c4$model <- 4
  c4$Estimate <- c(rep("Observed", 5), rep("Negative Control", 5))

  # Quantiles
  if (model5 == T){
  breaks <- quantile(cclong$exp_di0_5[cclong$exp_di0_5 >0])
  breaks <- c(0, breaks)

  cclong$diQ <- cut(cclong$exp_di0_5, breaks = breaks, include.lowest = T)
  cclong$diQNC<- cut(cclong$exp_di0_5NC, breaks = breaks, include.lowest = T)
  
  modQ <- clogit(case ~ diQ + diQNC + strata(id), cclong)

  c5 <- data.frame(summary(modQ)$conf.int)
  c5$name <- c("Q1", "Q2", "Q3", "Q4", "Q1", "Q2", "Q3", "Q4")
  c5$model <- 5
  c5$Estimate <- c(rep("Observed", 4), rep("Negative Control", 4))

  coefs <- rbind(c1,c2,c3,c4,c5)}
  else {
    coefs <- rbind(c1,c2,c3,c4)
  }
  
  
  return(coefs)
  
}

run_models_with_neg_control_2 <- function(cclong, subset = F, year_subset = NULL, model5 = F){
  # No 0-1km distance
  if (subset == T){
    cclong <- subset(cclong, year >= year_subset)
  }
  # run the conditional logistic regression 
  #mod0_1 <- clogit(case ~ exp_di0_1 + exp_di0_1NC + strata(id), cclong)
  mod0_2 <- clogit(case ~ exp_di0_2 + exp_di0_2NC +strata(id), cclong)
  mod0_3 <- clogit(case ~ exp_di0_3 + exp_di0_3NC + strata(id), cclong)
  mod0_4 <- clogit(case ~ exp_di0_4 + exp_di0_4NC + strata(id), cclong)
  mod0_5 <- clogit(case ~ exp_di0_5 + exp_di0_5NC + strata(id), cclong)
  
  c1 <- data.frame(rbind(#summary(mod0_1)$conf.int,
                         summary(mod0_2)$conf.int,
                         summary(mod0_3)$conf.int,
                         summary(mod0_4)$conf.int,
                         summary(mod0_5)$conf.int))
  c1$name <- c("0-2 km", "0-2 km", "0-3 km", "0-3 km",
               "0-4 km", "0-4 km", "0-5 km", "0-5 km")
  c1$model <- 1
  c1$Estimate <- rep(c("Observed", "Negative control"), 4)
  
  # continuous models
  mod_cont <- clogit(case ~ exp_di1_2 +
                       exp_di2_3 + exp_di3_4 + exp_di4_5 + 
                       exp_di1_2NC +
                       exp_di2_3NC + exp_di3_4NC + exp_di4_5NC +
                       strata(id), cclong)
  
  c2 <- data.frame(summary(mod_cont)$conf.int)
  c2$name <- rep(c("1-2 km", "2-3 km", "3-4 km", "4-5 km"),2)
  c2$model <- 2
  c2$Estimate <- c(rep("Observed", 4), rep("Negative Control", 4))
  
  # binary models
  # run the conditional logistic regression
  #mod0_1bin <- clogit(case ~ exp_di0_1bin + exp_di0_1binNC + strata(id), cclong)
  mod0_2bin <- clogit(case ~ exp_di0_2bin + exp_di0_2binNC + strata(id), cclong)
  mod0_3bin <- clogit(case ~ exp_di0_3bin + exp_di0_3binNC + strata(id), cclong)
  mod0_4bin <- clogit(case ~ exp_di0_4bin + exp_di0_4binNC + strata(id), cclong)
  mod0_5bin <- clogit(case ~ exp_di0_5bin + exp_di0_5binNC + strata(id), cclong)
  
  mod_contbin <- clogit(case ~  exp_di1_2bin + exp_di2_3bin +
                          exp_di3_4bin + exp_di4_5bin +
                           exp_di1_2binNC + exp_di2_3binNC +
                          exp_di3_4binNC + exp_di4_5binNC + strata(id), cclong)
  
  c3 <- data.frame(rbind(#summary(mod0_1bin)$conf.int,
                         summary(mod0_2bin)$conf.int,
                         summary(mod0_3bin)$conf.int,
                         summary(mod0_4bin)$conf.int,
                         summary(mod0_5bin)$conf.int))
  c3$name <- c("0-2 km", "0-2 km", "0-3 km", "0-3 km", 
               "0-4 km", "0-4 km", "0-5 km", "0-5 km")
  c3$model <- 3
  c3$Estimate <- rep(c("Observed", "Negative Control"), 4)
  
  c4 <- data.frame(summary(mod_contbin)$conf.int)
  c4$name <- rep(c("1-2 km", "2-3 km", "3-4 km", "4-5 km"),2)
  c4$model <- 4
  c4$Estimate <- c(rep("Observed", 4), rep("Negative Control", 4))
 
   # Quantiles
  if (model5 == T){
    breaks <- quantile(cclong$exp_di0_5[cclong$exp_di0_5 >0])
    breaks <- c(0, breaks)
    
    cclong$diQ <- cut(cclong$exp_di0_5, breaks = breaks, include.lowest = T)
    cclong$diQNC<- cut(cclong$exp_di0_5NC, breaks = breaks, include.lowest = T)
    
    modQ <- clogit(case ~ diQ + diQNC + strata(id), cclong)
    
    c5 <- data.frame(summary(modQ)$conf.int)
    c5$name <- c("Q1", "Q2", "Q3", "Q4", "Q1", "Q2", "Q3", "Q4")
    c5$model <- 5
    c5$Estimate <- c(rep("Observed", 4), rep("Negative Control", 4))
    
    coefs <- rbind(c1,c2,c3,c4,c5)}
  else {
    coefs <- rbind(c1,c2,c3,c4)}
  return(coefs)
}

waldtest <- function(df1, df2, model) {
  # info about wald test from here: https://stats.stackexchange.com/questions/143264/citation-for-statistical-test-for-difference-between-two-odds-ratios
  # Note: can adjust which distance is being considered, as needed. default = 0-5km
  if (model == "C") {
    mod0_5_1 <- clogit(case ~ exp_di0_5 + exp_di0_5NC + strata(id), df1)
    mod0_5_2 <- clogit(case ~ exp_di0_5 + exp_di0_5NC + strata(id), df2)}
  
  if (model == "B") {
    mod0_5_1 <- clogit(case ~ exp_di0_5bin + exp_di0_5binNC + strata(id), df1)
    mod0_5_2 <- clogit(case ~ exp_di0_5bin + exp_di0_5binNC + strata(id), df2)}
  
  out_1 <- summary(mod0_5_1)
  out_2 <- summary(mod0_5_2)
  B1 <- out_1$coefficients[1,1]
  B2 <- out_2$coefficients[1,1]  
  SE1 <- out_1$coefficients[1,3]
  SE2 <- out_2$coefficients[1,3]
  z <- (B2 - B1) / sqrt(SE2^2 + SE1^2)
  pval <- 2*pnorm(-abs(z))
  return(c(z, pval))
}

make_plots <- function(coefs, model_num, xlab, ystart = 0.2, h = 1){
  
  coeft1 <- coefs %>% subset(model == model_num)
  
  p1 <- ggplot(coeft1) +
    geom_point(aes(x = name, y = OR, color = Estimate), size = 2.5, position = position_dodge(width = 0.4)) +
    geom_errorbar(aes(x = name, ymin = lower, ymax = upper, color = Estimate), 
                  width = 0.25, size = 0.6, position = position_dodge(width = 0.4)) +
    geom_hline(aes(yintercept = 1), linetype = "dashed") + scale_x_discrete(limits = rev) +
    theme_minimal() + theme(axis.text=element_text(size=12), axis.title=element_text(size=14)) + 
    scale_color_manual(values = c("#ED0000FF", "#00468BFF")) + 
    #scale_color_paletteer_d(`"ggsci::lanonc_lancet"`) +
    xlab(paste0(xlab)) + ylab("Odds Ratio (OR)") +
    coord_flip() + theme(legend.position = "bottom")
  p1
  
  coeft1$OR95CI <- paste0(round(coeft1$OR, 3), " (", round(coeft1$lower, 3), ", ", round(coeft1$upper, 3), ")")
  coeft1sub <- coeft1 %>% arrange(name) %>% dplyr::select("OR95CI", "name", "Estimate")
  rownames(coeft1sub) <- NULL
  colnames(coeft1sub) <- c("OR (95% CI)", "Distance", "Estimate")
  
  t1 <- ggplot(coeft1sub) +
    geom_text(aes(x = 1, y = Distance, color = Estimate, label = `OR (95% CI)`), 
              position = position_dodge(width = 0.7)) +
    scale_color_manual(values = c("#ED0000FF", "#00468BFF")) +
    # scale_color_paletteer_d(`"ggsci::lanonc_lancet"`) +
    scale_y_discrete(limits = rev) +
    theme_void() + theme(legend.position = "none")
  t1
  #make the table
  
   g <- ggdraw() + 
     draw_plot(p1, x = 0, y = 0, width = 0.7, height = 1) +
     draw_plot(t1, x = 0.65, y = ystart, width = 0.3, height = h)
  
  #g <- plot_grid(p1, t1, align = "h", rel_widths = 0.7, 0.3)
  
  return(g)
}

make_plots_sameAxes <- function(coefs, model_num, xlab, ystart = 0.2, h = 1){
  
  coeft1 <- coefs %>% subset(model == model_num)
  
  p1 <- ggplot(coeft1) + ylim(c(0.7, 1.85)) + 
    geom_point(aes(x = name, y = OR, color = Estimate), size = 2.5, position = position_dodge(width = 0.4)) +
    geom_errorbar(aes(x = name, ymin = lower, ymax = upper, color = Estimate), 
                  width = 0.25, size = 0.6, position = position_dodge(width = 0.4)) +
    geom_hline(aes(yintercept = 1), linetype = "dashed") + scale_x_discrete(limits = rev) +
    theme_minimal() + theme(axis.text=element_text(size=12), axis.title=element_text(size=14)) + 
    scale_color_manual(values = c("#ED0000FF", "#00468BFF")) + 
    #scale_color_paletteer_d(`"ggsci::lanonc_lancet"`) +
    xlab(paste0(xlab)) + ylab("Odds Ratio (OR)") +
    coord_flip() + theme(legend.position = "bottom")
  p1
  
  coeft1$OR95CI <- paste0(round(coeft1$OR, 3), " (", round(coeft1$lower, 3), ", ", round(coeft1$upper, 3), ")")
  coeft1sub <- coeft1 %>% arrange(name) %>% dplyr::select("OR95CI", "name", "Estimate")
  rownames(coeft1sub) <- NULL
  colnames(coeft1sub) <- c("OR (95% CI)", "Distance", "Estimate")
  
  t1 <- ggplot(coeft1sub) +
    geom_text(aes(x = 1, y = Distance, color = Estimate, label = `OR (95% CI)`), 
              position = position_dodge(width = 0.7)) +
    scale_color_manual(values = c("#ED0000FF", "#00468BFF")) +
    # scale_color_paletteer_d(`"ggsci::lanonc_lancet"`) +
    scale_y_discrete(limits = rev) +
    theme_void() + theme(legend.position = "none")
  t1
  #make the table
  
  g <- ggdraw() + 
    draw_plot(p1, x = 0, y = 0, width = 0.7, height = 1) +
    draw_plot(t1, x = 0.65, y = ystart, width = 0.3, height = h)
  
  #g <- plot_grid(p1, t1, align = "h", rel_widths = 0.7, 0.3)
  
  return(g)
}

make_plots_strata <- function(coefs, model_num, xlab, ystart = 0.2, h = 1){
  
  coeft1 <- coefs %>% subset(model == model_num)
  
  p1 <- ggplot(coeft1) +
    geom_point(aes(x = name, y = OR, color = Estimate), size = 2, position = position_dodge(width = 0.4)) +
    geom_errorbar(aes(x = name, ymin = lower, ymax = upper, color = Estimate), 
                  width = 0.25, size = 0.6, position = position_dodge(width = 0.4)) +
    geom_hline(aes(yintercept = 1), linetype = "dashed") +
    theme_minimal() + theme(axis.text=element_text(size=12), axis.title=element_text(size=14)) + 
    scale_x_discrete(limits = rev) +
  # season: scale_color_manual(values = c("#A6611A", "#4DAC26", "#CA0020", "#0571B0")) +
# race/eth scale_color_manual(values = c("#286EB9", "black", "#ED6F57", "#FFB319")) +
# age: scale_color_manual(values = c("#008ECC","#134D9C","#0B2056")) +
 # land cover: scale_color_manual(values = c("#A65B00","#D30000")) + 
    # scale_color_paletteer_d(`"ggsci::lanonc_lancet"`) +
   # sex:  scale_color_manual(values = c("#AA336A","#009dc4")) + 
    xlab(paste0(xlab)) + ylab("Odds Ratio (OR)") +
    coord_flip() + theme(legend.position = "bottom")
  p1
  
  coeft1$OR95CI <- paste0(round(coeft1$OR, 3), " (", round(coeft1$lower, 3), ", ", round(coeft1$upper, 3), ")")
  coeft1sub <- coeft1 %>% arrange(name) %>% dplyr::select("OR95CI", "name", "Estimate")
  rownames(coeft1sub) <- NULL
  colnames(coeft1sub) <- c("OR (95% CI)", "Distance", "Estimate")
  
  t1 <- ggplot(coeft1sub) +
    geom_text(aes(x = 1, y = Distance, color = Estimate, label = `OR (95% CI)`), 
              position = position_dodge(width = 0.7), size = 4) +
 # season: scale_color_manual(values = c("#A6611A", "#4DAC26", "#CA0020", "#0571B0")) +
# race/eth: scale_color_manual(values = c("#286EB9", "black", "#ED6F57", "#FFB319")) +
# age:  scale_color_manual(values = c("#008ECC","#134D9C","#0B2056")) +
# land cover: scale_color_manual(values = c("#A65B00","#D30000")) + 
# sex: scale_color_manual(values = c("#AA336A","#009dc4")) + 
    #scale_color_paletteer_d(`"ggsci::lanonc_lancet"`) +
    scale_y_discrete(limits = rev) +
    theme_void() + theme(legend.position = "none")
  t1

   g <- ggdraw() + 
     draw_plot(p1, x = 0, y = 0, width = 0.7, height = 1) +
     draw_plot(t1, x = 0.65, y = ystart, width = 0.3, height = h)
  
  return(g)
}

make_plots_sensitivity <- function(coefs, model_num, xlab){
  
  coeft1 <- coefs %>% subset(model == model_num)
  
  p1 <- ggplot(coeft1) +
    geom_point(aes(x = name, y = OR, color = Estimate), size = 2, position = position_dodge(width = 0.4)) +
    geom_errorbar(aes(x = name, ymin = lower, ymax = upper, color = Estimate), 
                  width = 0.25, size = 0.6, position = position_dodge(width = 0.4)) +
    geom_hline(aes(yintercept = 1), linetype = "dashed") +
    theme_minimal() + theme(axis.text=element_text(size=12), axis.title=element_text(size=14)) + 
    scale_x_discrete(limits = rev) +
    scale_color_manual(values = c("#7851A9", "#73c91c","#658D1B","#ff595e", "#1982c4", "#0033A0", "#b394b9")) + 
    xlab(paste0(xlab)) + ylab("Odds Ratio (OR)") +
    coord_flip() + theme(legend.position = "bottom")
  p1
  
  coeft1$OR95CI <- paste0(round(coeft1$OR, 3), " (", round(coeft1$lower, 3), ", ", round(coeft1$upper, 3), ")")
  coeft1sub <- coeft1 %>% arrange(name) %>% dplyr::select("OR95CI", "name", "Estimate")
  rownames(coeft1sub) <- NULL
  colnames(coeft1sub) <- c("OR (95% CI)", "Distance", "Estimate")
  
  t1 <- ggplot(coeft1sub) +
    geom_text(aes(x = 1, y = Distance, color = Estimate, label = `OR (95% CI)`), 
              position = position_dodge(width = 0.7), size = 4) +
    scale_color_manual(values = c("#7851A9", "#73c91c","#658D1B","#ff595e", "#1982c4", "#0033A0", "#b394b9")) + 
    scale_y_discrete(limits = rev) +
    theme_void() + theme(legend.position = "none")
  t1
  
  g <- ggdraw() + 
    draw_plot(p1, x = 0, y = 0, width = 0.7, height = 1) +
    draw_plot(t1, x = 0.65, y = 0.15, width = 0.3, height = 0.85)
  
  return(g)
}

make_plots_obsvOnly <- function(coefs, model_num, xlab, ystart = 0.2, h = 1){
  
  coeft1 <- coefs %>% subset(model == model_num)
  coeft1 <- coeft1 %>% subset(Estimate == "Observed")
  
  p1 <- ggplot(coeft1) +
    geom_point(aes(x = name, y = OR), size = 3, position = position_dodge(width = 0.4)) +
    geom_errorbar(aes(x = name, ymin = lower, ymax = upper), 
                  width = 0.2, size = 0.75, position = position_dodge(width = 0.4)) +
    geom_hline(aes(yintercept = 1), linetype = "dashed") +
    theme_minimal() + theme(axis.text=element_text(size=12), axis.title=element_text(size=14)) + 
    scale_x_discrete(limits = rev) +
    scale_color_paletteer_d(`"ggsci::lanonc_lancet"`) +
    xlab(paste0(xlab)) + ylab("Odds Ratio (OR)") +
    coord_flip() + theme(legend.position = "bottom")
  p1
  
  coeft1$OR95CI <- paste0(round(coeft1$OR, 3), " (", round(coeft1$lower, 3), ", ", round(coeft1$upper, 3), ")")
  coeft1sub <- coeft1 %>% arrange(name) %>% dplyr::select("OR95CI", "name", "Estimate")
  rownames(coeft1sub) <- NULL
  colnames(coeft1sub) <- c("OR (95% CI)", "Distance", "Estimate")
  
  t1 <- ggplot(coeft1sub) +
    geom_text(aes(x = 1, y = Distance, label = `OR (95% CI)`), 
              position = position_dodge(width = 0.7), size = 4.5) +
    scale_color_paletteer_d(`"ggsci::lanonc_lancet"`) +
    scale_y_discrete(limits = rev) +
    theme_void() + theme(legend.position = "none")
  t1

  g <- ggdraw() + 
    draw_plot(p1, x = 0, y = 0, width = 0.7, height = 1) +
    draw_plot(t1, x = 0.65, y = ystart, width = 0.3, height = h)
  
  return(g)
}

make_plots_obsvOnly_sameAxes <- function(coefs, model_num, xlab, ystart = 0.2, h = 1){
  
  coeft1 <- coefs %>% subset(model == model_num)
  coeft1 <- coeft1 %>% subset(Estimate == "Observed")
  
  p1 <- ggplot(coeft1) + ylim(c(0.80,1.6)) + 
    geom_point(aes(x = name, y = OR), size = 3, position = position_dodge(width = 0.4)) +
    geom_errorbar(aes(x = name, ymin = lower, ymax = upper), 
                  width = 0.2, size = 0.75, position = position_dodge(width = 0.4)) +
    geom_hline(aes(yintercept = 1), linetype = "dashed") +
    theme_minimal() + theme(axis.text=element_text(size=12), axis.title=element_text(size=14)) + 
    scale_x_discrete(limits = rev) +
    scale_color_paletteer_d(`"ggsci::lanonc_lancet"`) +
    xlab(paste0(xlab)) + ylab("Odds Ratio (OR)") +
    coord_flip() + theme(legend.position = "bottom")
  p1
  
  coeft1$OR95CI <- paste0(round(coeft1$OR, 3), " (", round(coeft1$lower, 3), ", ", round(coeft1$upper, 3), ")")
  coeft1sub <- coeft1 %>% arrange(name) %>% dplyr::select("OR95CI", "name", "Estimate")
  rownames(coeft1sub) <- NULL
  colnames(coeft1sub) <- c("OR (95% CI)", "Distance", "Estimate")
  
  t1 <- ggplot(coeft1sub) +
    geom_text(aes(x = 1, y = Distance, label = `OR (95% CI)`), 
              position = position_dodge(width = 0.7), size = 4.5) +
    scale_color_paletteer_d(`"ggsci::lanonc_lancet"`) +
    scale_y_discrete(limits = rev) +
    theme_void() + theme(legend.position = "none")
  t1
  
  g <- ggdraw() + 
    draw_plot(p1, x = 0, y = 0, width = 0.7, height = 1) +
    draw_plot(t1, x = 0.65, y = ystart, width = 0.3, height = h)
  
  return(g)
}

