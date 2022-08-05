---
title: "Excerno analysis"
author: "Audrey Mitchell, Marco Ruiz, Soua Yang, Jaime Davila"
date: "8/3/2022"
output:
  html_document: 
    keep_md: yes
editor_options:
  chunk_output_type: console
---




# Set up

We start our analysis by loading our libraries and functions from other files


```r
library(NMF)
library(MutationalPatterns)
library(tidyverse)
library(ggrepel)
library(ROCR)
library(kableExtra)
library(gridExtra)
library(broom)
library(plotROC)
source("bayes_function.R")
source("get_classification_df.R")
source("simulation.R")
```

We also define the function `vector_to_matrix` which allows to switch between representations of mutational profiles, in particular to use `plot_96_profile` from `MutationalPatterns`


```r
vector_to_matrix <- function(mutations.vector){
  tab <- table(mutations.vector)/sum(table(mutations.vector))
  #transform table to a data frame and assign column names to mutations and frequencies
  mutations.df <- data.frame(tab)
  colnames(mutations.df) <- c("mutations", "frequencies")
  
  #defines a function that joins two data frames and replaces NA with 0
  left_join_NA <- function(x, y, by) {
    left_join(x = x, y = y, by = by) %>% 
      mutate_each(funs(replace(., which(is.na(.)), 0)))
  }
  
  #converts vector of 96 mutation types to a data frame
  mutation.types.df <- data.frame(mutations)
  #joins together sample probabilities with data frame of mutation types so that all 96 mutation types are included
  complete.mutations.df <- left_join_NA(mutation.types.df, mutations.df, by = "mutations")
  
  #converts data frame of mutation types and frequencies to a data frame and assigns rownames
  mutations.matrix <- as.matrix(complete.mutations.df$frequencies)
  rownames(mutations.matrix) = mutations
  
  return(mutations.matrix)
}
```

We start by loading the FFPE mutational signatures and the SBS mutational signatures from COSMIC:


```r
#load FFPE signature as a matrix
ffpe_signature <- as.matrix(load_old_mutational_matrix("supplied_results/ffpe.signature.txt"))

#create vector of the 96 mutation types
mutations <- rownames(ffpe_signature) 

#Load COSMIC signatures (version 3)
cosmic_signatures <- get_known_signatures(muttype = "snv")

#Create vector of the signature names
signature_names <- colnames(cosmic_signatures)
```

And we create a list with all the COSMIC mutational signatures


```r
#Create list of signature matrices with assigned rownames
signature_matrix_list <- vector ("list", length(signature_names))
for (i in 1:length(signature_names)){
  sig <- cosmic_signatures[,i]
  sig_mat <- as.matrix(sig)
  rownames(sig_mat) <- mutations
  colnames(sig_mat) <- signature_names[i]
  signature_matrix_list[[i]] <- sig_mat
}
```



# Initial testing

We created two simulated datasets as follows:

* We combine 500 mutations from SBS4 and 500 mutations from FFPE

* We combine 500 mutations from SBS6 and 500 mutations from FFPE


## First dataset: SBS4

Let's create and visualize our first dataset by simulating mutations from SBS4


```r
set.seed(500)
signature_4_mat <- as.matrix(cosmic_signatures[,"SBS4"])
rownames(signature_4_mat) <- mutations

sig4_500 <- create_signature_sample_vector(signature_4_mat, 500)
sig4_profile <- vector_to_matrix(sig4_500)

plot_96_profile(sig4_profile)+
  labs(title = "Signature 4 Simulated Sample", x = "Context")+
  theme(plot.title = element_text(hjust = 0.5))
```

![](paper_analysis_files/figure-html/unnamed-chunk-5-1.png)<!-- -->

Then, let's add the FFPE noise and visualize the resulting mutational profile


```r
addffpe_500 <- create_signature_sample_vector(ffpe_signature, 500)
sig4_ffpe <- vector_to_matrix(c(sig4_500, addffpe_500))

plot_96_profile(sig4_ffpe)+
  labs(title = "Simulated Sample with FFPE Noise", x = "Context")+
  theme(plot.title = element_text(hjust = 0.5))
```

![](paper_analysis_files/figure-html/unnamed-chunk-6-1.png)<!-- -->

Now, let's apply our method on the simulated dataset:


```r
sig4_ffpe_tbl <- get_classification_df(
  list(sig4_500, addffpe_500),
  c("SBS4", "FFPE"), 
  list(signature_4_mat, ffpe_signature)) %>%
  as_tibble()
```

And calculate the confusion matrix


```r
(conf_mat <- table(sig4_ffpe_tbl$truth, sig4_ffpe_tbl$classify))
```

```
##       
##        FFPE SBS4
##   FFPE  500    0
##   SBS4   52  448
```

Our method is very good at identifying mutations that originated from the FFPE process, although some mutations originated from SBS4 gets missclassified as FFPE. Let's calculate the performance metrics:



```r
(sig4_sensitivity <- (conf_mat[2,2] / (conf_mat[2,2]+ conf_mat[2,1])))
```

```
## [1] 0.896
```

```r
(sig4_specificity <- (conf_mat[1,1] / (conf_mat[1,1]+ conf_mat[1,2])))
```

```
## [1] 1
```

## Second dataset: SBS6

Let's repeat the same experiment, this time using SBS6 as our baseline signature


```r
signature_6_mat <- as.matrix(cosmic_signatures[,"SBS6"])
rownames(signature_6_mat) <- mutations
sig6_500 <- create_signature_sample_vector(signature_6_mat, 500)

sig6_ffpe_tbl <- get_classification_df(
  list(sig6_500, addffpe_500),
  c("SBS6", "FFPE"),
  list(signature_6_mat, ffpe_signature))%>%
  as_tibble()
```

And let's calculate our confusion matrix and our performance characteristics


```r
(conf_mat <- table(sig6_ffpe_tbl$truth, sig6_ffpe_tbl$classify))
```

```
##       
##        FFPE SBS6
##   FFPE  365  135
##   SBS6   83  417
```

```r
(sig6_sensitivity <- (conf_mat[2,2] / (conf_mat[2,2]+ conf_mat[2,1])))
```

```
## [1] 0.834
```

```r
(sig6_specificity <- (conf_mat[1,1] / (conf_mat[1,1]+ conf_mat[1,2])))
```

```
## [1] 0.73
```

It seems this time our method did not perform as well. Let's create mutational profile plots for:

* The original profile (SBS6)
* The original profile (SBS6) plus the FFPE noise
* The profile obtained after filtering FFPE using `excerno`



```r
# Original profile
sig6_profile <- vector_to_matrix(sig6_500)
sig6_profile_plot <- plot_96_profile(sig6_profile)+
  labs(title = "Signature 6 Simulated Sample", x = "Context")+
  theme(plot.title = element_text(hjust = 0.5))
# Profile + FFPE noise
sig6_ffpe <- vector_to_matrix(c(sig6_500, addffpe_500))
sig6_ffpe_plot <- plot_96_profile(sig6_ffpe)+
  labs(title = "Simulated Sample with FFPE Noise", x = "Context")+
  theme(plot.title = element_text(hjust = 0.5))
# Profile after excerno
sig6_ffpe_filtered <- sig6_ffpe_tbl %>%
  filter(classify != "FFPE")
sig6_recon_matrix <- vector_to_matrix(sig6_ffpe_filtered$mutations)

sig6_recon_plot <- plot_96_profile(sig6_recon_matrix)+
  labs(title = "Signature 6 Reconstructed", x = "Context")+
  theme(plot.title = element_text(hjust = 0.5))

(sig6_grid <- grid.arrange(sig6_profile_plot, sig6_ffpe_plot, sig6_recon_plot))
```

![](paper_analysis_files/figure-html/unnamed-chunk-12-1.png)<!-- -->

```
## TableGrob (3 x 1) "arrange": 3 grobs
##   z     cells    name           grob
## 1 1 (1-1,1-1) arrange gtable[layout]
## 2 2 (2-2,1-1) arrange gtable[layout]
## 3 3 (3-3,1-1) arrange gtable[layout]
```

```r
ggsave("paper_figures/sup1.png", sig6_grid, width = 175, height = 243, units = "mm")
```

It seems filtered SBS6 dataset differs from the original signature, particularly in C>T mutations preceded by an A or a T 


## ROC curve for both datasets

We would like to combine the information from both datasets (SBS4 and SBS6) in order to plot them together.


```r
sig4_tbl <- sig4_ffpe_tbl %>%
  select(truth, SBS4) %>%
  rename (bayes = SBS4) %>%
  mutate (dataset = "SBS4",
          truth = if_else (truth=="SBS4",1, 0))

sig6_tbl <- sig6_ffpe_tbl %>%
  select(truth, SBS6) %>%
  rename (bayes = SBS6) %>%
  mutate (dataset = "SBS6",
          truth = if_else (truth=="SBS6",1, 0))

allsim_tbl <- rbind (sig4_tbl, sig6_tbl)
```

And we generate the ROC curve for both datasets and calculate their AUC



```r
scale = seq (0, 1, by=0.2)

(roc_plot <- 
  ggplot(allsim_tbl, aes(d=truth, m=bayes, color=dataset))+
  geom_roc(cutoffs.at = c(0.5), cutoff.labels = c(0.5)) +
  labs(x="False Positive Rate",
       y="True Positive Rate",
       color="Signature")+
  scale_x_continuous(labels=scale, breaks=scale)+
  scale_y_continuous(labels=scale, breaks=scale))
```

![](paper_analysis_files/figure-html/unnamed-chunk-14-1.png)<!-- -->

```r
ggsave("paper_figures/fig2a.png", roc_plot, width = 150, height = 90, units = "mm")

calc_auc(roc_plot)
```

```
##   PANEL group dataset      AUC
## 1     1     1    SBS4 0.959482
## 2     1     2    SBS6 0.842970
```

# When is `excerno` a good classifier?

We hypothesize that the performance of our classifier depends on how different the base mutational signature is from the FFPE signature. 

## Initial Setup

A straightforward way of calculating the similarity across signatures the cosine similarity. So we first generate the cosine similarity of all the COSMIC signatures against the FFPE signature:


```r
ffpe_cosmic_cos_sims <- vector(length=length(signature_names))
for (i in 1:length(signature_names)){
  cos_sim <- cos_sim_matrix(signature_matrix_list[[i]], ffpe_signature)
  ffpe_cosmic_cos_sims[i] <- cos_sim[1,1]
}
```

Second, we will create a function that mixes mutations from FFPE and another given mutational signature and calls `excerno`


```r
test_mix_mutations <- function (cosmic_idx, num_cosmic, num_ffpe) {

  # Create baseline mutations
  baseline_sig <- signature_matrix_list[[cosmic_idx]]
  baseline_name <- signature_names[cosmic_idx]
  baseline_mut_vec <- create_signature_sample_vector(baseline_sig, num_cosmic)
  
  # Create FFPE mutations
  ffpe_mut_vec <- create_signature_sample_vector(ffpe_signature, num_ffpe)

  # Call excerno classification and format as a tibble
  get_classification_df (list(baseline_mut_vec, ffpe_mut_vec), 
                          c(baseline_name,"FFPE"), 
                          list(baseline_sig, ffpe_signature)) %>% 
    as_tibble()
}
```

And we try `test_mix_mutations` and plot the results to check it is working:


```r
set.seed(12345)
test_tbl <- test_mix_mutations(1,500,500)

test_tbl %>%
  select (mutations, truth, SBS1, FFPE) %>%
  pivot_longer(c(SBS1, FFPE)) %>%
  ggplot(aes(x=truth, y=value, fill=name)) +
    geom_boxplot() 
```

![](paper_analysis_files/figure-html/unnamed-chunk-17-1.png)<!-- -->

Now let's create a function that calculates the performance metrics (sensitivity, specificity, accuracy) on the output of `test_mix_mutations`


```r
calc_metrics_mix <- function (cosmic_idx, num_cosmic, num_ffpe) {
  
  baseline_name <- signature_names[cosmic_idx]
  test_tbl <- test_mix_mutations(cosmic_idx,num_cosmic,num_ffpe)
  
  tp <- test_tbl %>%
    filter(truth == baseline_name & classify == baseline_name) %>%
    count() %>%
    pull(n)
  
  fp <- test_tbl %>%
    filter(truth != baseline_name & classify == baseline_name) %>%
    count() %>%
    pull(n)
    
  fn <- test_tbl %>%
    filter(truth == baseline_name & classify == "FFPE") %>%
    count() %>%
    pull(n)
  
  tn <- test_tbl %>%
    filter(truth == "FFPE" & classify == "FFPE") %>%
    count() %>%
    pull(n)
  
  tibble(
    Sensitivity = tp / (tp +fn),
    Specificity = tn / (tn + fp),
    Accuracy = (tp +tn) /(tp+tn+fp+fn))
}
```

And we test our method:


```r
set.seed(12345)

calc_metrics_mix(1,500,500)
```

```
## # A tibble: 1 x 3
##   Sensitivity Specificity Accuracy
##         <dbl>       <dbl>    <dbl>
## 1       0.956       0.922    0.939
```

```r
calc_metrics_mix(4,500,500)
```

```
## # A tibble: 1 x 3
##   Sensitivity Specificity Accuracy
##         <dbl>       <dbl>    <dbl>
## 1       0.906           1    0.953
```

```r
calc_metrics_mix(6,500,500)
```

```
## # A tibble: 1 x 3
##   Sensitivity Specificity Accuracy
##         <dbl>       <dbl>    <dbl>
## 1        0.81       0.736    0.773
```

We create a function that calculates the average metrics over a fixed number of simulations `num_iter` and test it


```r
calc_metrics_mix_num <- function (cosmic_idx, num_iter, num_cosmic, num_ffpe) {

  baseline_name <- signature_names[cosmic_idx]
  
  1:num_iter %>% 
    map_dfr(~calc_metrics_mix(cosmic_idx,num_cosmic,num_ffpe))  %>%
    pivot_longer(1:3) %>% 
    group_by(name) %>%
    summarize(mean = mean(value)) %>%
    mutate(sbs = baseline_name)
}

calc_metrics_mix_num(1,10,500,500)
```

```
## # A tibble: 3 x 3
##   name         mean sbs  
##   <chr>       <dbl> <chr>
## 1 Accuracy    0.93  SBS1 
## 2 Sensitivity 0.948 SBS1 
## 3 Specificity 0.912 SBS1
```

```r
calc_metrics_mix_num(4,10,500,500)
```

```
## # A tibble: 3 x 3
##   name         mean sbs  
##   <chr>       <dbl> <chr>
## 1 Accuracy    0.955 SBS4 
## 2 Sensitivity 0.909 SBS4 
## 3 Specificity 1     SBS4
```

```r
calc_metrics_mix_num(6,10,500,500)
```

```
## # A tibble: 3 x 3
##   name         mean sbs  
##   <chr>       <dbl> <chr>
## 1 Accuracy    0.784 SBS6 
## 2 Sensitivity 0.826 SBS6 
## 3 Specificity 0.742 SBS6
```

## Simulations across all COSMIC signatures

We are ready to generate the data over all of COSMIC mutational signatures, using 10 iterations and 500 mutations for FFPE and 500 mutations from the baseline COSMIC signature



```r
sbs_num = length(signature_names)
base_num = 500
ffpe_num = 500
reps = 10 

set.seed(123456)
sim_tbl <- 1:sbs_num %>%
  map_dfr(~calc_metrics_mix_num(.x, reps, base_num, ffpe_num))
```

First we summarize the metrics across all the COSMIC mutational signatures


```r
sim_tbl <- sim_tbl %>%
  mutate (name = factor (name, 
                         levels = c("Sensitivity", "Specificity","Accuracy"))) 
sim_tbl %>%
  ggplot(aes(mean))+
  geom_histogram(binwidth = 0.025)+
  facet_grid(.~name) +
  labs(x="Metric across all COSMIC signatures")
```

![](paper_analysis_files/figure-html/unnamed-chunk-22-1.png)<!-- -->

```r
ggsave("paper_figures/supp2.png", dpi = 300, width = 175, height = 75, units = "mm")

summary_tbl <- sim_tbl %>%
  group_by(name) %>%
  summarize (median = median (mean),
           IQR = IQR(mean)) %>%
  mutate_if(is.numeric, format, digits=2)

summary_tbl %>%
  kbl() %>%
  kable_styling()
```

<table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> name </th>
   <th style="text-align:left;"> median </th>
   <th style="text-align:left;"> IQR </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Sensitivity </td>
   <td style="text-align:left;"> 0.89 </td>
   <td style="text-align:left;"> 0.111 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Specificity </td>
   <td style="text-align:left;"> 1.00 </td>
   <td style="text-align:left;"> 0.091 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Accuracy </td>
   <td style="text-align:left;"> 0.94 </td>
   <td style="text-align:left;"> 0.120 </td>
  </tr>
</tbody>
</table>


## Performance characteristics vs similarity against FFPE signature

Now we explore the relationship of the performance metrics against the similarity between the baseline signature and FFPE using a figure


```r
sim_tbl <- sim_tbl %>%
  mutate (ffpe_cos = rep(ffpe_cosmic_cos_sims, each=3))

sim_tbl %>%
  ggplot(aes(x=ffpe_cos, mean)) +
  geom_point()+
  geom_smooth(method="lm", formula = "y~x")+
  facet_grid(~name)
```

![](paper_analysis_files/figure-html/unnamed-chunk-23-1.png)<!-- -->

We define a function to represent the mutational signatures that fall below a threshold for a particular metric.


```r
plot_metrics_excerno <- function (tbl, metric, threshold) {

  # Subset the table to the specific metric  
  new_tbl <- filter(tbl, name==metric) 
  # Identify signature that fall below the threshold
  outliers <- filter( new_tbl, (mean < threshold) | sbs == "SBS4" | sbs == "SBS6")

  ggplot(new_tbl, aes(x = ffpe_cos, y = mean))+
    geom_point()+
    geom_smooth(method = "lm", se=FALSE)+
    geom_point(data = outliers, color = "red")+
    geom_hline(yintercept=threshold, linetype="dashed", color = "red")+
    geom_label_repel(data=outliers, 
                   aes(label = sbs),
                   max.overlaps = 10, 
                   label.size = .1)+
    labs(x = "Cosine similarity against FFPE signature", 
         y = metric )
}
```

And we plot all the metrics using a threshold of 80%


```r
plot_metrics_excerno(sim_tbl, "Accuracy", 0.8) +
  scale_x_continuous(labels=scale, breaks=scale)
```

![](paper_analysis_files/figure-html/unnamed-chunk-25-1.png)<!-- -->

```r
ggsave("paper_figures/fig2b.png", dpi = 300, width = 175, height = 108, units = "mm")

plot_metrics_excerno(sim_tbl, "Sensitivity", 0.8)+
    scale_x_continuous(labels=scale, breaks=scale)
```

![](paper_analysis_files/figure-html/unnamed-chunk-25-2.png)<!-- -->

```r
ggsave("paper_figures/supp3.png", dpi = 300, width = 175, height = 108, units = "mm")

plot_metrics_excerno(sim_tbl, "Specificity", 0.8)+
  scale_x_continuous(labels=scale, breaks=scale)
```

![](paper_analysis_files/figure-html/unnamed-chunk-25-3.png)<!-- -->

```r
ggsave("paper_figures/supp4.png", dpi = 300, width = 175, height = 108, units = "mm")
```

And we create auxiliar function for getting the coefficients of the linear model and their fit


```r
print_lm_coefs <- function (tbl, metric) {
  new_tbl <-  filter (tbl, name==metric)
  lm_metrics<- lm(mean ~ ffpe_cos, data = new_tbl)

  tidy(lm_metrics) %>%
    mutate_if(is.numeric, format, digits=2)%>%
    kbl() %>%
    kable_styling()
}

print_lm_fit <- function (tbl, metric) {
  new_tbl <-  filter (tbl, name==metric)
  lm_metrics<- lm(mean ~ ffpe_cos, data = new_tbl)

  glance(lm_metrics) %>%
    mutate_if(is.numeric, format, digits=2)%>%
    kbl() %>%
    kable_styling()
}
```

And print the information of the linear model for 

* Accuracy


```r
print_lm_coefs(sim_tbl, "Accuracy")
```

<table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> term </th>
   <th style="text-align:left;"> estimate </th>
   <th style="text-align:left;"> std.error </th>
   <th style="text-align:left;"> statistic </th>
   <th style="text-align:left;"> p.value </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> (Intercept) </td>
   <td style="text-align:left;"> 0.99 </td>
   <td style="text-align:left;"> 0.0053 </td>
   <td style="text-align:left;"> 186 </td>
   <td style="text-align:left;"> 3.3e-82 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ffpe_cos </td>
   <td style="text-align:left;"> -0.37 </td>
   <td style="text-align:left;"> 0.0149 </td>
   <td style="text-align:left;"> -24 </td>
   <td style="text-align:left;"> 2.9e-32 </td>
  </tr>
</tbody>
</table>

```r
print_lm_fit(sim_tbl, "Accuracy")
```

<table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> r.squared </th>
   <th style="text-align:left;"> adj.r.squared </th>
   <th style="text-align:left;"> sigma </th>
   <th style="text-align:left;"> statistic </th>
   <th style="text-align:left;"> p.value </th>
   <th style="text-align:left;"> df </th>
   <th style="text-align:left;"> logLik </th>
   <th style="text-align:left;"> AIC </th>
   <th style="text-align:left;"> BIC </th>
   <th style="text-align:left;"> deviance </th>
   <th style="text-align:left;"> df.residual </th>
   <th style="text-align:left;"> nobs </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> 0.91 </td>
   <td style="text-align:left;"> 0.91 </td>
   <td style="text-align:left;"> 0.028 </td>
   <td style="text-align:left;"> 600 </td>
   <td style="text-align:left;"> 2.9e-32 </td>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:left;"> 129 </td>
   <td style="text-align:left;"> -253 </td>
   <td style="text-align:left;"> -247 </td>
   <td style="text-align:left;"> 0.047 </td>
   <td style="text-align:left;"> 58 </td>
   <td style="text-align:left;"> 60 </td>
  </tr>
</tbody>
</table>

* Sensitivity


```r
print_lm_coefs(sim_tbl, "Sensitivity")
```

<table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> term </th>
   <th style="text-align:left;"> estimate </th>
   <th style="text-align:left;"> std.error </th>
   <th style="text-align:left;"> statistic </th>
   <th style="text-align:left;"> p.value </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> (Intercept) </td>
   <td style="text-align:left;"> 0.94 </td>
   <td style="text-align:left;"> 0.010 </td>
   <td style="text-align:left;"> 92 </td>
   <td style="text-align:left;"> 1.2e-64 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ffpe_cos </td>
   <td style="text-align:left;"> -0.31 </td>
   <td style="text-align:left;"> 0.029 </td>
   <td style="text-align:left;"> -11 </td>
   <td style="text-align:left;"> 2.5e-15 </td>
  </tr>
</tbody>
</table>

```r
print_lm_fit(sim_tbl, "Sensitivity")
```

<table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> r.squared </th>
   <th style="text-align:left;"> adj.r.squared </th>
   <th style="text-align:left;"> sigma </th>
   <th style="text-align:left;"> statistic </th>
   <th style="text-align:left;"> p.value </th>
   <th style="text-align:left;"> df </th>
   <th style="text-align:left;"> logLik </th>
   <th style="text-align:left;"> AIC </th>
   <th style="text-align:left;"> BIC </th>
   <th style="text-align:left;"> deviance </th>
   <th style="text-align:left;"> df.residual </th>
   <th style="text-align:left;"> nobs </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> 0.66 </td>
   <td style="text-align:left;"> 0.66 </td>
   <td style="text-align:left;"> 0.055 </td>
   <td style="text-align:left;"> 114 </td>
   <td style="text-align:left;"> 2.5e-15 </td>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:left;"> 90 </td>
   <td style="text-align:left;"> -175 </td>
   <td style="text-align:left;"> -168 </td>
   <td style="text-align:left;"> 0.17 </td>
   <td style="text-align:left;"> 58 </td>
   <td style="text-align:left;"> 60 </td>
  </tr>
</tbody>
</table>

* Specificity

```r
print_lm_coefs(sim_tbl, "Specificity")
```

<table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> term </th>
   <th style="text-align:left;"> estimate </th>
   <th style="text-align:left;"> std.error </th>
   <th style="text-align:left;"> statistic </th>
   <th style="text-align:left;"> p.value </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> (Intercept) </td>
   <td style="text-align:left;"> 1.04 </td>
   <td style="text-align:left;"> 0.014 </td>
   <td style="text-align:left;"> 75 </td>
   <td style="text-align:left;"> 2.1e-59 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ffpe_cos </td>
   <td style="text-align:left;"> -0.42 </td>
   <td style="text-align:left;"> 0.039 </td>
   <td style="text-align:left;"> -11 </td>
   <td style="text-align:left;"> 9.4e-16 </td>
  </tr>
</tbody>
</table>

```r
print_lm_fit(sim_tbl, "Specificity")
```

<table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> r.squared </th>
   <th style="text-align:left;"> adj.r.squared </th>
   <th style="text-align:left;"> sigma </th>
   <th style="text-align:left;"> statistic </th>
   <th style="text-align:left;"> p.value </th>
   <th style="text-align:left;"> df </th>
   <th style="text-align:left;"> logLik </th>
   <th style="text-align:left;"> AIC </th>
   <th style="text-align:left;"> BIC </th>
   <th style="text-align:left;"> deviance </th>
   <th style="text-align:left;"> df.residual </th>
   <th style="text-align:left;"> nobs </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> 0.67 </td>
   <td style="text-align:left;"> 0.67 </td>
   <td style="text-align:left;"> 0.074 </td>
   <td style="text-align:left;"> 120 </td>
   <td style="text-align:left;"> 9.4e-16 </td>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:left;"> 72 </td>
   <td style="text-align:left;"> -138 </td>
   <td style="text-align:left;"> -132 </td>
   <td style="text-align:left;"> 0.32 </td>
   <td style="text-align:left;"> 58 </td>
   <td style="text-align:left;"> 60 </td>
  </tr>
</tbody>
</table>

# What is the influence of FFPE percentage in the performance of `excerno`?

We are interested in the performance characteristics of excerno under different percentages of FFPE signature. 

## An initial simulation

So we will mimick our original scenario scenario of SBS4 and SBS6, however this time we will introduce varying proportions going from 5% to 100% every 5%.

Let's start with SBS4 



```r
set.seed(654321)
pct_ffpe = seq (10,90, by=5)
metrics_order = c("Sensitivity","Specificity", "Accuracy")

sbs4_dilution_tbl <- pct_ffpe %>%
  map_dfr(~calc_metrics_mix_num(4,10, (100-.x)*10, .x*10)) %>%
  mutate(pct_ffpe = rep(pct_ffpe, each=3),
         name = factor (name, levels = metrics_order))
```

And plot it


```r
ggplot(sbs4_dilution_tbl, 
       aes(x=pct_ffpe, y=mean)) +
  geom_point()+
  geom_smooth(se=FALSE)+
  facet_grid(name~.) + 
  geom_hline(yintercept=0.8, linetype="dashed",color = "red")
```

![](paper_analysis_files/figure-html/unnamed-chunk-31-1.png)<!-- -->

Seems the metrics stay above 80%, except for specificity at 10% which slightly dips below 80%

Let's do the same analysis with SBS6:


```r
sbs6_dilution_tbl <- pct_ffpe %>%
  map_dfr(~calc_metrics_mix_num(6,10, (100-.x)*10, .x*10)) %>%
  mutate(pct_ffpe = rep(pct_ffpe, each=3),
         name = factor (name, levels = metrics_order))

ggplot(sbs6_dilution_tbl, 
       aes(x=pct_ffpe, y=mean)) +
  geom_point()+
  geom_smooth(se=FALSE)+
  facet_grid(name~.) + 
  geom_hline(yintercept=0.8, linetype="dashed",color = "red")
```

![](paper_analysis_files/figure-html/unnamed-chunk-32-1.png)<!-- -->

Here we see how sensitivity dips below our threshold for all ffpe percentages *above* 60% ffpe and specificity drops below 80% for all the % ffpe *below* 50%. Accuracy stays mostly above 80% and it is best when there is only one class of mutations (either FFPE or a single SBS)

Let's combine both datasets in a single graph by using color and focus on sensitivity and specificity.


```r
dilution_tbl <- rbind (sbs4_dilution_tbl, sbs6_dilution_tbl)

scale_x = seq (10,90, by=10)
scale_y = seq (0.2, 1, by=0.2)

dilution_tbl %>%
  filter (name!="Accuracy") %>%
  ggplot(aes(x=pct_ffpe, y=mean, color=sbs)) +
    geom_point()+
    geom_smooth(se=FALSE)+
    facet_grid(name~.) + 
    geom_hline(yintercept=0.8, linetype="dashed",color = "red")+
    labs(x= "Percent of FFPE signature",
         y="",
         color = "Signature")+
    scale_x_continuous(breaks = scale_x, labels = scale_x)+
    scale_y_continuous(breaks = scale_y, labels = scale_y)
```

![](paper_analysis_files/figure-html/unnamed-chunk-33-1.png)<!-- -->

```r
ggsave("paper_figures/fig3a.png", dpi = 300, width = 150, height = 90, units = "mm")
```

Notice the trend is that for increasing percent of FFPE signature, the sensitivity decreases, while the specificity increases.

## A general approach

We are interested in assessing the trends in performance metrics across *all* of the signatures in COSMIC. In order to do that we need to create an auxiliary function `calc_dillution_metrics` that generates a tibble with this information for a particular mutational signature. We also create a function `plot_dillution_metrics` that plots this information


```r
calc_dillution_metrics <- function (ffpe_grid, cosmic_idx) {
  metrics_order = c("Sensitivity","Specificity", "Accuracy")
  sig_name=factor(signature_names[cosmic_idx], signature_names)

  ffpe_grid %>%
    map_dfr(~calc_metrics_mix_num(cosmic_idx,10, (100-.x)*10, .x*10)) %>%
    mutate(pct_ffpe = rep(ffpe_grid, each=3),
         sbs = sig_name,
         metric = factor (name, levels = metrics_order)) %>%
    select(-c(name))
}

plot_dillution_metrics <- function (sim_tbl) {
  scale_x = seq (10,90, by=10)
  scale_y = seq (0.2, 1, by=0.2)

  sim_tbl %>%
  ggplot(aes(x=pct_ffpe, y=mean, color=sbs)) +
    geom_point()+
    geom_smooth(se=FALSE)+
    facet_grid(metric~.) + 
    geom_hline(yintercept=0.8, linetype="dashed",color = "red")+
    labs(x= "Percent of FFPE signature",
         y="",
         color = "Signature")+
    scale_x_continuous(breaks = scale_x, labels = scale_x)+
    scale_y_continuous(breaks = scale_y, labels = scale_y)
}
```

We test our functions for SBS4 and SBS6:


```r
calc_dillution_metrics(pct_ffpe, 4) %>%
  plot_dillution_metrics()
```

![](paper_analysis_files/figure-html/unnamed-chunk-35-1.png)<!-- -->

```r
calc_dillution_metrics(pct_ffpe, 6) %>%
  plot_dillution_metrics()
```

![](paper_analysis_files/figure-html/unnamed-chunk-35-2.png)<!-- -->

Now we are ready to calculate the dillution metrics across all COSMIC mutational signatures (this takes 20min to run)


```r
set.seed(123456)
# Full dataset run
sim_dillution_tbl <- 1:sbs_num %>%
  map_dfr(~calc_dillution_metrics(pct_ffpe, .x))
```

And we plot the sensitivity and specificity across all levels of FFPE contamination:


```r
scale_x = seq (10,90, by=10)
scale_y = seq (0.2, 1, by=0.2)

sim_dillution_tbl %>%
  filter(metric != "Accuracy") %>%
  ggplot (aes(x=pct_ffpe, y=mean, color=metric)) +
    geom_point (alpha=0.3) +
    geom_smooth(method="lm", se=FALSE)+
    facet_grid(metric~.)+
    geom_hline(yintercept=0.8, linetype="dashed",color = "red")+
    scale_x_continuous(breaks = scale_x, labels = scale_x)+
    scale_y_continuous(breaks = scale_y, labels = scale_y)+
    labs(x= "Percent of FFPE signature",
         y="")+
    theme(legend.position = "none")
```

![](paper_analysis_files/figure-html/unnamed-chunk-37-1.png)<!-- -->

```r
ggsave("paper_figures/supp5.png", dpi = 300, width = 150, height = 90, units = "mm")
```




```r
sessionInfo()
```

```
## R version 4.1.0 (2021-05-18)
## Platform: x86_64-w64-mingw32/x64 (64-bit)
## Running under: Windows 10 x64 (build 19043)
## 
## Matrix products: default
## 
## locale:
## [1] LC_COLLATE=English_United States.1252 
## [2] LC_CTYPE=English_United States.1252   
## [3] LC_MONETARY=English_United States.1252
## [4] LC_NUMERIC=C                          
## [5] LC_TIME=English_United States.1252    
## 
## attached base packages:
## [1] stats4    parallel  stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] plotROC_2.3.0            broom_0.7.10             gridExtra_2.3           
##  [4] kableExtra_1.3.4         ROCR_1.0-11              ggrepel_0.9.1           
##  [7] forcats_0.5.1            stringr_1.4.0            dplyr_1.0.8             
## [10] purrr_0.3.4              readr_2.1.1              tidyr_1.2.0             
## [13] tibble_3.1.6             ggplot2_3.3.5            tidyverse_1.3.1         
## [16] MutationalPatterns_3.2.0 GenomicRanges_1.44.0     GenomeInfoDb_1.28.4     
## [19] IRanges_2.26.0           S4Vectors_0.30.2         NMF_0.23.0              
## [22] Biobase_2.52.0           BiocGenerics_0.38.0      cluster_2.1.2           
## [25] rngtools_1.5.2           pkgmaker_0.32.2          registry_0.5-1          
## 
## loaded via a namespace (and not attached):
##  [1] nlme_3.1-153           fs_1.5.0               bitops_1.0-7          
##  [4] lubridate_1.8.0        webshot_0.5.2          doParallel_1.0.16     
##  [7] RColorBrewer_1.1-3     httr_1.4.2             tools_4.1.0           
## [10] backports_1.4.0        bslib_0.3.1            utf8_1.2.2            
## [13] R6_2.5.1               mgcv_1.8-38            DBI_1.1.1             
## [16] colorspace_2.0-3       withr_2.5.0            tidyselect_1.1.2      
## [19] compiler_4.1.0         rvest_1.0.2            cli_3.2.0             
## [22] xml2_1.3.3             labeling_0.4.2         sass_0.4.0            
## [25] scales_1.2.0           systemfonts_1.0.3      digest_0.6.29         
## [28] svglite_2.0.0          rmarkdown_2.11         XVector_0.32.0        
## [31] pkgconfig_2.0.3        htmltools_0.5.2        highr_0.9             
## [34] dbplyr_2.1.1           fastmap_1.1.0          rlang_1.0.2           
## [37] readxl_1.3.1           rstudioapi_0.13        farver_2.1.0          
## [40] jquerylib_0.1.4        generics_0.1.2         jsonlite_1.8.0        
## [43] RCurl_1.98-1.5         magrittr_2.0.1         GenomeInfoDbData_1.2.6
## [46] Matrix_1.3-4           Rcpp_1.0.8.3           munsell_0.5.0         
## [49] fansi_1.0.3            lifecycle_1.0.1        stringi_1.7.6         
## [52] yaml_2.3.5             ggalluvial_0.12.3      zlibbioc_1.38.0       
## [55] plyr_1.8.7             grid_4.1.0             crayon_1.5.1          
## [58] lattice_0.20-45        splines_4.1.0          haven_2.4.3           
## [61] hms_1.1.1              knitr_1.36             pillar_1.7.0          
## [64] reshape2_1.4.4         codetools_0.2-18       reprex_2.0.1          
## [67] glue_1.6.2             evaluate_0.14          modelr_0.1.8          
## [70] vctrs_0.4.1            tzdb_0.2.0             foreach_1.5.1         
## [73] cellranger_1.1.0       gtable_0.3.0           assertthat_0.2.1      
## [76] xfun_0.23              gridBase_0.4-7         xtable_1.8-4          
## [79] pracma_2.3.6           viridisLite_0.4.0      iterators_1.0.13      
## [82] ellipsis_0.3.2
```

