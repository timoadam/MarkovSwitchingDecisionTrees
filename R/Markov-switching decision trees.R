## Load required packages
library(rpart)
library(dplyr)
library(rpart.plot)
library(vip)

## Load data
data_raw <- read.csv("v5.csv")
data <- data_raw %>% filter(play_type %in% c("pass", "run")) #, posteam %in% c("NE", "SEA", "HOU", "CLE"))
data <- data[c("posteam", "game_id", "play_type", "ydstogo", "down", "score_differential", "shotgun", "qtr", "posteam_type",
               "game_date")]

## Create factors
data$down <- as.factor(data$down)
data$shotgun <- as.factor(data$shotgun)
data$qtr <- as.factor(data$qtr)

# unique(data$posteam)
data$game_date <- as.Date(data$game_date)
data_train <- data %>% filter(game_date < as.Date("2018-09-09"))
data_test <- data %>% filter(game_date >= as.Date("2018-09-09"))

## Fit baseline model
mod0 <- rpart(play_type ~ ydstogo + as.factor(down) + score_differential + as.factor(shotgun) + as.factor(qtr), data = data, method = "class", control = list(minbucket = 100, cp = 0.001))
rpart.plot(mod0, type = 5)

## fit_msdt: function that fits Markov-switching decision trees
## Inptuts:
##  data: a data frame
##  N: number of states
##  weights0: initial weights
##  stat: initial distribution = stationary distribution?
##  max_iter: maximum number of iterations
##  conv_tol: converegence tolerance
## Output:
##  a list containing various objects
fit_msdt <- function(data, N = 2, weights0 = NULL, stat = FALSE, max_iter = 1000, conv_tol = 1e-03) { 
  # Formula for decision trees
  formula = play_type ~ ydstogo + down + score_differential + shotgun + qtr
  # Initialisation
  delta = NULL
  gamma = matrix(c(0.95, 0.05, 0.05, 0.95), ncol = 2)
  if(is.null(weights0)) {
    weights0 = matrix(NA, nrow = nrow(data), ncol = N)
    weights0[, 1] = runif(nrow(weights0), 0, 1)
    weights0[, 2] = 1 - weights0[, 1]
  }
  mod = list()
  term = FALSE
  old = 0
  allprobs = matrix(NA, nrow(data), N)
  lalpha = lbeta = matrix(NA, N, nrow(data))
  for(i in 1:N) {
    mod[[i]] = rpart(formula, weights = weights0[, i], data = data, method = "class", control = list(minbucket = 100, cp = 0.001))
    foo1 <- predict(mod[[i]], data = data, type = "prob")
    foo2 <- rep(1, times = nrow(data))
    foo2[data$play_type == "run"] = 2
    foo2[foo2 == 1] = foo1[which(foo2 == 1), 1]
    foo2[foo2 == 2] = foo1[which(foo2 == 2), 2]
    allprobs[, i] = foo2
  }
  allprobs = ifelse(!is.na(allprobs), allprobs, 1)
  # Loop
  while(term == FALSE) {
    for(i in 1:max_iter) {
      delta_next = delta
      gamma_next = gamma
      if(is.null(delta)) {
        delta = solve(t(diag(N) - gamma + 1), rep(1, N))
      }
      # E-step
      foo = delta * allprobs[1,]
      sumfoo = sum(foo)
      lscale = log(sumfoo)
      foo = foo / sumfoo
      lalpha[, 1] = log(foo) + lscale
      for(j in 2:nrow(data)) {
        if(data$game_id[j] != data$game_id[j - 1]) {
          foo = foo %*% matrix(rep(delta, 2), ncol = 2, byrow = TRUE) * allprobs[j,] # Checken, ob das mit dem delta so Sinn ergibt
          sumfoo = sum(foo)
          lscale = lscale + log(sumfoo)
          foo = foo / sumfoo
          lalpha[, j] = log(foo) + lscale
        } else {
          foo = foo %*% gamma * allprobs[j,]
          sumfoo = sum(foo)
          lscale = lscale + log(sumfoo)
          foo = foo / sumfoo
          lalpha[, j] = log(foo) + lscale
        }
      }
      weights = NULL
      foo = rep(1 / N, N)
      lbeta[, nrow(data)] = rep(0, N)
      foo = foo / sum(foo)
      lscale = log(N)
      for(j in (nrow(data) - 1):1) {
        if(data$game_id[j + 1] != data$game_id[j]) {
          foo = matrix(rep(delta, 2), ncol = 2, byrow = TRUE) %*% (allprobs[j + 1,] * foo) # Checken, ob das mit dem delta so Sinn ergibt
          lbeta[, j] = log(foo) + lscale
          sumfoo = sum(foo)
          foo = foo / sumfoo
          lscale = lscale + log(sumfoo)
        } else {
          foo = gamma %*% (allprobs[j + 1,] * foo)
          lbeta[, j] = log(foo) + lscale
          sumfoo = sum(foo)
          foo = foo / sumfoo
          lscale = lscale + log(sumfoo)
        }    
      }                   
      lallprobs = log(allprobs)
      llh = max(lalpha[, nrow(data)]) + log(sum(exp(lalpha[, nrow(data)] - max(lalpha[, nrow(data)]))))
      weights = matrix(NA, N, nrow(data))
      for(j in 1:nrow(data)) {
        weights[,j] = exp(lalpha[,j] + lbeta[, j] - llh)
      }
      # M step
      for(j in 1:N) {
        for(k in 1:N) {
          gamma_next[j, k] = gamma[j, k] * sum(exp(lalpha[j, 1:(nrow(data) - 1)] + lallprobs[2:nrow(data), k] + lbeta[k, 2:nrow(data)] - llh))
        }
      }
      gamma_next = gamma_next / apply(gamma_next, 1, sum)
      if(stat == TRUE) {
        delta_next = solve(t(diag(N) - gamma_next + 1), rep(1, N))
      }else{
        delta_next = exp(lalpha[, 1] + lbeta[, 1] - llh)
        delta_next = delta_next / sum(delta_next)
      }
      for(j in 1:N){
        mod[[j]] = rpart(formula, weights = weights[j,], data = data, method = "class", control = list(minbucket = 100, cp = 0.001))
        foo1 <- predict(mod[[j]], data = data, type = "prob")
        foo2 <- rep(1, times = nrow(data))
        foo2[data$play_type == "run"] = 2
        foo2[foo2 == 1] = foo1[which(foo2 == 1), 1]
        foo2[foo2 == 2] = foo1[which(foo2 == 2), 2]
        allprobs[, j] = foo2
      }
      # Check convergence
      cat("Iteration = ", i, ", log-likelihood = ", round(llh, 3), "\r", sep = "")
      conv_crit = abs(llh - old)
      if(conv_crit < conv_tol | i == max_iter) {
        if(i == max_iter) {
          print(paste("No convergence within", max_iter, "iterations"))
        }else{
          print(paste("Convergence after", i, "iterations, log-likelihood =", round(llh, 3)))
        }
        term = TRUE
        break
      }  
      delta = delta_next
      gamma = gamma_next
      old = llh
    }
  }
  # Return output
  return(list(mod = mod, delta = delta_next, gamma = gamma_next, llh = llh, state_probs = weights, N = N))
}

## viterbi: function that decodes the states
## Inputs:
##  mod: a fitted Markov-switching decision tree
##  data: a data frame
## Output:
##  the most likely state sequence
viterbi <- function(mod, data) {
  n = nrow(data)
  allprobs = matrix(NA, n, mod$N)
  for(j in 1:mod$N) {
    foo1 <- predict(mod$mod[[j]], data = data, type = "prob")
    foo2 <- rep(1, times = n)
    foo2[data$play_type == "run"] = 2
    foo2[foo2 == 1] = foo1[which(foo2 == 1), 1]
    foo2[foo2 == 2] = foo1[which(foo2 == 2), 2]
    allprobs[, j] = foo2
  }
  xi = matrix(0, n, mod$N)
  foo = mod$delta * allprobs[1,]
  xi[1,] = foo / sum(foo)
  for(i in 2:n){
    foo = apply(xi[i - 1,] * mod$gamma, 2, max) * allprobs[i,]
    xi[i,] = foo / sum(foo)
  }
  iv = numeric(n)
  iv[n] = which.max(xi[n,])
  for(i in (n - 1):1){
    iv[i] <- which.max(mod$gamma[,iv[i + 1]] * xi[i,])
  }
  # Return most likely state sequence
  return(iv)
}

## Fit models with different initial weights
teams_considered <- c("NE", "CLE") # c("NE", "SEA", "HOU", "CLE")
nruns <- 100
llks <- list()
mods_all <- list()

for(j in 1:length(teams_considered)){
  mods_team <- list()
  foo_team <- rep(NA, nruns)
  sub <- filter(data_train, posteam == teams_considered[j])
  for(i in 1:nruns){
    set.seed(i)
    print(paste0("Team ", j, "/", length(teams_considered), "; ","Run = ", i, "/", nruns))
    
    weights0 <- matrix(NA, nrow = nrow(sub), ncol = 2)
    weights0[, 1] <- runif(nrow(weights0), 0, 1)
    weights0[, 2] <- 1 - weights0[, 1]
    mods_team[[i]] <- fit_msdt(data = sub, weights0 = weights0)
    foo_team[i] <- mods_team[[i]]$llh
  }
  mods_all[[j]] <- mods_team
  llks[[j]] <-  foo_team
}

load("multiple_starting_values.RData")

## Get seed that results in the model with the highest log-likelihood
llk_max_idx <- sapply(llks, which.max)
mod_NE <- mods_all[[1]][llk_max_idx[1]][[1]]
# mod_SEA <- mods_all[[2]][llk_max_idx[2]][[1]]
# mod_HOU <- mods_all[[3]][llk_max_idx[3]][[1]]
mod_CLE <- mods_all[[2]][llk_max_idx[2]][[1]]

data_NE <- filter(data_train, posteam == "NE")
# data_SEA <- filter(data_train, posteam == "SEA")
# data_HOU <- filter(data_train, posteam == "HOU")
data_CLE <- filter(data_train, posteam == "CLE")

# load("final_models.RData")

## Get estimated initial/state transition probabilities
mod_NE$delta
mod_NE$gamma

## Decode the states
states <- viterbi(mod = mod_NE, data = data_NE)

## Plot results four teams
pdf(file = "msdt_results.pdf",   # The directory you want to save the file in
    width = 12, # The width of the plot in inches
    height = 12) # The height of the plot in inches

par(mfrow = c(2, 4), oma = c(0, 2, 2, 0))
## NE
rpart.plot(mod_NE$mod[[1]], type = 5, main = "State 1")
rpart.plot(mod_NE$mod[[2]], type = 5, main = "State 2")
#mtext("State 1 State 2", outer = TRUE, cex = 1.333)
#mtext("NE", side = 1, line = 0, cex = 1.333, las = 3, outer = TRUE)
## SEA
rpart.plot(mod_SEA$mod[[1]], type = 5, main = "State 1")
rpart.plot(mod_SEA$mod[[2]], type = 5, main = "State 2")
#mtext("State 1 State 2", outer = TRUE, cex = 1.333)
#mtext("SEA", side = 2, line = 0, cex = 1.333, las = 3, outer = TRUE)
## HOU
rpart.plot(mod_HOU$mod[[1]], type = 5, main = "State 1")
rpart.plot(mod_HOU$mod[[2]], type = 5, main = "State 2")
## CLE
rpart.plot(mod_CLE$mod[[1]], type = 5, main = "State 1")
rpart.plot(mod_CLE$mod[[2]], type = 5, main = "State 2")

dev.off()


## Plot results two teams
pdf(file = "msdt_results_2teams.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 8) # The height of the plot in inches

par(mfrow = c(2, 2), oma = c(0, 2, 2, 0))
## NE
rpart.plot(mod_SEA$mod[[1]], type = 5, main = "State 1", cex = 0.45)
rpart.plot(mod_SEA$mod[[2]], type = 5, main = "State 2", cex = 0.45)
#mtext("State 1 State 2", outer = TRUE, cex = 1.333)
#mtext("NE", side = 1, line = 0, cex = 1.333, las = 3, outer = TRUE)
## CLE
rpart.plot(mod_HOU$mod[[1]], type = 5, main = "State 1", cex = 0.45)
rpart.plot(mod_HOU$mod[[2]], type = 5, main = "State 2", cex = 0.45)

dev.off()

# Example time series -----------------------------------------------------

data_raw %>% filter(posteam == "SEA") %>% 
  group_by(game_id) %>% 
  mutate(nr_points = max(total_home_score) + max(total_away_score)) %>% 
  summarise(score = mean(nr_points)) %>% 
  arrange(-score)

example_playcalls <- data_raw %>% filter(posteam == "SEA", game_id == 2015111511, play_type %in% c("pass", "run")) %>%
                       pull(play_type)

example.data <- data.frame(play = 1:length(example_playcalls), type = example_playcalls)

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
timeline_plot <- ggplot(example.data, aes(x = play, y = 0, col = factor(type))) + 
  scale_color_manual(name = "", values = cbbPalette, labels = c("Pass", "Run")) + geom_hline(yintercept=0, 
                                                       color = "black", size = 0.3) +
  geom_point(aes(y = 0), size = 4) + theme_classic() + ylim(-0.1, 0.1) +
  theme(axis.line.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x =element_blank(),
        axis.ticks.x =element_blank(),
        axis.line.x =element_blank(),
        legend.box = "horizontal",
        legend.direction = "horizontal", 
        #legend.position = "bottom"
        legend.position = c(0.5, 0.35)
  ) +
  geom_text(data = example.data %>% filter(play %in% c(1, 10, 20, 30, 40)), aes(x = play, y = -0.01, 
                                     label = play),
            size = 2.5, vjust = 0.5, color = 'black') 

timeline_plot

cowplot::save_plot("example_ts.pdf", timeline_plot)


# Viterbi -----------------------------------------------------------------

NE_ids <- data_train %>% filter(posteam == "NE") %>% pull(game_id) %>% unique

nr_switches <- rep(NA, length(NE_ids))
for(i in 1:length(NE_ids)){
  cur_vit <- viterbi(mod_NE, data_train %>% filter(posteam == "NE", game_id == NE_ids[i]))
  nr_switches[i] <- diff(cur_vit) %>% unique %>% length
}
which(nr_switches == max(nr_switches))
switch_ids <- NE_ids[which(nr_switches == max(nr_switches))]
data_raw %>% filter(posteam == "NE") %>% 
  filter(game_id %in% switch_ids) %>% 
  group_by(game_id) %>% 
  mutate(nr_points = max(total_home_score) + max(total_away_score)) %>% 
  summarise(avg_score = mean(nr_points)) %>% 
  arrange(-avg_score)
data_train %>% filter(posteam == "NE", game_id == 2012121613)

data_viterbi <- data_train %>% filter(posteam == "NE", game_id == 2012121613)
data_viterbi$viterbi <- viterbi(mod_NE, data_viterbi)

ggplot(data_viterbi, aes(x = 1:nrow(data_viterbi), y = play_type, colour = factor(viterbi))) + 
  geom_point()


# Variable importance plots -----------------------------------------------

library(ggplot2)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# SEA
df <- data.frame(imp = mod_SEA$mod[[1]]$variable.importance)
df2_state1 <- df %>% 
  tibble::rownames_to_column() %>% 
  dplyr::rename("variable" = rowname) %>% 
  dplyr::arrange(imp) %>%
  dplyr::mutate(variable = forcats::fct_inorder(variable)) %>% 
  mutate(state = "State 1")
# vip_plot_NE_state1 <- ggplot(df2) +
#   geom_col(aes(x = variable, y = imp),
#            col = "black", show.legend = FALSE) +
#   coord_flip() +
#   scale_fill_grey() +
#   theme_minimal()

df <- data.frame(imp = mod_SEA$mod[[2]]$variable.importance)
df2_state2 <- df %>% 
  tibble::rownames_to_column() %>% 
  dplyr::rename("variable" = rowname) %>% 
  dplyr::arrange(imp) %>%
  dplyr::mutate(variable = forcats::fct_inorder(variable)) %>% 
  mutate(state = "State 2")
df_plot_NE <- bind_rows(df2_state1, df2_state2)
vip_plot_NE <- ggplot(df_plot_NE) +
  geom_col(aes(x = variable, y = imp, fill = forcats::fct_relevel(state, "State 2")), position = "dodge") +
  #geom_col(aes(x = variable, y = imp, fill = state), position = "dodge") +
  coord_flip() +
  ylab("variable importance") +
  scale_fill_manual(name = "", values = cbPalette, breaks = c("State 1", "State 2")) +
  theme_minimal() + ggtitle("Seattle Seahawks") +
  theme(legend.position = "none")
vip_plot_NE

# HOU
df <- data.frame(imp = mod_HOU$mod[[2]]$variable.importance)
df2_state1 <- df %>% 
  tibble::rownames_to_column() %>% 
  dplyr::rename("variable" = rowname) %>% 
  dplyr::arrange(imp) %>%
  dplyr::mutate(variable = forcats::fct_inorder(variable)) %>% 
  mutate(state = "State 1")

df <- data.frame(imp = mod_HOU$mod[[1]]$variable.importance)
df2_state2 <- df %>% 
  tibble::rownames_to_column() %>% 
  dplyr::rename("variable" = rowname) %>% 
  dplyr::arrange(imp) %>%
  dplyr::mutate(variable = forcats::fct_inorder(variable)) %>% 
  mutate(state = "State 2")
df_plot_CLE <- bind_rows(df2_state1, df2_state2)
vip_plot_CLE <- ggplot(df_plot_CLE) +
  geom_col(aes(x = variable, y = imp, fill = forcats::fct_relevel(state, "State 2")), position = "dodge") +
  #geom_col(aes(x = variable, y = imp, fill = state), position = "dodge") +
  coord_flip() +
  ylab("variable importance") +
  xlab("") +
  scale_fill_manual(name = "", values = cbPalette, breaks = c("State 1", "State 2")) +
  theme_minimal() + ggtitle("Houston Texans") +
  theme(legend.position = c(-.6, .6))
vip_plot_CLE

vip_plot <- cowplot::plot_grid(vip_plot_NE, vip_plot_CLE, ncol = 2)
cowplot::save_plot("vip_plot.pdf", vip_plot)

# Prediction --------------------------------------------------------------

## Predictions with standard decision tree
mod_NE_standard <- rpart(play_type ~ ydstogo + down + score_differential + shotgun + qtr, 
                         data = filter(data_train, posteam == "NE"), 
                         method = "class", control = list(minbucket = 100, cp = 0.001))
mod_CLE_standard <- rpart(play_type ~ ydstogo + down + score_differential + shotgun + qtr, 
                          data = filter(data_train, posteam == "CLE"), 
                          method = "class", control = list(minbucket = 100, cp = 0.001))
mod_SEA_standard <- rpart(play_type ~ ydstogo + down + score_differential + shotgun + qtr, 
                          data = filter(data_train, posteam == "SEA"), 
                          method = "class", control = list(minbucket = 100, cp = 0.001))
mod_HOU_standard <- rpart(play_type ~ ydstogo + down + score_differential + shotgun + qtr, 
                          data = filter(data_train, posteam == "HOU"), 
                          method = "class", control = list(minbucket = 100, cp = 0.001))

predictions_NE_standard <- predict(mod_NE_standard, newdata = filter(data_test, posteam == "NE"), type = "class")
predictions_CLE_standard <- predict(mod_CLE_standard, newdata = filter(data_test, posteam == "CLE"), type = "class")
predictions_SEA_standard <- predict(mod_SEA_standard, newdata = filter(data_test, posteam == "SEA"), type = "class")
predictions_HOU_standard <- predict(mod_HOU_standard, newdata = filter(data_test, posteam == "HOU"), type = "class")

## Predictions with MS-DT

# Create list with data frames for each new time series
data_test_NE <- filter(data_test, posteam == "NE") %>% split(f = .[["game_id"]])
data_test_SEA <- filter(data_test, posteam == "SEA") %>% split(f = .[["game_id"]])
data_test_CLE <- filter(data_test, posteam == "CLE") %>% split(f = .[["game_id"]])
data_test_HOU <- filter(data_test, posteam == "HOU") %>% split(f = .[["game_id"]])

msdt_forecast <- function(xf, h = 1, dataset, mod, first_obs = FALSE) {
  n <- nrow(dataset)
  true_response <- dataset$play_type
  nxf <- length(xf)
  dxf <- matrix(0, nrow = h, ncol = nxf)
  
  prob_state1 <- predict(mod$mod[[1]], newdata = dataset[1,], type = "prob")
  prob_state1 <- prob_state1[dimnames(prob_state1)[[2]] == "pass"]
  prob_state2 <- predict(mod$mod[[2]], newdata = dataset[1,], type = "prob")
  prob_state2 <- prob_state2[dimnames(prob_state2)[[2]] == "pass"]
  # 2 = run; 1 = pass
  
  if(first_obs){
    prob_pass_t1 <- mod$delta[1] * prob_state1 + mod$delta[2] * prob_state2
    dxf[1, ] <- c(prob_pass_t1, 1 - prob_pass_t1)
  }
  else{
    ## check t=1
    if(true_response[1] == "run"){
      foo <- mod$delta * c(1 - prob_state1, 1 - prob_state2)
    }
    else{
      foo <- mod$delta * c(prob_state1, prob_state2)
    }
    #foo <- mod$delta * c(prob_state1, prob_state2) # dpois(x[1], mod$lambda)
    sumfoo <- sum(foo)
    lscale <- log(sumfoo)
    foo <- foo / sumfoo
    for (i in 2:(n - 1)) {
      prob_state1 <- predict(mod$mod[[1]], newdata = dataset[i,], type = "prob")
      prob_state1_pass <- prob_state1[dimnames(prob_state1)[[2]] == "pass"]
      prob_state1_run <- prob_state1[dimnames(prob_state1)[[2]] == "run"]
      prob_state2 <- predict(mod$mod[[2]], newdata = dataset[i,], type = "prob")
      prob_state2_pass <- prob_state2[dimnames(prob_state2)[[2]] == "pass"]
      prob_state2_run <- prob_state2[dimnames(prob_state2)[[2]] == "run"]
      if(true_response[i] == "run"){
        foo <- foo %*% mod$gamma * c(prob_state1_run, prob_state2_run)
      }
      else{
        foo <- foo %*% mod$gamma * c(prob_state1_pass, prob_state2_pass)
      }
      
      #foo <- foo %*% mod$gamma * dpois(x[i], mod$lambda)
      sumfoo <- sum(foo)
      lscale <- lscale + log (sumfoo)
      foo <- foo / sumfoo
    }
    
    for (i in 1:h) {
      foo <- foo %*% mod$gamma
      for (j in 1:mod$N){
        prob_state_j <- predict(mod$mod[[j]], newdata = dataset[n, ], type = "prob")
        prob_state_j_pass <- prob_state_j[dimnames(prob_state_j)[[2]] == "pass"]
        prob_state_j_run <- prob_state_j[dimnames(prob_state_j)[[2]] == "run"]
        
        dxf[i ,] <- dxf[i ,] + foo[j] * c(prob_state_j_pass, prob_state_j_run) # dpois(xf, mod$lambda[j])
      }
    }
  }
  return(dxf)
}


for(i in 1:length(data_test_NE)) {
  cur_match <- data_test_NE[[i]]
  all_preds <- rep(NA, nrow(cur_match))
  # t = 1
  cur_predictions_01 <- msdt_forecast(xf = 0:1, h = 1, dataset = cur_match[1, ], mod = mod_NE, first_obs = TRUE)
  cur_predictions <- apply(cur_predictions_01, 1, which.max)
  cur_predictions[cur_predictions == 1] <- "pass"
  cur_predictions[cur_predictions == 2] <- "run"
  all_preds[1] <- cur_predictions
  for(t in 2:nrow(cur_match)){
    cur_predictions_01 <- msdt_forecast(xf = 0:1, h = 1, dataset = cur_match[1:t, ], mod = mod_NE, first_obs = FALSE)
    cur_predictions <- apply(cur_predictions_01, 1, which.max)
    cur_predictions[cur_predictions == 1] <- "pass"
    cur_predictions[cur_predictions == 2] <- "run"
    all_preds[t] <- cur_predictions
  }
  data_test_NE[[i]]$pred_MSDT <- all_preds
}

for(i in 1:length(data_test_CLE)) {
  cur_match <- data_test_CLE[[i]]
  all_preds <- rep(NA, nrow(cur_match))
  # t = 1
  cur_predictions_01 <- msdt_forecast(xf = 0:1, h = 1, dataset = cur_match[1, ], mod = mod_CLE, first_obs = TRUE)
  cur_predictions <- apply(cur_predictions_01, 1, which.max)
  cur_predictions[cur_predictions == 1] <- "pass"
  cur_predictions[cur_predictions == 2] <- "run"
  all_preds[1] <- cur_predictions
  for(t in 2:nrow(cur_match)){
    cur_predictions_01 <- msdt_forecast(xf = 0:1, h = 1, dataset = cur_match[1:t, ], mod = mod_CLE, first_obs = FALSE)
    cur_predictions <- apply(cur_predictions_01, 1, which.max)
    cur_predictions[cur_predictions == 1] <- "pass"
    cur_predictions[cur_predictions == 2] <- "run"
    all_preds[t] <- cur_predictions
  }
  data_test_CLE[[i]]$pred_MSDT <- all_preds
}

data_test_NE_df <- bind_rows(data_test_NE)
data_test_CLE_df <- bind_rows(data_test_CLE)
mean(data_test_NE_df$pred_MSDT != data_test_NE_df$play_type)
mean(predictions_NE_standard != data_test_NE_df$play_type)

mean(data_test_CLE_df$pred_MSDT != data_test_CLE_df$play_type)
mean(predictions_CLE_standard != data_test_CLE_df$play_type)

for(i in 1:length(data_test_SEA)) {
  cur_match <- data_test_SEA[[i]]
  all_preds <- rep(NA, nrow(cur_match))
  # t = 1
  cur_predictions_01 <- msdt_forecast(xf = 0:1, h = 1, dataset = cur_match[1, ], mod = mod_SEA, first_obs = TRUE)
  cur_predictions <- apply(cur_predictions_01, 1, which.max)
  cur_predictions[cur_predictions == 1] <- "pass"
  cur_predictions[cur_predictions == 2] <- "run"
  all_preds[1] <- cur_predictions
  for(t in 2:nrow(cur_match)){
    cur_predictions_01 <- msdt_forecast(xf = 0:1, h = 1, dataset = cur_match[1:t, ], mod = mod_SEA, first_obs = FALSE)
    cur_predictions <- apply(cur_predictions_01, 1, which.max)
    cur_predictions[cur_predictions == 1] <- "pass"
    cur_predictions[cur_predictions == 2] <- "run"
    all_preds[t] <- cur_predictions
  }
  data_test_SEA[[i]]$pred_MSDT <- all_preds
}
data_test_SEA_df <- bind_rows(data_test_SEA)

mean(data_test_SEA_df$pred_MSDT != data_test_SEA_df$play_type)
mean(predictions_SEA_standard != data_test_SEA_df$play_type)

for(i in 1:length(data_test_HOU)) {
  cur_match <- data_test_HOU[[i]]
  all_preds <- rep(NA, nrow(cur_match))
  # t = 1
  cur_predictions_01 <- msdt_forecast(xf = 0:1, h = 1, dataset = cur_match[1, ], mod = mod_HOU, first_obs = TRUE)
  cur_predictions <- apply(cur_predictions_01, 1, which.max)
  cur_predictions[cur_predictions == 1] <- "pass"
  cur_predictions[cur_predictions == 2] <- "run"
  all_preds[1] <- cur_predictions
  for(t in 2:nrow(cur_match)){
    cur_predictions_01 <- msdt_forecast(xf = 0:1, h = 1, dataset = cur_match[1:t, ], mod = mod_HOU, first_obs = FALSE)
    cur_predictions <- apply(cur_predictions_01, 1, which.max)
    cur_predictions[cur_predictions == 1] <- "pass"
    cur_predictions[cur_predictions == 2] <- "run"
    all_preds[t] <- cur_predictions
  }
  data_test_HOU[[i]]$pred_MSDT <- all_preds
}
data_test_HOU_df <- bind_rows(data_test_HOU)

mean(data_test_HOU_df$pred_MSDT != data_test_HOU_df$play_type)
mean(predictions_HOU_standard != data_test_HOU_df$play_type)

## Predictions of the first match ------------------------------------------
mod_NE_standard <- rpart(play_type ~ ydstogo + down + score_differential + shotgun + qtr, 
                         data = filter(data_train, posteam == "NE"), 
                         method = "class", control = list(minbucket = 100, cp = 0.001))
mod_CLE_standard <- rpart(play_type ~ ydstogo + down + score_differential + shotgun + qtr, 
                          data = filter(data_train, posteam == "CLE"), 
                          method = "class", control = list(minbucket = 100, cp = 0.001))

data_test_NE <- filter(data_test, posteam == "NE", game_id == 2018090905)
data_test_CLE <- filter(data_test, posteam == "CLE", game_id == 2018090901)

predictions_NE_standard <- predict(mod_NE_standard, newdata = data_test_NE, type = "class")
predictions_CLE_standard <- predict(mod_CLE_standard, newdata = data_test_CLE, type = "class")

mean(data_test_NE$play_type != predictions_NE_standard)
mean(data_test_CLE$play_type != predictions_CLE_standard)

## Predictions with MS-DT for NE
cur_match <- data_test_NE
all_preds <- rep(NA, nrow(cur_match))
# t = 1
cur_predictions_01 <- msdt_forecast(xf = 0:1, h = 1, dataset = cur_match[1, ], mod = mod_NE, first_obs = TRUE)
cur_predictions <- apply(cur_predictions_01, 1, which.max)
cur_predictions[cur_predictions == 1] <- "pass"
cur_predictions[cur_predictions == 2] <- "run"
all_preds[1] <- cur_predictions
for(t in 2:nrow(cur_match)){
  cur_predictions_01 <- msdt_forecast(xf = 0:1, h = 1, dataset = cur_match[1:t, ], mod = mod_NE, first_obs = FALSE)
  cur_predictions <- apply(cur_predictions_01, 1, which.max)
  cur_predictions[cur_predictions == 1] <- "pass"
  cur_predictions[cur_predictions == 2] <- "run"
  all_preds[t] <- cur_predictions
}
data_test_NE$pred_MSDT <- all_preds

mean(data_test_NE$play_type != data_test_NE$pred_MSDT)

## Predictions with MS-DT for CLE
cur_match <- data_test_CLE
all_preds <- rep(NA, nrow(cur_match))
# t = 1
cur_predictions_01 <- msdt_forecast(xf = 0:1, h = 1, dataset = cur_match[1, ], mod = mod_CLE, first_obs = TRUE)
cur_predictions <- apply(cur_predictions_01, 1, which.max)
cur_predictions[cur_predictions == 1] <- "pass"
cur_predictions[cur_predictions == 2] <- "run"
all_preds[1] <- cur_predictions
for(t in 2:nrow(cur_match)){
  cur_predictions_01 <- msdt_forecast(xf = 0:1, h = 1, dataset = cur_match[1:t, ], mod = mod_CLE, first_obs = FALSE)
  cur_predictions <- apply(cur_predictions_01, 1, which.max)
  cur_predictions[cur_predictions == 1] <- "pass"
  cur_predictions[cur_predictions == 2] <- "run"
  all_preds[t] <- cur_predictions
}
data_test_CLE$pred_MSDT <- all_preds

mean(data_test_CLE$play_type != data_test_CLE$pred_MSDT)
