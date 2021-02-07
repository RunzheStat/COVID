library("ggpubr")
#' simuCompeting
#' for one seed, generate results for all methods (including all ws)
#' 
#' @return a list with each entry for one method. If Pareto, then a list of RC = matrix(NA, 2, n_city), each for a weight. 

prepareDataForSIR_est_with_SEIR <- function(trajectories) {
  n_city = length(trajectories)
  t0 = length(trajectories[[1]]$action + 1)
  data_SIR = list()
  for (i in 1:n_city) {
    data_i = list()
    data_i[["action"]] = c(trajectories[[i]]$action, NA)
    data_i[["infected"]] = trajectories[[i]]$state[4, ]
    data_i[["date"]] = trajectories[[i]][["date"]]
    data_i[["gdp"]] = trajectories[[i]]$gdp
    data_i[["population"]] = trajectories[[i]]$popu
    data_SIR[[i]] = data_i
  }
  return(data_SIR)
}

f_GSEIR <- function(X_lt, A, theta = c(0.15, 1 / 7, 0.78), seed = 1){
  # theta = [gamma, sigma, beta1, beta2, beta3]
  set.seed(seed)
  # get current state variables and parameters
  J = 3
  XS = X_lt[1]
  XE = X_lt[2]
  XI = X_lt[3]
  XR = X_lt[4]
  popu = sum(X_lt)
  
  gamma = theta[1]
  sigma = theta[1]
  beta = theta[2 + A]
  
  eS = rpois(1, max(beta * XS * XI / popu, 0))
  eE = rbinom(1, max(as.integer(XE), 0), sigma)
  e_R = rbinom(1, max(as.integer(XI), 0), gamma)
  
  XS_new = XS -  eS
  XE_new = XE +  eS - eE
  XR_new = XR + e_R
  XI_new = popu - XS_new - XR_new - XE_new
  
  
  r = c(XS = as.integer(XS_new)
        , XE = as.integer(XE_new)
        , XI = as.integer(XI_new)
        , XR = as.integer(XR_new))
  names(r) = c("XS", "XE", "XI","XR")
  return(r)
}


simuCompeting <- function(seed = 1, model = NA, data = NA, ws = exp(seq(-1, 6.7, 1)),
                          u_prior = c(1, .2, .1)
                          , city_index = NA
                          , n_city = 6, 
                          T_begin = 10, T_end = 120, 
                          theta = NA, brief = T
                          , catious_c = 1
                          , decision.t = 7
                          , MC = F
                          , stochastic = T
                          , generate.belife.every.step = F
                          , M = 100
                          , lag = 9
                          , methods = NA
                          , thress = NA
                          # , thre.Ts = NA
                          , thre.Ts.2 = NA
                          , thre.Ts.3 = NA
){

  rs = list()

  if ("thre_XN" %in% methods) {
    once <-function(thres){
      return(simuOnce(seed = seed, data = data, 
                      w = 1, model = model, method = "thre_XN"
                      , thres = thres
                      , lag = lag
                      , u_prior = u_prior
                      , stochastic = stochastic
                      , generate.belife.every.step = generate.belife.every.step
                      , n_city =  n_city, city_index = city_index
                      , M = M, T_begin =  T_begin, T_end =  T_end, 
                      theta =  theta)
      )
    }
    if (MC) {
      r_threshold = mclapply(thress, once)
    }else{
      r_threshold = lapply(thress, once)  
    }
    rs[["thre_XN"]] = r_threshold
    cat("<----------------------", "thre_XN DONE!", "----------------------> \n")
  }
  
  if ("thre_T3" %in% methods) {
    once <-function(thre.Ts.3){
      return(simuOnce(seed = seed, data = data, 
                      w = 1, model = model, method = "thre_T3"
                      , lag = lag
                      , u_prior = u_prior
                      , stochastic = stochastic
                      , generate.belife.every.step = generate.belife.every.step
                      , n_city =  n_city, city_index = city_index
                      , M = M, T_begin =  T_begin, T_end =  T_end
                      , thre.Ts.3 = thre.Ts.3
                      , theta =  theta)
      )
    }
    if (MC) {
      r_threshold = mclapply(thre.Ts.3, once)
    }else{
      r_threshold = lapply(thre.Ts.3, once)  
    }
    rs[["thre_T3"]] = r_threshold
    cat("<----------------------", "thre_T3 DONE!", "----------------------> \n")
  }
  
  if ("thre_T2" %in% methods) {
    once <-function(thre.Ts.2){
      return(simuOnce(seed = seed, data = data, 
                      w = 1, model = model, method = "thre_T2"
                      , lag = lag
                      , u_prior = u_prior 
                      , stochastic = stochastic
                      , generate.belife.every.step = generate.belife.every.step
                      , n_city =  n_city, city_index = city_index
                      , M = M, T_begin =  T_begin, T_end =  T_end
                      , thre.Ts.2 = thre.Ts.2
                      , theta =  theta)
      )
    }
    if (MC) {
      r_threshold = mclapply(thre.Ts.2, once)
    }else{
      r_threshold = lapply(thre.Ts.2, once)  
    }
    rs[["thre_T2"]] = r_threshold
    cat("<----------------------", "thre_T2 DONE!", "----------------------> \n")
  }
  return(rs)
}



#' simuOnce
#' for each policy (e.g., a fixed w), generate one replication for all cities
#' 
#' @param policy "ours", "oracle", "over", "mediation"
#' @param
#' @return RC = matrix(NA, 2, n_city)


simuOnce <-function(data, 
                    model, 
                    seed = 1
                    , method = "ours", w = 1
                    , u_prior =c(1, .2, .1) 
                    , city_index = NA
                    , n_city = 9, T_begin = 9, T_end = 61 
                    , theta = c(0.12, 0.29, 0.11, 0.04) , brief = T
                    , MC = F
                    , lag = 9
                    , M = 100
                    , use.model.for.ours = T
                    , stochastic = T
                    , generate.belife.every.step = F
                    , catious_c = 0.8
                    , decision.t = 7
                    , thres = c(10, 30)
                    , thre.Ts.2 = NA
                    , thre.Ts.3 = NA){
  
  
  ## Initialization ----
  set.seed(seed) 
  trajectories = list()
  
  s_seed = seed * 123
  
  ## Initialize trajectories
  for (i in 1:n_city) {
    city_data = data[[city_index[i]]]
    city = list()
    # generate ture states [think!]
    init_state = matrix(NA, 4, (lag + 1))
    init_state[4, ] = city_data$infected[(T_begin - lag):T_begin]
    init_state[3, ] = city_data$infected[(T_begin):(T_begin + lag)] - init_state[4, ]
    init_state[2, ] = init_state[3, ]
    init_state[1, ] = city_data$population[1] - init_state[2, ] - init_state[3, ] - init_state[4, ]
    city[["state"]] = init_state
    city[["action"]] = city_data$action[(T_begin - 1 - (lag - 1)):(T_begin - 1)]
    city[["popu"]] = city_data$population[1]
    city[["gdp"]] = (city_data$gdp /1e9 / 365)[1]
    city[["date"]] = city_data$date[(T_begin - lag):T_begin]
    city[["R"]] = 0
    city[["C"]] = 0
    city[["init.flag"]] = F
    trajectories[[i]] = city
  }
  
  ## all cities begin to move on together, following our workflow

  for (t0 in T_begin:(T_end-1)) {
    
    data.SIR = prepareDataForSIR_est_with_SEIR(trajectories)
    hat.theta.t0 = getSIRMAP(data = data.SIR, region_index = c(1:n_city), 
                             t0 = t0, lag = lag,  u = u_prior, 
                             for_simu = T)[1:8]
    hat.theta.t0.mean = c(hat.theta.t0[1] / (hat.theta.t0[1] + hat.theta.t0[2]), 
                          hat.theta.t0[3] / hat.theta.t0[6],
                          hat.theta.t0[4] / hat.theta.t0[7],
                          hat.theta.t0[5] / hat.theta.t0[8])
    
    if(t0 %% 20 == 0)cat(t0, "|| (seed, w) = ", seed, log(w), 
                         "|| R0 = ", round(hat.theta.t0.mean[2:4] / hat.theta.t0.mean[1], 2),"\n")
    
    s_seed = s_seed + 1
    
    if (use.model.for.ours) {
      search.method = "pre"
    }else{
      search.method = "grid"
    }
    ### 
    updateCity_i <- function(city_i, is.A.point = T){
      
      if (is.A.point) {
        if (method == "ours") {
          A = policySearch(theta = hat.theta.t0, t0 = t0, T_end = T_end,
                           XR_t0_lag = city_i$state[4, (t0 - lag):t0 ], # city_46, chengdu, 2:10
                           A_t0_1_lag = city_i$action[(t0 - lag):(t0-1)], 
                           popu = city_i$popu, 
                           GDP_daily = city_i$gdp,
                           w = w
                           , lag = lag
                           , M = M
                           , only_action = T
                           , stochastic = stochastic
                           , generate.belife.every.step = generate.belife.every.step
                           , J = 3, seed = s_seed, echo = F,
                           method = search.method, model = model
                           , decision.t = decision.t
                           , catious_c = catious_c
          )
        }else if (method == "DQN"){
          init_state = generateBeliefState(hat.theta.t0, pre_XRs = city_i$state[4, (t0 - lag):t0]
                                           , pre_As = city_i$action[(t0 - lag):(t0-1)]
                                             , popu = city_i$popu
                                           , stochastic = stochastic)
          init_state = round(init_state)
          Xl = city_i$gdp
          # 1. Agent disappears
          # 2. Multiple regions
          if (t0 == T_begin) {
            gsir = GSIR(init_state, hat.theta.t0, w, T_end, Xl, t0)
            agent = Agent(gsir, timesteps_per_episode = T_end, batch_size = 64)
            A = agent$train(100) + 1
          }else{
            gsir = GSIR(init_state, hat.theta.t0, w, T_end, Xl, t0)
            agent = Agent(gsir, timesteps_per_episode = T_end, batch_size = 64)
            # agent$update(init_state, t0, hat.theta.t0)
            A = agent$train(100) + 1
          }
        }else if (method == "oracle") {
          A = policySearch(theta = theta, t0 = t0, T_end = T_end,
                           XR_t0_lag = city_i$state[4, (t0 - lag):t0 ], # city_46, chengdu, 2:10
                           A_t0_1_lag = city_i$action[(t0 - lag):(t0-1)], 
                           S_t0 = city_i$state[c(1, 3, 4), t0], 
                           is.state.oracle = T, 
                           popu = city_i$popu, 
                           GDP_daily = city_i$gdp,
                           w = w
                           , only_action = T
                           , stochastic = stochastic
                           , generate.belife.every.step = generate.belife.every.step
                           , lag = lag
                           , M = M
                           , J = 3, seed = s_seed, echo = F,
                           method = search.method, model = model
                           , decision.t = decision.t
                           , catious_c = catious_c
          )
        }#else 
      }else{
        A = city_i$action[length(city_i$action)] # last action
      }
      
      if (method == "thre_XN") {
        
        # daily.new = city_i$state[3, t0] - city_i$state[3, t0 - 7]
        daily.new = city_i$state[4, t0] - city_i$state[4, t0 - 1]
        if (daily.new > thres[2]){ 
          A = 3
        }else if(daily.new > thres[1]){
          A = 2
        }else{
          A = 1
        }
      }else if (method == "thre_T2") {
        
        any.new.recent =  (city_i$state[4, t0] - city_i$state[4, max(t0 - thre.Ts.2,1)] > 0)
        
        
        if (any.new.recent) {
          A = 2
        }else{
          A = 1
        }
      }else if (method == "thre_T3") {
        T.first.case.city.i = min(which(city_i$state[4, 1:t0] > 0))
        
        no.any.new.recent.thre.T = (city_i$state[4, t0] - city_i$state[4, max(t0 - thre.Ts.3, 1)] == 0)
        
        no.any.new.recent.thre.2T = (city_i$state[4, t0] - city_i$state[4, max(t0 - 2 *  thre.Ts.3, 1)] == 0)
        
        if (t0 - T.first.case.city.i > thre.Ts.3) {
          if (t0 - T.first.case.city.i > 2  * thre.Ts.3) {
            if (no.any.new.recent.thre.2T) {
              A = 1
            }else if(no.any.new.recent.thre.T){
              A = 2
            }else{
              A = 3
            }
          }else{
            A = 2
          }
          
        }else{
          A = 1
        }
        
      }
      
      city_i$action = c(city_i$action, A)
      city_i$date = city_data$date[(T_begin  - lag): (t0 + 1)]
      ## observe next state and R/C
      # collect true current state and sample the next state  (and thus reward/cost)
      X_lt = city_i$state[, t0]
      X_lt_old = X_lt
      
      X_lt = f_GSEIR(X_lt = X_lt, A = A, theta = theta, seed = s_seed)
      city_i$state = cbind(city_i$state, X_lt)
      YRC_one_step = g(X_lt_old, A, X_lt, city_i$gdp, w, seed = s_seed)
      
      city_i$R = city_i$R + YRC_one_step[2]
      city_i$C = city_i$C + YRC_one_step[3]
      
      return(city_i)
      
    }
    
    if (t0 %% decision.t == (T_begin %% decision.t))is.A.point = T
    
    if (MC) {
      trajectories = mclapply(trajectories, updateCity_i, is.A.point = is.A.point)
    }else{
      trajectories = lapply(trajectories, updateCity_i, is.A.point = is.A.point)
    }
    
    is.A.point = F
    
  }
  
  ## output
  
  RC = matrix(NA, 2, n_city)
  As = list()
  for (i in 1:n_city) {
    RC[1, i] = -trajectories[[i]]$R
    RC[2, i] = trajectories[[i]]$C
    As[[i]] = trajectories[[i]]$action
  }
  
  As = do.call(rbind, As)
  
  if (brief) {
    # r = RC
    r = list(RC = RC, As = As)
  }else{
    for (j in 1:length(trajectories)) {
      trajectories[[j]] = rbind(trajectories[[j]]$state[3:4, 1:(T_end-1)], trajectories[[j]]$action)
      trajectories[[j]] =  as.data.frame(trajectories[[j]])
      names(trajectories[[j]]) = c(1:(length(trajectories[[1]]) - 1))
      rownames(trajectories[[j]]) = c("I", "R", "A")
    }
    r = list(RC = RC, trajectories = trajectories)
  }
  
  return(r)
  
}

