rm(list=ls())
################################################################################################################
################################################################################################################
version.name = "H1N1"

u_prior = c(1, .2, .1)
rep = 100
# ws = exp(seq(-4, 9, 1))
ws = exp(seq(-2, 6, 1))
T_begin = 12
decision.t = 7



# parameters for COMPETINGS

vec.thre.2 = seq(5, 50, 5)
thress = list()
for (i in c(1:length(vec.thre.2))) { # why problem in ggplot?
  thress[[i]] = c(0, vec.thre.2[i])
}

# thre.Ts.2 = seq(4, 14, 2)
# thre.Ts.3 = seq(4, 10)
thre.Ts.2 = seq(4, 12, 2)
thre.Ts.3 = seq(4, 10, 2)

###########

stochastic = F; generate.belife.every.step = T
use.model.for.ours = F; catious_c =  1
competing.methods = c("thre_XN", "thre_T2", "thre_T3") # "over", "mediation", 
M = 100
T_end = 120
city_index = c(1,3,4,5,6,7); n_city = length(city_index)
lag = 9
t0 = 61

################################################################################################################
################################################################################################################

# PC
main_path = "./"

main_path = "/Users/mac/Desktop/CHIL_code"
## Paths
out_path = paste(main_path,"out/",sep = "")
simu.res.path = paste(main_path, "out/simu_res/", version.name, ".RDS", sep = "")
simu.plot.res.path = paste(main_path, "out/simu_res/", version.name, "_plotres.RDS", sep = "")
plot.filename = paste(out_path, "simu_fig/", version.name, ".pdf",sep = "")
simu.monitor.path = paste(main_path, "out/simu_monitor/", version.name, "__monitor.txt", sep = "")
simu.competing.monitor.path = paste(main_path, "out/simu_monitor/", version.name, "_competing__monitor.txt", sep = "")

## Code
miceadds::source.all(paste(main_path,"code in package/",sep = ""))
miceadds::source.all(paste(main_path,"simu",sep = ""))
## Reading Data
data = data_cn = readRDS(paste(main_path,"/data.rds",sep = ""))
name = names(data_cn)
#### model 
model = NA

################################################################################################################


n_core = parallel::detectCores();n_core

theta = getSIRMAP(data=data_cn, MAP.only = F, t0 = t0, region_index = city_index,# brief = brief, 
                  lag = lag, alpha = 0.1, u=u_prior, echo = F, brief = F)[1:7, 1]; theta


##########################################################
########## Do something to theta here! ###################
theta[2:7] = theta[2:7] / 2.3179 * 1.4

##########################################################

g_C <-function(A, GDP = NA, 
               vec_c = c(0, 0.368, 0.484), sigma_c = c(0, 0.239, 0.181),
               deterministic = F, seed = 1,
               ratio = 0.1){
  if (deterministic)return(vec_c[A] * GDP * ratio)
  set.seed(seed)
  r = rnorm(1, vec_c[A], sigma_c[A])
  if(is.na(r)){
    cat(A, r)
  }
  return(r * GDP * ratio)
}

################################################################################################################
## COMPETING METHODS
################################################################################################################


oneSeedCompeting <- function(seed){
  r = simuCompeting(seed = seed, 
                    model = model, data = data_cn
                    , T_begin = lag + 1, T_end = T_end, lag = 9
                    , n_city = n_city, city_index = city_index
                    , ws = ws
                    , u_prior = u_prior, theta = theta
                    , stochastic = stochastic
                    , generate.belife.every.step = generate.belife.every.step
                    , decision.t = decision.t, catious_c = catious_c, M = M
                    , thress = thress
                    , thre.Ts.2 = thre.Ts.2
                    , thre.Ts.3 = thre.Ts.3
                    , methods = competing.methods
                    , MC = F)
  return(r)
}

# quick competing methods
cl <- parallel::makeCluster(n_core, outfile = simu.competing.monitor.path)
doParallel::registerDoParallel(cl)
oneSettingMultiSeeds = foreach(i = 1:rep, .packages = c("matrixStats")) %dopar% {
  oneSeedCompeting(i)
}   ## LLL: longtime to check

saveRDS(oneSettingMultiSeeds, simu.res.path)
oneSettingMultiSeeds = readRDS(simu.res.path)
stopCluster(cl)

r_thre_XN = getParetoCities("thre_XN", oneSettingMultiSeeds, thress)
r_thre_T2 = getParetoCities("thre_T2", oneSettingMultiSeeds, thre.Ts.2)
r_thre_T3 = getParetoCities("thre_T3", oneSettingMultiSeeds, thre.Ts.3)

print(r_thre_XN[[1]])
print(r_thre_T2[[1]])
print(r_thre_T3[[1]])



################################################################################################################
## OURS
################################################################################################################

onceSeedWs <-function(i = 1, method =  "ours"){
  # iter over weights first
  seed = i %/% length(ws)
  w = ws[i %% length(ws) + 1]
  r = simuOnce(seed = seed, data = data, 
               w = w, model = model, method = method, 
               u_prior =  u_prior
               , use.model.for.ours = use.model.for.ours
               , n_city =  n_city , city_index = city_index
               , M = M
               , stochastic = stochastic , generate.belife.every.step = generate.belife.every.step
               , lag = lag
               , T_begin =  lag + 1, T_end =  T_end
               , theta =  theta
               , brief = T
               , catious_c = catious_c
               , decision.t = decision.t
  )
  # cat("\n", "<--------------------------- seed = ", seed, "----|| log(omega) = ", log(w), "||---- i = " , i , "---------------------------->\n")
  cat("\n", "<--------------------------- (i, seed, log(omega) = ", i, seed, log(w) , "---------------------------->\n")
  cat("R:", round(r$RC[1, ]),"\n")
  cat("C:", round(r$RC[2, ]),"\n")
  cat("As:", round(colMeans(r$As)), "\t ----------------------------> \n")
  return(r)
}

reps.ours = (rep * length(ws))
cat("simu.monitor.path:", simu.monitor.path, "\n")
print(Sys.time()); start_time = Sys.time() # output: a (len-rep) list of len-n_method lists. If Pareto, then a len_n_weight list of RC = matrix(NA, 2, n_city)
cl <- parallel::makeCluster(n_core, outfile = simu.monitor.path)
doParallel::registerDoParallel(cl)
r.ours = foreach(i = 1:reps.ours, .packages = c("matrixStats")) %dopar% {
  onceSeedWs(i, method = "ours")
}
time.cost = as.character(Sys.time() - start_time)
print(Sys.time() - start_time)
stopCluster(cl)


oneSettingMultiSeeds_backup = oneSettingMultiSeeds
oneSettingMultiSeeds = putOurs_into_oneSettingMultiSeeds(oneSettingMultiSeeds,
                                                         r.ours, method =  "ours")
saveRDS(oneSettingMultiSeeds, simu.res.path)
r_ours = getParetoCities("ours", oneSettingMultiSeeds, ws)

plot.res = list(r_thre_XN = r_thre_XN, r_thre_T3 = r_thre_T3, r_thre_T2 = r_thre_T2, r_ours = r_ours)
saveRDS(plot.res, simu.plot.res.path)


################################################################################################################
## Plotting: get results for each method  ------------------------------------------------------------

# r_oracle = getParetoCities("oracle", oneSettingMultiSeeds, ws)
simu.plot.res.path = "/Users/mac/Desktop/CHIL_code/H1N1/H1N1_plotres.RDS"
simu.res.path = "/Users/mac/Desktop/CHIL_code/H1N1/H1N1.RDS"

plot.res = readRDS(simu.plot.res.path)

r_thre_XN_R = plot.res$r_thre_XN
r_thre_T3_R = plot.res$r_thre_T3
r_thre_T2_R = plot.res$r_thre_T2
r_ours_R = plot.res$r_ours

inlog = T
origin_zero = F
point_size = 2

# get sd
oneSettingMultiSeeds = readRDS(simu.res.path)
r_thre_XN_ud = getParetoCities("thre_XN", oneSettingMultiSeeds, thress, sd=TRUE, inlog = inlog)
r_thre_T2_ud = getParetoCities("thre_T2", oneSettingMultiSeeds, thre.Ts.2, sd=TRUE, inlog = inlog)
r_thre_T3_ud = getParetoCities("thre_T3", oneSettingMultiSeeds, thre.Ts.3, sd=TRUE, inlog = inlog)
r_ours_ud = getParetoCities("ours", oneSettingMultiSeeds, ws, sd=TRUE, inlog = inlog)

select_w = c(1:length(ws))

################################################################################################################
width=10;height=6



# FIGURE 3 in paper
ggs = list()
l=1
for (l in 1:n_city) { 
  r_ours_R[[l]][, select_w][2, ] = r_ours_R[[l]][, select_w][2, ] * 10
  r_thre_XN_R[[l]][2, ] = r_thre_XN_R[[l]][2, ] * 10
  r_thre_T3_R[[l]][2, ] = r_thre_T3_R[[l]][2, ] * 10
  r_thre_T2_R[[l]][2, ] = r_thre_T2_R[[l]][2, ] * 10
  
  
  ggs[[l]] = plotParetoFrontierforSimu(ours = r_ours_R[[l]][, select_w] # 9
                                       , point_size = point_size
                                       , line_size = .8
                                       # , behav = NA #evaluateBehavPolicySimu(data, l, T_begin, 120, lag)
                                       # , oracle  = r_oracle[[l]][, 1:select_w]
                                       , add_sd = FALSE
                                       , sd_ours = r_ours_ud$sd[[l]]/10 # 9
                                       , sd_thre_XN = r_thre_XN_ud$sd[[l]]/10 # 10
                                       , sd_thre_T3 = r_thre_T3_ud$sd[[l]]/10 # 4
                                       , sd_thre_T2 = r_thre_T2_ud$sd[[l]]/10 # 5
                                       , r_thre_XN = r_thre_XN_R[[l]] # 10
                                       , r_thre_T3 = r_thre_T3_R[[l]] # 4
                                       , r_thre_T2 = r_thre_T2_R[[l]] # 5
                                       , city_name = name[l]
                                       , origin_zero = origin_zero, inlog = inlog)
  print(ggs[[l]])
}

# gg = ggs[[1]]
# gg + guides(fill = guide_legend(override.aes = list(linetype = 0)) )
# ggsave(filename = "./out/simu_fig/try.png", plot = gg, units = c("in"), dpi = 500)



combinePlots <-
  function(ggs,
           plot.filename,
           width = 6,
           height = 8) {
    g = ggarrange(
      ggs[[1]],
      ggs[[2]],
      ggs[[3]],
      ggs[[4]],
      ggs[[5]],
      ggs[[6]],
      ncol = 2,
      nrow = 3,
      legend = "top",
      label.x = 0,
      label.y = 0,
      font.label = list(size = 20, face = "bold"),
      common.legend = T,
      align = "v"
    )
    yvar =  c("Loss of GDP (billion Chinese Yuan)")
    xvar = c(expression(
      paste("Count of cumulative infected (in ", log[10], " scale)", sep = "")
    ))
    g1 = annotate_figure(
      g,
      left = text_grob(
        yvar,
        just = "centre",
        size = 12,
        rot = 90,
        hjust = .5,
        vjust = .5
      ),
      bottom = text_grob(
        xvar,
        just = "centre",
        size = 12,
        rot = 0,
        hjust = .48,
        vjust = 0
      )
    )
    print(g1)
    ggsave(
      filename = plot.filename,
      plot = g1,
      width = width,
      height = height,
      units = c("in"),
      dpi = 500
    )
    
  }


figure_3_filename = paste(main_path,"/", "Pareto_H1N1.png", sep = "")
combinePlots(ggs, figure_3_filename, 6 / 1.2, 8 / 1.2)


# FIGURE 1 in supplementary


ggsd = list()
l=1
for (l in 1:n_city) { 
  ggsd[[l]] = plotParetoFrontierforSimu(ours = r_ours_R[[l]][, select_w] # 9
                                       , point_size = point_size
                                       , line_size = .8
                                       , behav = evaluateBehavPolicySimu(data, city_index[l], T_begin, 120, lag)
                                       # , oracle  = r_oracle[[l]][, 1:select_w]
                                       , add_sd = TRUE
                                       , sd_ours = r_ours_ud$sd[[l]]/10 # 9
                                       , sd_thre_XN = r_thre_XN_ud$sd[[l]]/10 # 10
                                       , sd_thre_T3 = r_thre_T3_ud$sd[[l]]/10 # 4
                                       , sd_thre_T2 = r_thre_T2_ud$sd[[l]]/10 # 5
                                       , r_thre_XN = r_thre_XN_R[[l]] # 10
                                       , r_thre_T3 = r_thre_T3_R[[l]] # 4
                                       , r_thre_T2 = r_thre_T2_R[[l]] # 5
                                       , city_name = name[city_index[l]]
                                       , origin_zero = origin_zero, inlog = inlog)
  print(ggsd[[l]])
}


figure_1_filename = paste(main_path,"out/simu_fig/", "figure_1_in_supplementary.pdf", sep = "")

combinePlots(ggsd, figure_1_filename, 10/1.2, 5.5/1.2)


# ----
#####  T3

# thre.Ts.3 = seq(2, 10)
# select.entrys = seq(4, 10, 2) - 1
# # select.entrys = seq(3, 10, 2) - 1
# 
# for (i in c(1:6)) {
#   r_thre_T3_R[[i]] = r_thre_T3_R[[i]][1:2, select.entrys]
#   r_thre_T3_ud$sd[[i]] = r_thre_T3_ud$sd[[i]][1:2, select.entrys]
# }
# thre.Ts.3 = thre.Ts.3[select.entrys]
# 
# 
# ##### T2
# thre.Ts.2 = seq(2, 14)
# select.entrys = seq(4, 12, 2) - 1
# # select.entrys = seq(3, 12, 3) - 1
# 
# for (i in c(1:6)) {
#   r_thre_T2_R[[i]] = r_thre_T2_R[[i]][1:2, select.entrys]
#   r_thre_T2_ud$sd[[i]] = r_thre_T2_ud$sd[[i]][1:2, select.entrys]
# }
# thre.Ts.2 = thre.Ts.2[ select.entrys]


##### thre
# vec.thre.2 = seq(5, 50, 5)
# thress = list()
# for (i in c(1:length(vec.thre.2))) { # why problem in ggplot?
#   thress[[i]] = c(0, vec.thre.2[i])
# }
# 
# select.entrys = seq(1,10,2)
# for (i in c(1:6)) {
#   r_thre_XN[[i]] = r_thre_XN[[i]][1:2, select.entrys]
# }
# thress = thress[select.entrys]

# ws = exp(seq(-4, 9, 1))
# select_w = c(4:10); ws = ws[select_w]; cat(log(ws))
