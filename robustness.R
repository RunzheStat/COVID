version.name = "final" # "0529_120_t0_12"


rep = 100
ws = exp(seq(-2, 6, 2)) / 10
T_begin = 12
decision.t = 7

M = 100
T_end = 120
n_city = 6
lag = 9
t0 = 61


################################################################################################################
################################################################################################################
organize_data <- function(main_path, sub_path){
  l_ws = 5
  ws = exp(seq(-2, 6, 2)) / 10
  
  save.path = paste(main_path, "_", sub_path, ".txt", sep = "")
  r.ours = readRDS(save.path)
  
  res = list()
  for (i in 1:125){
    seed = i %/% l_ws 
    w.index = i %% l_ws + 1
    w = ws[w.index]
    if(i <= 5){
      res[[w.index]] = r.ours[[i]]$RC
    }else{
      res[[w.index]] = res[[w.index]] + r.ours[[i]]$RC
    }
    
  }
  
  
  for (w in 1:5){
    res[[w]] = res[[w]] / 25
    res[[w]][1,] = log10(res[[w]][1,])
  }
  return(res)
}

# res = list()
# for (lag in 7:11){
#   res[[lag]] = organize_data("/Users/wrunzhe/Google Drive/Confident/AAAI_COVID/robust_lag/0828_tuning_result", sub_path = lag)
# }

para = 1
res = list()

for (effect2 in seq(5, 20, 5)) {
  for (effect3_effect2 in seq(5, 20, 5)) {
    effect3 = effect3_effect2 + effect2
    res[[para]] = organize_data("./robust_effect/0828_tuning_effect_result", sub_path = paste(effect3, effect2, sep = "_"))
    para = para + 1
  }
}
para.all = para - 1



################################################################################################################
################################################################################################################
## Paths
main_path = "./"
miceadds::source.all(paste(main_path, "codes", sep = ""))
out_path = paste(main_path, "output/", sep = "")
simu.res.path = paste(main_path, "output/", version.name, ".RDS", sep = "")
simu.plot.res.path = paste(main_path, "output/", version.name, "_plotres.RDS", sep = "")
plot.filename = paste(out_path, "output/", version.name, ".pdf", sep = "")
simu.monitor.path = paste(main_path, "output/", version.name, "__monitor.txt", sep = "")
simu.competing.monitor.path = paste(main_path,
                                    "output/",
                                    version.name,
                                    "_competing__monitor.txt",
                                    sep = "")
## Code
miceadds::source.all(paste(main_path, "codes", sep = ""))
## Reading Data
data = data_cn = readRDS(paste(main_path, "data.rds", sep = ""))
name = names(data_cn)


################################################################################################################
## Plotting: get results for each method  ------------------------------------------------------------

simu.plot.res.path = "./output/0606_plotres.RDS"

# r_oracle = getParetoCities("oracle", oneSettingMultiSeeds, ws)
plot.res = readRDS(simu.plot.res.path)

r_thre_XN_R = plot.res$r_thre_XN
r_thre_T3_R = plot.res$r_thre_T3
r_thre_T2_R = plot.res$r_thre_T2
r_ours_R = plot.res$r_ours

inlog = T
origin_zero = F
point_size = 2

#######################################################

################################################################################################################
################################################################################################################

# lag, w, (2, cities)

plotParetoFrontierforSensitivityAnalysis <-function(ours, paras, point_size = 1, line_size = 1 
                                     , add_sd
                                     , sd_ours
                                     , sd_thre_XN
                                     , sd_thre_T3
                                     , sd_thre_T2
                                     , r_thre_XN
                                     , r_thre_T3
                                     , r_thre_T2
                                     , RC_fixed, behav, city_name, origin_zero = F, inlog = F, echo = F
                                     , adjust_10 = T # at first, there is a 10 factor in our definition, and afterwards we removed it. to avoid re-running everything, we adjust it here.
                                     ){
  figure_size = point_size
  citytype = c(rep(city_name, length(C)))
  theme_set(theme_bw())
  cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
            "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  # cbp1 <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
  #           "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  colors <- c( "Pareto" = cbp1[4], "oracle" = cbp1[8], "threshold" = cbp1[2], "RC_fixed"=cbp1[3],
               "mitigation" = cbp1[6], "suppresion" = cbp1[1], "behaviour" = cbp1[7])
  shapes <- c("Pareto" = 16, "oracle" = 1
              , "threshold" = 5,"mitigation"= 2, "suppresion" = 0
              , "behaviour" = 13)
  
  ##########################################################################
  # ours
  ##########################################################################
  # scale_color_manual(name="Action",labels=c(1,2,3), values=colors, drop=FALSE) +
  # scale_shape_manual(name="Action",labels=c(1,2,3), values=c(15, 16, 17), drop=FALSE) +
  # , color=Action, shape=Action
  gg = ggplot()
  R_ours = c()
  C_ours = c()
  for (i in 1:length(paras)) {
    
    R = ours[[paras[i]]][1,]
    
    # if (inlog)R = log10(R) # already
    C = ours[[paras[i]]][2,]
    
    R_ours = c(R_ours, R)
    C_ours = c(C_ours, C)
    citytype = c(rep(city_name, length(C)))
    dat = data.frame(R = R, C = C, citytype)
    
    dat$policy = "Pareto"
    gg = gg +
      geom_point(data = dat, aes(x=R, y=C, colour=policy, shape=policy), size = figure_size * 1.3) + # 1.3
      geom_line(data = dat, aes(x=R, y=C, colour=policy, shape=policy), show.legend = FALSE, alpha = 1, size = line_size) +
      xlab("Count of cumulative infected (in log scale)") + ylab("Loss of GDP (billion Chinese Yuan)")+
      # ggtitle(city_name) +
      facet_wrap(~citytype, scale="free_y", dir="v") +
      scale_colour_manual(name="Policy",values = colors, drop=F)+
      scale_shape_manual(name="Policy",values = shapes, drop=F)+
      scale_fill_manual(name="Policy",values = colors, drop=F)+
      labs(colour  = "Policy", shape = "Policy", fill="Policy", drop=F) +
      # theme(plot.title = element_text(hjust = 0.5)) +
      theme(plot.caption = element_text(color = "#D6604D", face = "italic")   ) 
    if(add_sd){gg = gg + geom_rect(data=dat, aes(xmin=R-sd_ours[1,]*1.96, xmax=R+sd_ours[1,]*1.96, ymin=C-sd_ours[2,]*1.96, 
                                                 ymax=C+sd_ours[2,]*1.96, fill=policy, colour=policy),  alpha=0.1) }
    if(echo){print(gg)}
  }
  
  
  ##########################################################################
  # competing
  ##########################################################################
  
  if (!missing(r_thre_T2)) {
    R_thre_T2 = r_thre_T2[1,]
    if (inlog)R_thre_T2 = log10(R_thre_T2)
    C_thre_T2 = r_thre_T2[2,]
    if(adjust_10)C_thre_T2 = C_thre_T2 * 10
    dat_thre_T2 = data.frame(R = R_thre_T2, C = C_thre_T2)
    dat_thre_T2$policy = "mitigation"
    if(add_sd){
      dat_thre_T2$sdk_x = sd_thre_T2[1,]*1.96;
      dat_thre_T2$sdk_y = sd_thre_T2[2,]*1.96;
    }
    
    gg =   gg + geom_point(data = dat_thre_T2, aes(x = R, y = C, color= policy, shape=policy), size=figure_size) + 
      geom_line(data = dat_thre_T2, aes(x = R, y = C, color= policy, shape=policy)
                , alpha = 1, show.legend = FALSE
                , size = line_size) 
    if(add_sd){
      gg = gg + geom_rect(data = dat_thre_T2, aes(xmin=R-sdk_x, xmax=R+sdk_x, ymin=C-sdk_y, 
                                                  ymax=C+sdk_y, fill=policy, colour=policy),  alpha=0.1)
    }
    if(echo){print(gg)}
  }
  
  if (!missing(r_thre_T3)) {
    R_thre_T3 = r_thre_T3[1,]
    if (inlog)R_thre_T3 = log10(R_thre_T3)
    C_thre_T3 = r_thre_T3[2,]
    if(adjust_10)C_thre_T3 = C_thre_T3 * 10
    dat_thre_T3 = data.frame(R = R_thre_T3, C = C_thre_T3)
    dat_thre_T3$policy = "suppresion"
    gg = gg + 
      geom_point(data = dat_thre_T3, aes(x = R, y = C, color= policy, shape=policy), size=figure_size) + 
      geom_line(data = dat_thre_T3, aes(x = R, y = C, color= policy, shape=policy)
                , alpha = 1, show.legend = FALSE
                , size = line_size)
    if(add_sd){
      gg = gg + geom_rect(data = dat_thre_T3, aes(xmin=R-sd_thre_T3[1,]*1.96, xmax=R+sd_thre_T3[1,]*1.96, ymin=C-sd_thre_T3[2,]*1.96, 
                                                  ymax=C+sd_thre_T3[2,]*1.96, fill=policy, colour=policy),  alpha=0.1)
    }
    if(echo){print(gg)}
  }
  
  if (!missing(behav)) {
    if(inlog) behav[1] = -log10(-behav[1])
    dat_behav = data.frame(R=-behav[1], C=behav[2], policy="behaviour")
    
    if(add_sd){
      gg = gg + geom_point(dat_behav, mapping = aes( x = R, y = C , color= policy, shape=policy, fill=policy), size= figure_size + 1) # , shape = 5
    } else{
      gg = gg + geom_point(dat_behav, mapping = aes( x = R, y = C , color= policy, shape=policy), size= figure_size + 1) # , shape = 5
    }
    if(echo){print(gg)}
  }  
  # dat_behav = data.frame(x =  -behav[1], y = behav[2])
  
  if (!missing(r_thre_XN)) {
    R_thre_XN = r_thre_XN[1,]
    if (inlog)R_thre_XN = log10(R_thre_XN)
    C_thre_XN = r_thre_XN[2,]
    if(adjust_10)C_thre_XN = C_thre_XN * 10
    dat_thre_XN= data.frame(R = R_thre_XN, C = C_thre_XN)
    dat_thre_XN$policy = "threshold"
    gg = gg + 
      geom_point(data = dat_thre_XN, aes(x = R, y = C
                                         , color = policy, shape = policy) , size=figure_size) + 
      geom_line(data = dat_thre_XN, aes(x = R, y = C, color = policy, shape = policy)
                ,  alpha = 1, show.legend = FALSE
                , size = line_size)
    if (add_sd){
      gg = gg + geom_rect(data = dat_thre_XN, aes(xmin=R-sd_thre_XN[1,]*1.96, xmax=R+sd_thre_XN[1,]*1.96, ymin=C-sd_thre_XN[2,]*1.96, 
                                                  ymax=C+sd_thre_XN[2,]*1.96, fill=policy, colour=policy),  alpha=0.1)
    }
    
    if(echo){print(gg)}
  }
  
  
  ##########################################################################
  # limit
  ##########################################################################
  if(echo){print(gg)}
  xs = c(R_ours, -behav[1], R_thre_XN, R_thre_T3, R_thre_T2)
  ys = c(C_ours, behav[2], C_thre_XN, C_thre_T3, C_thre_T2)
  if (add_sd) {
    sds = cbind(sd_ours, c(0,0), sd_thre_XN, sd_thre_T3, sd_thre_T2) *2
    sdsx = sds[1,];  sdsy = sds[2,]
  }else{
    sdsx = 0
    sdsy = 0
  }
  
  const = 1
  xsmin = const * ( min(xs) - max(sdsx) )
  xsmax = const * ( max(xs) + max(sdsx) )
  ysmin = const * (min(ys) - max(sdsy))
  ysmax = const * (max(ys) + max(sdsy))
  if(origin_zero){
    gg = gg + xlim(0, xsmax) +
      ylim(0, ysmax)
    
  }else{
    gg= gg + xlim(xsmin, xsmax) +  ylim(ysmin, ysmax)
  }
  gg = gg +  theme(axis.title=element_blank()) + expand_limits(x=c(1,4))
  return(gg)
}



################################################################################################################
################################################################################################################
# lag, w, (2, cities)

# range_tuning_para = seq(7, 11, 1)
range_tuning_para = c(1:para.all)
range_w = c(1:length(ws))

# FIGURE 3 in paper
ggs = list()

for (l in 1:n_city) {
  ours = list()
  for (para in range_tuning_para) {
    ours[[para]] = matrix(NA, 2, length(range_w))
    for (w in range_w){
      ours[[para]][, w] = res[[para]][[w]][, l]
    }
  }
  ### we need to provide a list (len-range_tuning_para) of RC = [2, range_w]
  
  ggs[[l]] = plotParetoFrontierforSensitivityAnalysis(
    ours = ours,
    paras = range_tuning_para, 
    point_size = point_size,
    line_size = .8,
    behav = evaluateBehavPolicySimu(data, l, T_begin, 120, lag),
    add_sd = F,
    #sd_ours = r_ours_ud$sd[[l]] / 10, # / sqrt{100}
    #sd_thre_XN = r_thre_XN_ud$sd[[l]] / 10,
    #sd_thre_T3 = r_thre_T3_ud$sd[[l]] / 10,
    # sd_thre_T2 = r_thre_T2_ud$sd[[l]] / 10,
    r_thre_XN = r_thre_XN_R[[l]],
    r_thre_T3 = r_thre_T3_R[[l]],
    r_thre_T2 = r_thre_T2_R[[l]],
    city_name = name[l],
    origin_zero = origin_zero,
    inlog = inlog
  )
  print(ggs[[l]])
}

############################################################
combinePlots <-
  function(ggs,
           plot.filename,
           width = 6,
           height = 8) {
    library(ggpubr)
    
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

combinePlots(ggs, 
             plot.filename = paste(out_path, "robust_effect.png", sep = "")
             , 6 / 1.2, 8 / 1.2)
