#customize theme
theme_new <- theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line=element_line(colour='black'),
        strip.background = element_rect(color = 'white', fill = 'white'),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "top")
#color configuration
mycolors =  c('#969696','#252525') 
mycolors3 = c('#d9d9d9','#969696','#252525') 
mycolors4 = c('#252525','#636363','#cccccc','#969696') 
mycolors5 = c('#d9d9d9', '#bdbdbd', '#969696', '#636363', '#252525')


colorSet2 <- scale_color_manual(values=c("#999999", "#E69F00")) 
colorSet3 <- scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) 
colorSet4 <- scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9", "#1a9641")) 
colorSet <- scale_color_manual(values=c("#d7191c", "#fdae61", "#a6d96a", "#1a9641")) 
colorSet5 <- scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9", "#1a9641", '#e6550d')) 
colorSet6 <- scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9", "#1a9641", '#7fcdbb', '#e6550d')) 

myshapevalues = c(5, 16, 17)
myshapevalues7 = c(10, 5, 3, 8, 11,9,15 )
