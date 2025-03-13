#percentage
#Percentage_16_ECMHomeoMac.xlsx
library(tidyverse)
library(ggsci)

Tissue <- c("Normal","Tumor","Carotid artery proximal adjacent (PA)","Carotid artery atherosclerotic core (AC)")
Count_16_ECMHomeoMac <- c("429", "9347", "0", "2302") %>% as.numeric()
Count_Others <- c("74553","269757","989","8491") %>% as.numeric()

df <- data.frame(Tissue, Count_16_ECMHomeoMac,Count_Others      )
df$Percentage <- df$Count_16_ECMHomeoMac/(df$Count_16_ECMHomeoMac+df$Count_Others) *100
df$Percentage <- round(df$Percentage,2)

df$Tissue <- factor( df$Tissue,levels =  c("Normal","Tumor","Carotid artery proximal adjacent (PA)","Carotid artery atherosclerotic core (AC)")  )
df

#                                    Tissue Count_16_ECMHomeoMac Count_Others Percentage
#1                                   Normal                  429        74553       0.57
#2                                    Tumor                 9347       269757       3.35
#3    Carotid artery proximal adjacent (PA)                    0          989       0.00
#4 Carotid artery atherosclerotic core (AC)                 2302         8491      21.33

p = ggplot(df, aes(x = Tissue, y = Percentage, fill = Tissue)) +
  geom_bar(stat = "identity", show.legend = TRUE) + # 显示图例
  scale_fill_manual(values = c( "#5494ce","#CEA876","#0d898a","#A47053"   )) +
  scale_y_continuous(limits = c(0, 23)) +
  theme_bw(base_size = 16) + # 设置基本字体大小为16
  labs(title = "Percentage of Cluster16 across different tissues",
       x = "Tissue Type", 
       y = "Percentage (%)") +
  theme(axis.text.x = element_blank(),  # 旋转X轴标签
        axis.text = element_text(size = 14),  # 增加轴标签字体大小
        axis.title = element_text(size = 16), # 增加轴标题字体大小
        plot.title = element_text(size = 18, face = "bold"))+ # 增加标题字体大小并加粗
        
  #theme(plot.margin = unit(c(0.1,0.1,0.1,1.5), "inches")) +
  geom_text(aes(label = paste0(Percentage, "%")), vjust = -0.5, size = 5)

ggsave(p,file = "plots/Percentage_BarPlot.jpeg",height = 8,width = 8)
