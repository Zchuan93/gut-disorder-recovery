setwd('C:\\Users\\zhang\\Desktop\\2维PCoA绘图')
# Load package
library(vegan)
library(ggplot2)
library(ggthemes)
library(tidyverse)
# Load data
otu <- read.table('otu.txt',row.names = 1,sep="\t",header = T)
group <- read.table('group.txt',header = T)
# creat data
group$bacteria <- runif(128,0,20)
#pcoa
# vegdist函数，计算距离；method参数，选择距离类型
distance <- vegdist(otu, method = 'bray')
# 对加权距离进行PCoA分析
pcoa <- cmdscale(distance, k = (nrow(otu) - 1), eig = TRUE)

## plot data
# 提取样本点坐标
plot_data <- data.frame({pcoa$point})[1:2]

# 提取列名，便于后面操作。
plot_data$ID <- rownames(plot_data)
names(plot_data)[1:2] <- c('PCoA1', 'PCoA2')

# eig记录了PCoA排序结果中，主要排序轴的特征值（再除以特征值总和就是各轴的解释量）
eig = pcoa$eig

#为样本点坐标添加分组信息
plot_data <- merge(plot_data, group, by = 'ID', all.x = TRUE)
head(plot_data)

plot_data <- extract(plot_data, group, c("Group", "Timepoint"), regex = "([A-Z]+)(\\d+)")
  

# figure1
ggplot(data = plot_data, aes(x=PCoA1, y=PCoA2)) +
  geom_point(size = 3, aes(shape = Group, color = Timepoint)) +
  scale_color_manual(values = c("#f8766d", "#00ba38","#cacaca","#619cff"))+
  scale_shape_manual(values = 14:24) +
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""))+
  geom_hline(yintercept=0, linetype=4) +    
  geom_vline(xintercept=0 ,linetype=4)+          
  theme_bw() +
  guides(fill = guide_legend(override.aes=list(shape=14)))
ggsave('pcoa1.pdf',width = 6,height = 4)

# figure2
ggplot(data = plot_data, aes(x=PCoA1, y=PCoA2, fill=bacteria)) +
  geom_point(shape = 21,color = 'black',size=4) +
  scale_fill_gradient(low = '#f2fe32',high = '#180f7c')+
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""))+
  geom_hline(yintercept=0, linetype=4) +    
  geom_vline(xintercept=0 ,linetype=4)+          
  theme_few()+
  theme(legend.title = element_blank(),
        legend.position = c(0.8, 0.15),
        legend.direction = "horizontal")
ggsave('pcoa2.pdf',width = 4.5,height = 4)
        
sample_info