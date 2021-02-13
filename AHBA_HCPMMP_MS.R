library(stringr)
library(broom)
library(dplyr)
library(multcomp)
library(ppcor)
library(knitr)
library(kableExtra)
library(tidyverse)
library(ggplot2)
library(tidyr)
library(ggrepel)
library(ggpmisc)
library(magrittr)
library(data.table)

# load AHBA gene
df_AHBA=read.csv(
  'D:\\cye_py_code\\connectome\\AHBA_HCPMMP379_nc006.csv',
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

mat_AHBA = as.matrix(df_AHBA)

# load decoupling
main_path='D:\\Image_data\\XWMS\\GSP_csv\\'
df_sub=read.csv(str_c(main_path,'\\participants.csv'),header = TRUE,row.names=1) 
df_sub[,'edss'] <- as.numeric(df_sub[,'edss'])
df_sub[,'dd'] <- as.numeric(df_sub[,'dd'])
df_sub[,'pasat_3'] <- as.numeric(df_sub[,'pasat_3'])
df_sub[,'pasat_2'] <- as.numeric(df_sub[,'pasat_2'])
df_sub[,'head_t2_lesion_volume'] <- as.numeric(df_sub[,'head_t2_lesion_volume'])
df_sub[,'spinal_t2_lesion_volume'] <- as.numeric(df_sub[,'spinal_t2_lesion_volume'])
df_sub <- data.frame(names = row.names(df_sub), df_sub)

region_type = '' # choose between '' or 'cortical360_'
ICN_type = 'Yeo' # choose one of the following: 'Yeo', 'Cole'
activity_type = 'BOLD' # choose one of the following: 'BOLD', 'zALFF'
Yeo_overlap = '0.6'
FDR_correction = TRUE
Tract_type = 'meanlength' # choose one of the following: 'level-participant_connectome', 'meanlength'
csv_path = str_c(main_path,'CIS_MS_NC_NMO\\BOLD_4D_',Tract_type,'_',Yeo_overlap)

df_sub <- df_sub %>% dplyr::filter(!names %in% 'sub-nc011') 
df_sub <- df_sub %>% dplyr::filter(!names %in% 'sub-nc039') 
df_sub <- df_sub %>% dplyr::filter(!names %in% 'sub-cis002') 
df_sub <- df_sub %>% dplyr::filter(!names %in% 'sub-cis010') 
df_sub <- df_sub %>% dplyr::filter(!names %in% 'sub-cis015')
df_sub <- df_sub %>% dplyr::filter(!names %in% 'sub-ms015') 
df_sub <- df_sub %>% dplyr::filter(!names %in% 'sub-nmo019') 
df_sub <- df_sub %>% dplyr::filter(!names %in% 'sub-nmo002') 
csv_path = str_c(main_path,'CIS_MS_NC_NMO\\BOLD_4D_',Tract_type,'_',Yeo_overlap)
df_parcel=read.csv(
  str_c(csv_path, '\\s_deCoupIdx_node_individual_df.csv'),
  header = TRUE,
  row.names = 1,
  check.names = TRUE
)

df_sub = df_sub %>% arrange(group, age) %>% filter(!(group %in% 'NMO' & age >51))
df_parcel = merge(df_sub, df_parcel, by.x = 'names', by.y = 'subname')
df_parcel_NC_MS = df_parcel %>% dplyr::filter(group %in% c('NC','MS'))

df_parcel_individual_NC = df_parcel %>% dplyr::filter(group == 'NC') %>% dplyr::select('Visual1.04_L.Ctx':'Brain.Stem') 
df_parcel_individual_MS = df_parcel %>% dplyr::filter(group == 'MS') %>% dplyr::select('Visual1.04_L.Ctx':'Brain.Stem') 
df_parcel_individual_NMO = df_parcel %>% dplyr::filter(group == 'NMO') %>% dplyr::select('Visual1.04_L.Ctx':'Brain.Stem')

# T-statistics
vec_T_MS = vector(length = dim(df_parcel_individual_NC)[2])

for (i in 1:dim(df_parcel_individual_NC)[2]) {
  T = t.test(df_parcel_individual_MS[,i], df_parcel_individual_NC[,i],  var.equal=FALSE)
  vec_T_MS[i] = T$statistic
}

# F-statistics
ROInames= colnames(df_parcel_individual_NC)
vec_F_MS = vector(length = dim(df_parcel_individual_NC)[2])
count_c = 0
for (i in ROInames) {
  count_c = count_c+1
  stat_dist = aov(get(i) ~ group , data = df_parcel_NC_MS)
  vec_F_MS[count_c] <-  unlist(summary(stat_dist))["F value1"]
}

# if only for cortical, make sure #parcel=360
if (dim(mat_AHBA)[1] == 360){
  vec_T_MS = vec_T_MS[1:360]
  vec_F_MS = vec_F_MS[1:360]
}

# remove unstable nodes
df_node_rand=read.csv(str_c(csv_path,'\\s_deCoupIdx_node_rand_df.csv'), header = TRUE, row.names=1)
df_node_empi=read.csv(str_c(csv_path,'\\s_deCoupIdx_node_empirical_df.csv'), header = TRUE, row.names=1)
# if only for cortical, make sure #parcel=360
if (dim(mat_AHBA)[1] == 360){
  df_node_rand = df_node_rand[,1:360]
  df_node_empi = df_node_empi[,1:360]
}
df_node_rand2=df_node_rand %>% gather(key = nodelabel,value = value)
df_node_empi2=df_node_empi %>% gather(key = nodelabel,value = value)
node_quantile_ls = list()
for (i in colnames(df_node_empi)) {
  percent_name = str_c(i,'_percent')
  node_quantile_ls <- c(node_quantile_ls, ecdf(df_node_rand[,i])(df_node_empi[1,i]))
}

df_node_quantile = as.data.frame(node_quantile_ls, col.names = colnames(df_node_empi))
node_fail_surrogate <- colnames(df_node_quantile)[df_node_quantile >= 0.05]

# for (i in node_fail_surrogate) {
#   ROInames <<- ROInames[ROInames != i]  
# }

mat_AHBA = mat_AHBA[df_node_quantile >= 0.05,]
vec_T_MS = vec_T_MS[df_node_quantile >= 0.05]
vec_F_MS = vec_F_MS[df_node_quantile >= 0.05]



# spatial correlation
mat_correlation = matrix(nrow = dim(mat_AHBA)[2], ncol = 2)
for (i in 1:dim(mat_AHBA)[2]) {
  pcor_stats = cor.test(x = vec_T_MS,
                        y = mat_AHBA[, i],
                        method = "spearman")
  mat_correlation[i,1] = pcor_stats$estimate
  mat_correlation[i,2] = pcor_stats$p.value
}

df_correlation = as.data.frame(mat_correlation)
colnames(df_correlation) = c('rho', 'p-value')
rownames(df_correlation) = colnames(mat_AHBA)

df_correlation$`p-value` = p.adjust(df_correlation$`p-value`, method = "bonferroni")
df_correlation_sig = df_correlation %>% dplyr::filter(rho > 0.45 | rho < -0.45)

# histogram plot
ggplot(df_correlation, aes(x=rho)) + 
  geom_histogram(color="black", fill="white",binwidth=0.01)+ 
  geom_vline(aes(xintercept=0.48),color="blue", linetype="dashed", size=1)+ 
  labs(x = "correlation coefficient", y = 'frequency') + theme_light()
ggsave(file="gene_histo.png", width=6, height=4, dpi=300)

# correlation plot
df_plot = data.frame(decoupling = vec_T_MS, gene = mat_AHBA[,'CNR1'])
prog_col="#619CFF"
ggplot(df_plot, aes(x = decoupling, y = gene)) +
  geom_point(color = prog_col, size = 2) +
  geom_smooth(aes(color = prog_col), method = lm, se = TRUE)+
  scale_color_manual(values = prog_col)  +
  scale_fill_manual(values = prog_col)+ labs(x = "Decoupling Index", y = expression("gene expression")) +
  ggtitle("Transcriptome:CNR1") + theme_light() +
  theme(plot.title = element_text(hjust = 0.5),text = element_text(size=20),legend.position = "none") +ylim(0.4, 0.85)
ggsave(file="CNR1_corr.png", width=6, height=4, dpi=300)
