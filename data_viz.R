library(cowplot)
library(rstatix)
library(tidyverse)

meta_fixer <- function(meta_fp){
  mycols <- read_tsv(meta_fp)
  mycols <- colnames(mycols)
  meta <- read_tsv(meta_fp, skip = 2)
  colnames(meta) <- mycols
  colnames(meta)[1] <- 'sampleid'
  return(meta)
}

diet_labs = c('HF/HF' = 'High Fat/High Fiber',
              'HF/LF' = 'High Fat/Low Fiber',
              'LF/HF' = 'Low Fat/High Fiber',
              'LF/LF' = 'Low Fat/Low Fiber')

meta <- meta_fixer('./lacto_sum_metadata.tsv')
faith <- read_tsv('./cda/faith_pd.tsv')
uu <- read_tsv('./cda/unwe`')

meta %>% 
  filter(diet != 'Chow') %>% 
  mutate(diet = diet_labs[diet]) %>% 
  left_join(faith) -> meta

meta %>% 
  group_by(day_post_inf, diet) %>% 
  mutate(counts = n()) %>% 
  filter(counts > 1) %>% 
  ungroup() %>% 
  group_by(day_post_inf) %>% 
  kruskal_test(faith_pd ~ diet) %>% 
  adjust_pvalue()

meta %>% 
  filter(!is.na(faith_pd)) %>% 
  group_by(day_post_inf, diet) %>% 
  mutate(counts = n()) %>% 
  filter(counts > 1) %>% 
  ungroup() %>% 
  group_by(day_post_inf) %>% 
  dunn_test(faith_pd ~ diet) -> dunn_faith

meta %>% 
  filter(!is.na(faith_pd)) %>% 
  group_by(day_post_inf, diet) %>% 
  summarise(avg_faith = mean(faith_pd)) -> diet_summs


dunn_faith %>% 
  merge(diet_summs, by.x = c('group1', 'day_post_inf'),
                    by.y = c('diet', 'day_post_inf')) %>% 
  rename(avg_faith1 = avg_faith) %>% 
  merge(diet_summs, by.x = c('group2', 'day_post_inf'),
                    by.y = c('diet', 'day_post_inf')) %>% 
  rename(avg_faith2 = avg_faith) %>% 
  mutate(faith_perc_diff = (avg_faith1 - avg_faith2)/avg_faith1) -> dunn_faith


meta %>%
  ggplot(aes(x = day_post_inf, y = faith_pd)) +
    geom_boxplot(alpha = 0, outlier.shape = NA, 
                 aes(group = day_post_inf)) +
    geom_jitter(width = 0.1, height = 0, alpha = 0.6) +
    geom_smooth(se=FALSE) +
    facet_wrap(~diet, ncol = 1, strip.position = 'left') +
    theme_bw(base_size = 18) +
    theme(strip.text.y.left = element_text(angle = 0)) + 
    labs(x = '',
         y = "Faith's PD") -> alpha_plot


dunn_faith %>% 
  unite(groups, group1, group2, sep = ' : ') %>% 
  ggplot(aes(x = day_post_inf, y = groups, fill = faith_perc_diff)) +
    geom_tile(color = 'black') +
    geom_text(aes(label = p.adj.signif)) +
    scale_fill_gradient2() +
    theme_bw(base_size=16) +
    labs(x = 'Days Post Infection',
         y = 'Groups') -> dunn_plots


plot_grid(alpha_plot, dunn_plots,
          ncol = 1,
          axis = 'tblr',
          align = 'hv',
          rel_heights = c(1, 0.35))
