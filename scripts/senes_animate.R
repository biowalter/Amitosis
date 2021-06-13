library(ggplot2)
library(gganimate)
library(ggsci)
library(transformr)

pal = colorRampPalette(pal_jama()(7))(length(levels(reduced_all$IRS0)))
## set y limit to improve visualization (for coord_cartesian)
roof = reduced_all %>% 
  as_tibble() %>%
  filter(Alleles_n > 0 & Alleles_n < k)
roof = max(roof$P_dist) * 1.2

## change default animation par
options(gganimate.nframes = length(unique(reduced_all$GEN)) + 1)
options(gganimate.fps = 10)

## animated faceted plot
anim <- ggplot(reduced_all, aes(x = Alleles_n, P_dist, color = IRS0)) +
  geom_area(aes(group = GEN), fill = "gray") +
  xlim(c(0, k)) +
  scale_color_manual(values = pal) +
  coord_cartesian(ylim=c(0, roof)) +
  facet_wrap(~ IRS0, ncol = 3) +
  labs(title="", x = expression("x copies"), 
       y = expression("P"['(x)'])) +
  theme_bw(base_size = 8) + 
  theme(plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = 'white smoke'),
        axis.title.x = element_text(color = "black", size = 12, face = "plain", vjust=-1),
        axis.title.y = element_text(color = "black", size = 12, face = "plain", vjust= 0.5, 
                                    angle = 90, margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text.x = element_text(color = "black", size = 8, face = "plain"),
        axis.text.y = element_text(color = "black", size = 8, face = "plain"),
        legend.position = "none") +
  # gganimate specific block
  transition_states(
    GEN,
    transition_length = 2,
    state_length = 2
  ) + 
  ggtitle('Generation {frame} of {nframes}')

animate(anim, height = 2.5, width = 5, units = "in", res = 150)

## Save
anim_save("senes.gif", path = "~/1_Workspace/1_PostDOC/Isabel_Lab/presentation/images/")