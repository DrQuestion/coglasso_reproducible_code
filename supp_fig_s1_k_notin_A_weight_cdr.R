# This is the code to reproduce the Supplementary Figure S1 of the c-dependent
# weight of the contribution of the variables k \notin A in the coordinate 
# descent rule.

# Correct the following line accordingly. 
setwd("Your/output/path")

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# END OF REQUIRED INPUT FROM USER #
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

ggplot() + 
  xlim(0, 11) + ylim(-1, 1) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
  geom_function(fun = function(c) (1 - c)/(1 + c), linewidth = 1.2, color = "red") + 
  scale_x_continuous(limits = c(-0, 10), breaks = seq(0, 10)) +
  theme_bw() +
  theme(axis.line = element_blank(),                       
        panel.grid.major = element_line(color = "gray90"), 
        panel.grid.minor.x = element_blank()  ) +
  ylab("weight") + xlab("c value") 
ggsave("c_value_weight.svg")
