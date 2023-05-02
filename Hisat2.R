Flagstat_df = data.frame("Sample" = c("SRR6853319","SRR6853320","SRR6853321","SRR6853322","SRR6853323",
                                      "SRR6853324","SRR6853325","SRR6853326","SRR6853327","SRR6853328","SRR6853329","SRR6853330",
                                      "SRR6853331","SRR6853332","SRR6853333","SRR6853334","SRR6853335","SRR6853336"),
                         "Percent_mapped" = c(73.59, 75.36, 70.02, 65.99, 65.09, 71.43, 73.35, 75.11, 74.62, 72.42, 68.66, 66.51, 75.61,
                                              66.42, 70.33, 70.88, 75.50, 76.00),
                         "Percent_properly_paired" = c(62.47, 66.19, 61.36, 58.12, 57.72, 61.63, 64.47, 64.65, 65.45, 63.75,
                                                       59.21, 59.00, 63.05, 55.53, 57.24, 57.27, 64.62, 65.89))

library(lattice)
install.packages("ggplot2")
library(ggplot2)
install.packages("gridExtra")
library(gridExtra)

Mapped_reads_plot = ggplot(Flagstat_df, aes(y=Percent_mapped,x=Sample)) + 
  theme(axis.text.x = element_text(angle = 90)) +
  geom_point() +
  ylim(55,80) +
  ylab("Reads Mapped (%)") + 
  ggtitle("HiSAT2 Mapped Reads")

Paired_reads_plot = ggplot(Flagstat_df, aes(y=Percent_properly_paired,x=Sample)) + 
  theme(axis.text.x = element_text(angle = 90)) +
  geom_point() +
  ylim(55,80) +
  ylab("Properly Paired Reads (%)") + 
  ggtitle("HiSAT2 Paired Reads")

grid.arrange(Mapped_reads_plot, Paired_reads_plot)
