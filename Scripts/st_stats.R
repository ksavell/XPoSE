# sample tag assignment stats

# Subset to nuclei assigned to SampleTag08, which is rat HC-1
cells <- WhichCells(all, expression = Sample_tag == "SampleTag08_mm")

# Compare correct vs incorrect tag reads
correct <- FetchData(all, vars = "SampleTag08_reads")[cells, 1]
incorrect <- rowMeans(FetchData(all, vars = c("SampleTag02_reads", "SampleTag04_reads", "SampleTag06_reads"))[cells, ])

wilcox.test(correct, incorrect, paired = TRUE, alternative = "greater")

# Subset to nuclei assigned to SampleTag04, which is rat HC-2
cells <- WhichCells(all, expression = Sample_tag == "SampleTag04_mm")

# Compare correct vs incorrect tag reads
correct <- FetchData(all, vars = "SampleTag04_reads")[cells, 1]
incorrect <- rowMeans(FetchData(all, vars = c("SampleTag02_reads", "SampleTag06_reads", "SampleTag08_reads"))[cells, ])

wilcox.test(correct, incorrect, paired = TRUE, alternative = "greater")


# Subset to nuclei assigned to SampleTag02, which is rat HC-3
cells <- WhichCells(all, expression = Sample_tag == "SampleTag02_mm")

# Compare correct vs incorrect tag reads
correct <- FetchData(all, vars = "SampleTag02_reads")[cells, 1]
incorrect <- rowMeans(FetchData(all, vars = c("SampleTag04_reads", "SampleTag06_reads", "SampleTag08_reads"))[cells, ])

wilcox.test(correct, incorrect, paired = TRUE, alternative = "greater")

# Subset to nuclei assigned to SampleTag06, which is rat HC-4
cells <- WhichCells(all, expression = Sample_tag == "SampleTag06_mm")

# Compare correct vs incorrect tag reads
correct <- FetchData(all, vars = "SampleTag06_reads")[cells, 1]
incorrect <- rowMeans(FetchData(all, vars = c("SampleTag02_reads", "SampleTag04_reads", "SampleTag08_reads"))[cells, ])

wilcox.test(correct, incorrect, paired = TRUE, alternative = "greater")
