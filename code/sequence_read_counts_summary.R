library("here")
library("tidyverse")
library("scales")
# library("readxl")
# library("writexl")

################################
## get read counts

    reads <- list()
    
    # prior to demultiplexing
    
        reads$pre_demultiplexing <- list(
            `16S` = read_delim(here("out", "sequence_counts", "pre_demultiplexing", "16S.txt"), col_names = c("file","reads"), delim = " "),
            `ITS` = read_delim(here("out", "sequence_counts", "pre_demultiplexing", "ITS.txt"), col_names = c("file","reads"), delim = " ")
        )
        
    # after demultiplexing
        
        reads$post_demultiplexing <- list(
            `16S` = read_delim(here("out", "sequence_counts", "post_demultiplexing", "16S.txt"), col_names = c("sample","reads"), delim = " "),
            `ITS` = read_delim(here("out", "sequence_counts", "post_demultiplexing", "ITS.txt"), col_names = c("sample","reads"), delim = " ")
        )
        
    # after demultiplexing and filtering
        
        reads$post_demultiplexing_noNs <- list(
            `16S` = read_delim(here("out", "sequence_counts", "post_demultiplexing_noNs", "16S.txt"), col_names = c("sample","reads"), delim = " "),
            `ITS` = read_delim(here("out", "sequence_counts", "post_demultiplexing_noNs", "ITS.txt"), col_names = c("sample","reads"), delim = " ")
        )
        
    # after cutadapt
    
        reads$post_cutadapt <- list(
            `16S` = read_delim(here("out", "sequence_counts", "post_cutadapt", "16S.txt"), col_names = c("sample","reads"), delim = " "),
            `ITS` = read_delim(here("out", "sequence_counts", "post_cutadapt", "ITS.txt"), col_names = c("sample","reads"), delim = " ")
        )
        
    # post filter and trim
    
        reads$post_filterAndTrim <- list(
            `16S` = read_delim(here("out", "sequence_counts", "post_filterAndTrim", "16S.txt"), col_names = c("sample","reads"), delim = " "),
            `ITS` = read_delim(here("out", "sequence_counts", "post_filterAndTrim", "ITS.txt"), col_names = c("sample","reads"), delim = " ")
        )
        
    # post dada
        
        reads$post_dada <- list(
            `16S` = read_delim(here("out", "sequence_counts", "post_dada", "16S.txt"), delim = " "),
            `ITS` = read_delim(here("out", "sequence_counts", "post_dada", "ITS.txt"), delim = " ")
        )
        
    # post chimeras
        
        reads$post_remove_chimeras <- list(
            `16S` = read_delim(here("out", "sequence_counts", "post_chimeras", "16S.txt"), delim = " "),
            `ITS` = read_delim(here("out", "sequence_counts", "post_chimeras", "ITS.txt"), delim = " ")
        )
        
    # post export to phyloseq (should be the same)
        
        reads$post_export_phyloseq <- list(
            `16S` = read_delim(here("out", "sequence_counts", "post_phyloseq", "16S", "phyloseq.txt"), delim = " "),
            `16S_cleaned` = read_delim(here("out", "sequence_counts", "post_phyloseq", "16S", "phyloseq_cleaned.txt"), delim = " "),
            `ITS` = read_delim(here("out", "sequence_counts", "post_phyloseq", "ITS", "phyloseq.txt"), delim = " ")
        )
    
    # post normalize (should be the same)
    
        reads$post_export_phyloseq_normalized <- list(
            `16S` = read_delim(here("out", "sequence_counts", "post_normalize", "16S.txt"), delim = " "),
            `ITS` = read_delim(here("out", "sequence_counts", "post_normalize", "ITS.txt"), delim = " ")
        )

##################################
# code below needs to be revised #
##################################
## collate data
        
    stage.order <- c("before demultiplexing", "after demultiplexing", "after filtering Ns", "after cutadapt", "after filterAndTrim", "after dada", "after removing chimeras", "after exporting to phyloseq", "after normalizing")
        
    # gather reads by sample
        
        reads_by_sample <- bind_rows(
            list(
                "after demultiplexing" = bind_rows(list(`16S` = reads$post_demultiplexing$`16S`, `ITS` = reads$post_demultiplexing$`ITS`), .id = "dataset"),
                "after filtering Ns" = bind_rows(list(`16S` = reads$post_demultiplexing_noNs$`16S`, `ITS` = reads$post_demultiplexing_noNs$`ITS`), .id = "dataset"),
                "after cutadapt" = bind_rows(list(`16S` = reads$post_cutadapt$`16S`, `ITS` = reads$post_cutadapt$`ITS`), .id = "dataset"),
                "after filterAndTrim" = bind_rows(list(`16S` = reads$post_filterAndTrim$`16S`, `ITS` = reads$post_filterAndTrim$`ITS`), .id = "dataset"),
                "after dada" = bind_rows(list(`16S` = reads$post_dada$`16S`, `ITS` = reads$post_dada$`ITS`), .id = "dataset"),
                "after removing chimeras" = bind_rows(list(`16S` = reads$post_remove_chimeras$`16S`, `ITS` = reads$post_remove_chimeras$`ITS`), .id = "dataset"),
                # should be the same
                    "after exporting to phyloseq" = bind_rows(list(`16S` = reads$post_export_phyloseq$`16S`, `ITS` = reads$post_export_phyloseq$`ITS`), .id = "dataset"),
                # should be the same
                    "after normalizing" = bind_rows(list(`16S` = reads$post_export_phyloseq_normalized$`16S`, `ITS` = reads$post_export_phyloseq_normalized$`ITS`), .id = "dataset")
            ), .id = "stage") %>%
        select(dataset, everything())
        
    # calculate totals and add pre-demultiplexing data

        read_totals <- bind_rows(
                reads$pre_demultiplexing %>% mutate(stage = "before demultiplexing"),
                reads_by_sample %>%
                    group_by(dataset, stage) %>%
                    summarize(reads = sum(reads), .groups = "drop")) %>%
            mutate(stage = factor(stage, levels = stage.order)) %>% 
            arrange(dataset, stage)
            
        # add percentage of reads remaining
        read_totals <- read_totals %>% group_by(dataset) %>%
            mutate(fraction = reads/reads[1]) %>%
            ungroup() %>%
            select(dataset, stage, reads, fraction)
        
        read_totals

        # plot total reads in each dataset through the pipeline
        read_totals %>% ggplot(aes(x = stage, y= reads)) + geom_bar(stat = "identity") + facet_grid(~ dataset) +
            theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust=0)) +
            labs (x = "Pipeline stage", y= "Total remaining reads") +
            ggtitle("Total reads per dataset") +
            scale_y_continuous(labels = label_number(suffix = " M", scale = 1e-6)) + # show y-axis in millions
            theme(plot.margin = margin(t = 10, r = 50, b = 10, l = 10, unit = "pt"))
        ggsave(here("figures", "reads_per_dataset_total.pdf"), width = 6, height = 5)

        # plot fraction of reads in each dataset through the pipeline
        read_totals %>% ggplot(aes(x = stage, y= fraction)) + geom_bar(stat = "identity") + facet_grid(~ dataset) +
            theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust=0)) +
            labs (x = "Pipeline stage", y= "Fraction of initial reads remaining") +
            ggtitle("Total reads per dataset") +
            theme(plot.margin = margin(t = 10, r = 50, b = 10, l = 10, unit = "pt"))
        ggsave(here("figures", "reads_per_dataset_fraction.pdf"), width = 6, height = 5)
        
    # reads_by_sample
        
        # plot total reads per sample through the pipeline
        reads_by_sample %>%
            mutate(stage = factor(stage, levels = stage.order)) %>% 
            arrange(dataset, stage) %>%
            ggplot(aes(x = stage, y= reads)) + geom_boxplot() + facet_grid( ~ dataset) +
            theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust=0)) +
            labs (x = "Pipeline stage", y= "Remaining reads per sample") +
            scale_y_continuous(labels = label_number(suffix = " K", scale = 1e-3)) + # show y-axis in thousands
            ggtitle("Total reads per sample") +
            theme(plot.margin = margin(t = 10, r = 50, b = 10, l = 10, unit = "pt"))
        ggsave(here("figures", "reads_per_sample_total.pdf"), width = 6, height = 5)
        
        # plot per-sample read dropout through the pipeline
        reads_by_sample %>%
            mutate(stage = factor(stage, levels = stage.order)) %>% 
            group_by(dataset, sample) %>%
            mutate(read_loss_fraction = 1 - reads/dplyr::lag(reads, n=1, order_by=stage)) %>%
            ggplot(aes(x = stage, y= read_loss_fraction)) + geom_boxplot() + facet_grid( ~ dataset) +
            theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust=0)) +
            labs (x = "Pipeline stage", y= "Fraction of remaining reads dropping out") +
            ggtitle("Read dropout per sample") +
            theme(plot.margin = margin(t = 10, r = 50, b = 10, l = 10, unit = "pt"))
        ggsave(here("figures", "reads_per_sample_dropout_fraction.pdf"), width = 6, height = 5)
