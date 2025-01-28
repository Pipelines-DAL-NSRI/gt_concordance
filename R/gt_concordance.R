#' Concordance analysis
#' @description
#' This package identifies overlapping data (usually samples) between two files 
#' and calculates the total number of matching data (usually genotypes). 
#' It creates both an output csv file and a barplot. 
#' @param file1 is either an xlsx or csv file
#' @param file2 is either an xlsx or csv file
#' @examples concordance("genotypes1.xlsx", "genotypes2.xlsx")
#' @export

concordance <- function(file1, file2){
   
   # check if file exists
   # read in the file
   if(!file.exists(file1)){
      report::report("First file does not exist in the working directory")
      stop()
   } else {
      if(tools::file_ext(file1) == "csv"){
         file1 <- readr::read_csv(file1)
      } else if(tools::file_ext(file1) == "xlsx"){
         file1 <- readxl::read_excel(file1)
      }
   }
      
   if(!file.exists(file2)){
      report::report("Second file does not exist in the working directory")
      stop()
   } else {
      if(tools::file_ext(file2) == "csv"){
         file2 <- readr::read_csv(file2)
      } else if(tools::file_ext(file2) == "xlsx"){
         file2 <- readxl::read_excel(file2)
      }
   }
   
   file_list <- list(file1, file2)
   
   # check intersecting values
   # samples are in the first column
   
   overlaps <- intersect(file1[,1], file2[,1])
   overlaps <- as.character(unlist(overlaps))
   
   # transpose
   file_list <- lapply(
      file_list,
      function(x){
         library(dplyr)
         data.frame(t(x)) %>%
            janitor::row_to_names(row_number = 1)
      }
   )
   
   #subset
   file_list <- lapply(
      file_list,
      function(x){
            subset(x, select = overlaps)
            ID <- row.names(x)
         data.frame(ID, x)
            data.frame(t(x))
            markers <- rownames(x)
            data.frame(markers, x)
      }
   )
   
   # merge 
   merged <- file_list %>% purrr::reduce(full_join, by= "markers")
   ID <- merged$markers
   
   merged <- lapply(
      merged, 
      function(x){
         gsub(pattern = "|", replacement = "/", x = x, fixed = TRUE)}
      )
   
   merged <- as.data.frame(merged)
   names(merged) <-  sub('^X', '', names(merged))
   
   # subset based on those with .x and .y
   # first
   overlap1 <- merged %>% select(dplyr::ends_with(".x"))
   overlap2 <- merged %>% select(dplyr::ends_with(".y"))
   
   for_conc <- QurvE::zipFastener(overlap1, overlap2, along = 2)
   for_conc2 <- data.frame(ID, for_conc)
   for_conc2[is.na(for_conc2)] <- "NA"
   names(for_conc2) <-  sub('^X', '', names(for_conc2))
   
   concordance <- bind_cols(for_conc2 %>%
                               tidyr::gather(var, val, -matches("(.x$|ID)")) %>%
                               select(ID,val), for_conc2 %>%
                               tidyr::gather(var2, val2, -matches("(*.y$|ID)")) %>%
                               select(val2)) %>%
      add_count(ID) %>%
      group_by(ID) %>%
      summarise(
         Total = paste((ncol(for_conc2) - 1)/2),
         Concordant = paste(sum(val == val2)),
         Discordant = paste(((ncol(for_conc2) - 1)/2)- sum(val == val2))
      ) %>%
      left_join(for_conc2, by = c("ID" = "ID")) 
   
   
   readr::write_csv(concordance, file = "concordance.csv")
   
   pivot <- concordance[,1]
   pivot2 <- concordance[,3:4]
   pivot <- data.frame(pivot, pivot2)
   pivot <- pivot %>%
      tidyr::pivot_longer(!ID,
                   names_to = 'Condition',
                   values_to = 'Count'
      )
   
   # select 
   Count <- as.integer(pivot$Count)
   rsID <- pivot$ID
   Condition <- pivot$Condition
   
   visual <- data.frame(rsID, Count, Condition)
   
   library(ggplot2)
   
   visual %>%
      arrange(Count) %>%
      arrange(Condition) %>%
      mutate(rsID = forcats::fct_inorder(rsID)) %>%
      ggplot(aes(fill = Condition, x= rsID, y=Count)) +
      geom_bar(position = "stack", stat = "identity") +
      theme(
         axis.text.x = element_text(
            angle = 90,
            vjust = .3), 
         panel.background = element_blank()) +
      scale_fill_manual(values= c("#1ca7ec",
                                  "#DC143C"
      ))
   
   ggsave(filename = "concordance_plot.png", width = 12, height = 4, dpi = 600)
   
      }
