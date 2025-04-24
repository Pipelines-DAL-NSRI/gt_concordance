#' Concordance analysis
#' @description
#' This package identifies overlapping data (usually samples) between two files 
#' and calculates the total number of matching data (usually genotypes). 
#' It creates both an output csv file and a barplot. 
#' @param file1 is either an xlsx or csv file
#' @param file2 is either an xlsx or csv file
#' @param haplotypes (required) if set to TRUE, will no longer change the allele arrangement. 
#' @examples concordance("genotypes1.xlsx", "genotypes2.xlsx", haplotypes = FALSE)
#' @import readr
#' @import readxl
#' @import tools
#' @import dplyr
#' @import janitor
#' @import purrr
#' @import tibble
#' @import tidyselect
#' @import QurvE
#' @import tidyr
#' @import ggplot2
#' @import forcats
#' @export

concordance <- function(file1, file2, haplotypes = FALSE){
   
   if (!require("pacman")){
      install.packages("pacman")
   }
   
   pacman::p_load(readr, readxl, tools, dplyr, janitor, purrr, tibble, tidyselect, QurvE, tidyr, ggplot2, forcats, install = TRUE)
   
   # check if file exists
   # read in the file
   if(!file.exists(file1)){
      stop("First file does not exist in the working directory")
   } else {
      if(tools::file_ext(file1) == "csv"){
         file1 <- readr::read_csv(file1)
      } else if(tools::file_ext(file1) == "xlsx"){
         file1 <- readxl::read_excel(file1)
      }
   }
      
   if(!file.exists(file2)){
      stop("Second file does not exist in the working directory")
   } else {
      if(tools::file_ext(file2) == "csv"){
         file2 <- readr::read_csv(file2)
      } else if(tools::file_ext(file2) == "xlsx"){
         file2 <- readxl::read_excel(file2)
      }
   }
   
   file1 <- file1 %>%
      rename(Ind = 1)
   file2 <- file2 %>%
      rename(Ind = 1)
   
   file_list <- list(file1, file2)
   
   # check intersecting values
   # samples are in the first column
   overlaps <- as.list(intersect(file1$Ind, file2$Ind))

   
   #subset
   file_list2 <- lapply(
      file_list,
      function(x){
         x[x$Ind %in% overlaps, ]
      }
   )
   
   # transpose
   file_list3 <- lapply(
      file_list2,
      function(x){
         library(dplyr)
         data.frame(t(x)) %>%
            janitor::row_to_names(row_number = 1) %>%
            tibble::rownames_to_column(., var = "markers") 
      }
   )
   
   # rearrange column order
   overlaps <- as.character(overlaps)
   markers1 = file_list3[[1]]$markers
   markers2 = file_list3[[2]]$markers
   
   file_list4 <- lapply(
      file_list3,
      function(x){
         #markers = x$markers
         relocate(x, any_of(overlaps)) 
      }
   )

   
   #re-add marker column
   file_list4[[1]]$markers = markers1
   file_list4[[2]]$markers = markers2
 
   
   # merge 
   merged <- file_list4 %>% purrr::reduce(full_join, by= "markers")
   ID <- merged$markers
   
   merged <- lapply(
      merged, 
      function(x){
         gsub(pattern = "|", replacement = "/", x = x, fixed = TRUE)}
   )
   
   merged <- as.data.frame(merged)
   
   if(haplotypes == TRUE){
      
      print("Assuming the data are haplotypes.")

   } else if(haplotypes == FALSE){
      merged <- merged %>% mutate(across(tidyselect::everything(), ~ case_when(
         . == "A" ~ "A/A",
         . == "T" ~ "T/T",
         . == "C" ~ "C/C",
         . == "G" ~ "G/G",
         . == "T/C" ~ "C/T",
         . == "T/G" ~ "G/T",
         . == "T/A" ~ "A/T",
         . == "G/A" ~ "A/G",
         . == "C/G" ~ "G/C",
         . == "C/A" ~ "A/C",
         TRUE ~ .x)))
   } else {
      stop("Parameter haplotype is required.")
   }
   
   
   # subset based on those with .x and .y
   # rearrange based on order
   overlap1 <- merged %>% select(dplyr::ends_with(".x"))
   overlap2 <- merged %>% select(dplyr::ends_with(".y"))
   
   for_conc <- QurvE::zipFastener(overlap1, overlap2, along = 2)
   for_conc2 <- data.frame(ID, for_conc)
   for_conc2[is.na(for_conc2)] <- "N"
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
         Incomparable = paste(sum(val == "N") + sum(val2 == "N")),
         Concordant = paste(sum(val == val2)),
         Discordant = paste(((ncol(for_conc2) - 1)/2)- sum(val == val2) - (sum(val == "N") + sum(val2 == "N")))
      ) %>%
      left_join(for_conc2, by = c("ID" = "ID")) 
   
   
   readr::write_csv(concordance, file = "concordance.csv")
   
   pivot <- concordance[,1]
   pivot2 <- concordance[,3:5]
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
                                  "#fb7a8e",
                                  "#1f2f98"
      ))
   
   ggsave(filename = "concordance_plot.png", width = 12, height = 4, dpi = 600)
   
      }
