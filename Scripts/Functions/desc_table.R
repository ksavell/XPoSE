# dependencies for desc_table
# library(janitor) 
# library(flextable)

#' Function to create a descriptive table
#'
#' @param var1 rows
#' @param var2 columns
#' @param rowname 
#' @param colname 
#' @data data frame with variables 1 and 2
#'
#' @return
#' @export
desc_table <- function(var1,var2,rowname,colname,data){
  
  store <- data %>%
    tabyl({{var1}},{{var2}},show_na = TRUE) 
  
  #chires <- chisq.test(store) # used to calculate p-values between var1 and var2
  
  store <- store %>%      # tabulate if participant has been homeless and proportions by gender category
    adorn_totals("col") %>% # total in each row
    adorn_totals("row") %>% # total in each row
    adorn_percentages("row") %>%
    adorn_pct_formatting(digits = 1) %>%
    adorn_ns(position = "front") %>% 
    adorn_title(                          # adjust titles
      row_name = rowname,
      col_name = colname,
      placement = "combined")  # this is necessary to print as image 
  
  # store$`p-value` <- c(rep('',length(store$Total)-1),round(chires$p.value,5)) # this is for p-values
  store %>% flextable() %>% width(width=7) %>% theme_zebra()
  
}
