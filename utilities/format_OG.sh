# Load the required packages
library(dplyr)
library(tidyr)

# You should remove from Orthogroups.txt the unassigned Orthogroups

a <- read.table("assigned_ortho.txt", sep = " ", stringsAsFactors= FALSE, fill = TRUE)


# check how many columns a has, exclude the first column
cols_to_pivot <- paste0("V", 2:529)



# Use tidyr::pivot_longer() function to transform the data frame
df_transformed <- tidyr::pivot_longer(a, cols = all_of(cols_to_pivot), values_drop_na = TRUE) %>%
  select(ortho = V1, genes = value)


df_transformed <- distinct(df_transformed)

# Convert the resulting tibble to a data frame with the desired column names
df_transformed <- as.data.frame(df_transformed)



write.table(df_transformed, file ="final_map.txt", row.names = F, quote = F, sep ="\t")



#### OUtSIDE R
#### Remove empty second columns, remove colon after OG and transform back the ___ to :
awk '$2 != ""' final_map.txt | sed 's/:\t/\t/g' | sed 's/___/:/g' > final
mv final final_map.txt
