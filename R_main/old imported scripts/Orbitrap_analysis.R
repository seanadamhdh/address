library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)


#formula attribution for A12TH from icbm
#setwd("~/Desktop/orbit trap RZ/Processed_data")

A12_TH_orbitrap <- read.csv("~/Documents/GitHub/ADDRESS-adit_drainage_solute_source_control/data/Formula_attribution_A12TH_Feb2024.csv")

# fixing names: all sample Names have to be named Sample_NAME_1.csv; NAME must not contain "_"
names(A12_TH_orbitrap)[ which(str_detect(names(A12_TH_orbitrap),"Sample")&
                                !str_detect(names(A12_TH_orbitrap),"_1.csv"))]<-
  str_replace(names(A12_TH_orbitrap)[ which(str_detect(names(A12_TH_orbitrap),"Sample")&
                                              !str_detect(names(A12_TH_orbitrap),"_1.csv"))],
              ".csv","_1.csv")

# renaming (pulling sample "NAME" from original "Sample_NAME_1.csv")
names(A12_TH_orbitrap)[ which(str_detect(names(A12_TH_orbitrap),"Sample"))]<-
  strsplit(names(A12_TH_orbitrap)[which(str_detect(names(A12_TH_orbitrap),"Sample"))],"_")%>%
  map(2)%>%unlist

# long format -> 1 col sample id, 1 col intensities
A12_TH_orbitrap%>%pivot_longer(
  # Find all columns that start with 'A12'
  cols=names(A12_TH_orbitrap)[which(str_detect(names(A12_TH_orbitrap),"A12."))],
  names_to = "Sample_ID",
  values_to = "Intensity"
)->A12_TH_orbitrap_long

  

# Selecting the subset of relevant columns
A12_TH_orbitrap_subset <- A12_TH_orbitrap_long[c('Sample_ID', 'id', 'mz', 'diff', 'reference', 'formula', 'H.C', 'O.C', 'C', 'H', 'O', 'N', 'S', 'P', 'AI', 'AI.mod', 
                     'DBE', 'Aromatic', 'Highly.unsaturated', 'Unsaturated', 'Saturated', 'Intensity')]


# Calculate relative intensity
# grouped sums for each unique Sample_ID
A12_TH_orbitrap_subset%>%group_by(Sample_ID)%>%
  summarise(Intensity_sums=sum(Intensity,na.rm = T))->Intensity_sums

# joining by Sample_ID
left_join(A12_TH_orbitrap_subset,Intensity_sums,by="Sample_ID")%>%
  mutate(relative_Intensity=Intensity/Intensity_sums)->data_subset # rename to fit your workflow
 


#### END OF SEANS MODIFICATIONS ####
####################################

# Calculate NOSC for each row
data_subset$NOSC <- 4 - ((4 * data_subset$C + data_subset$H - 3 * data_subset$N - 2 * data_subset$O) / data_subset$C)

# Multiply and create new columns
cols_to_multiply <- c('mz', 'H.C', 'O.C', 'AI', 'AI.mod', 'DBE', 'Aromatic', 'Highly.unsaturated', 'Unsaturated', 'NOSC')
for(col in cols_to_multiply) {
  new_col_name <- paste(col, 'weighted', sep = '_')
  data_subset[[new_col_name]] <- data_subset[[col]] * data_subset$relative_intensity
}







#rename modified dataframes for revised code
modified_dataframes_A12TH <-   modified_dataframes

# Iterate over each dataframe in the list
for(i in seq_along(modified_dataframes)) {
  # Extract the name of the dataframe
  df_name <- names(modified_dataframes)[i]
  
  # Extract the dataframe
  df <- modified_dataframes[[i]]
  
  # Create a filename using the dataframe name
  filename <- paste0(df_name, ".csv")
  
  # Save the dataframe as a CSV file
  write.csv(df, filename, row.names = FALSE)
}

# Loop through the list of modified data frames
for(col_name in names(modified_dataframes)) {
  # Extract the data frame
  df <- modified_dataframes[[col_name]]

  
  # Plot the Van Krevelen diagram
  p <- ggplot(df, aes(x = O.C, y = H.C)) +
    geom_point(alpha = 0.5) +
    xlim(c(0, max(df$O.C, na.rm = TRUE))) +
    ylim(c(0, max(df$H.C, na.rm = TRUE))) +
    labs(title = paste("Van Krevelen Diagram for", col_name),
         x = "O/C Ratio",
         y = "H/C Ratio") +
    theme_bw()
  
  
  # Save the plot
  ggsave(filename = paste0("Van_Krevelen_", col_name, ".png"), width = 8, height = 6)
}

# Now, modified_dataframes list contains a dataframe for each 'Sample..' column

#looking at all peaks

# Loop through each dataframe in the list
for(i in seq_along(modified_dataframes)) {
  # Get the name of the current sample
  sample_name <- names(modified_dataframes)[i]
  
  # Get the current dataframe
  current_dataframe <- modified_dataframes[[sample_name]]
  
  # Create the ggplot using normalized relative intensity
  p <- ggplot(current_dataframe, aes(x = mz, y = relative_intensity)) +
    geom_line() +
    theme_minimal() +
    labs(title = paste("Mass Spectrum -", sample_name),
         x = "m/z", y = "Relative Intensity")
  
  # Convert to an interactive plotly object
  interactive_plot <- ggplotly(p)
  
  # Print the interactive plot
  print(interactive_plot)
  
  # Optionally, save the interactive plot to an HTML file
  html_file <- paste0("spectrum_", sample_name, ".html")
  htmlwidgets::saveWidget(interactive_plot, html_file)
}


# Load A12_all_top regular data to combine with orbitrap data of A12 top
setwd("~/Desktop/R_Anita_clean_Aug2023")

A12_all_top <- read_excel("A12_A13_allregsampling.xlsx",sheet="A12_all_top (2)")

# Initialize an empty dataframe to store the combined results
full_combined_dataframe <- data.frame()

# Loop through each dataframe in modified_dataframes
for(i in seq_along(modified_dataframes)) {
  # Get the current dataframe
  current_dataframe <- modified_dataframes[[i]]
  
  # Perform the left join with A12_all_top
  joined_df <- left_join(current_dataframe, A12_all_top, by = "Sample_ID")
  
  # Combine the joined dataframe with the full_combined_dataframe
  full_combined_dataframe <- bind_rows(full_combined_dataframe, joined_df)
}

# Now, full_combined_dataframe contains all joined dataframes combined into one and filtered out NAs
# all peaks combining A2 top all biogeochem data with orbitrap data
full_combined_dataframe_A12top <-   full_combined_dataframe %>% filter(!is.na(Intensity))

write.csv(full_combined_dataframe_A12top, "full_combined_dataframe_A12top.csv")


#reviewing combined data

ggplot(full_combined_dataframe_A12topPeak, aes(y = relative_intensity, x = H.C)) +
  geom_point() +
  labs(title = "Relative Intensity vs H.C",
       y = "Relative Intensity",
       x = "H.C")

ggplot(full_combined_dataframe_A12topPeak, aes(y = relative_intensity, x = doc_mgL)) +
  geom_point() +
  labs(title = "",
       x = "DOC (mg/L)",
       y = "Relative Intensity")

ggplot(full_combined_dataframe_A12top, aes(x = Sample_ID, y = relative_intensity)) +
  geom_boxplot() +
  #ylim(0,0.1)+
  labs(title = "Distribution of Relative Intensities per Sample",
       x = "Sample ID",
       y = "Relative Intensity")

##investigating relationships
# Select columns of interest
selected_columns <- full_combined_dataframe_A12topPeak[, c("relative_intensity", "pH", "EC", "doc_mgL", "dic_mgL", "Sr", "Suva254", 
                                                           "Fe_mgL", "Pb_mgL", "Zn_mgL", "Mn_mgL", "Ni_mgL","Al_mgL", "As_mgL", "Cu_mgL", "Cd_mgL", "AI.mod_weighted", "DBE_weighted")]
selected_columns <- full_combined_dataframe_A12topPeak[, c("relative_intensity", "C1_Fmax", "C2_Fmax", "C3_Fmax")]

# Create the GGally plot
ggpairs_plot <- ggpairs(selected_columns)

# Print the plot
print(ggpairs_plot)

####for creating second data frame with spearman rank correlation coefficients 

# Access the dataframe directly by its name
df <- modified_dataframes[["A12_TH1_21"]]
df <- df %>% filter(!is.na(Intensity))

# Now you can perform operations on df
# For example, calculate the Spearman rank correlation for O.C and H.C ratios
spearman_test <- cor.test(df$`H.C`, df$`O.C`, method = "spearman", use = "complete.obs")
spearman_cor <- spearman_test$estimate

# Output the Spearman correlation coefficient
print(spearman_cor)

# Add a new column with the Spearman correlation coefficient
df$Spearman_Cor <- rep(spearman_cor, times = nrow(df))

# You can now plot or perform other analyses with df as needed





######A12 bottom################################################################################
#formula attribution for A12BH from icbm
setwd("~/Desktop/orbit trap RZ/Processed_data")
A12_BH_orbitrap <- read.csv("Formula_attribution_A12BH_Jan2024.csv")

names(A12_BH_orbitrap)[names(A12_BH_orbitrap) == "Sample1_A12.BH.32_1.csv"] <- "A12_BH1_32"
names(A12_BH_orbitrap)[names(A12_BH_orbitrap) == "Sample2_A12.BH.34_1.csv"] <- "A12_BH1_34"
names(A12_BH_orbitrap)[names(A12_BH_orbitrap) == "Sample3_A12.BH.35_1.csv"] <- "A12_BH1_35"
names(A12_BH_orbitrap)[names(A12_BH_orbitrap) == "Sample4_A12.BH1.36.csv"] <- "A12_BH1_36"
names(A12_BH_orbitrap)[names(A12_BH_orbitrap) == "Sample5_A12.Bot33_1.csv"] <- "A12_BH1_33"
names(A12_BH_orbitrap)[names(A12_BH_orbitrap) == "Sample6_A12.bottomhose22_1.csv"] <- "A12_BH1_22"
names(A12_BH_orbitrap)[names(A12_BH_orbitrap) == "Sample7_A12.BH1.23_1.csv"] <- "A12_BH1_23"
names(A12_BH_orbitrap)[names(A12_BH_orbitrap) == "Sample8_A12.BH1.25_1.csv"] <- "A12_BH1_25"
names(A12_BH_orbitrap)[names(A12_BH_orbitrap) == "Sample9_A12.BH1.26_1.csv"] <- "A12_BH1_26"
names(A12_BH_orbitrap)[names(A12_BH_orbitrap) == "Sample10_A12.BH1.27_1.csv"] <- "A12_BH1_27"
names(A12_BH_orbitrap)[names(A12_BH_orbitrap) == "Sample11_A12.BH1.28_1.csv"] <- "A12_BH1_28"
names(A12_BH_orbitrap)[names(A12_BH_orbitrap) == "Sample12_A12.BH1.29_1.csv"] <- "A12_BH1_29"
names(A12_BH_orbitrap)[names(A12_BH_orbitrap) == "Sample13_A12.BH1.31_1.csv"] <- "A12_BH1_31"


# Find all columns that start with 'A12'
sample_cols <- grep("^A12", names(A12_BH_orbitrap), value = TRUE)

# Initialize a list to store each modified dataframe
modified_dataframes <- list()

# Loop through each 'Sample..' column
for(col_name in sample_cols) {
  # Rename the 'Sample..' column to 'Intensity' and create 'Sample_ID'
  names(A12_BH_orbitrap)[names(A12_BH_orbitrap) == col_name] <- 'Intensity'
  A12_BH_orbitrap$Sample_ID <- col_name
  
  # Select the relevant columns
  data_subset <- A12_BH_orbitrap[c('Sample_ID','id', 'mz','diff','reference', 'formula', 'H.C', 'O.C','C','H','O','N','S','P','AI', 'AI.mod', 
                                   'DBE','Aromatic','Highly.unsaturated','Unsaturated', 'Saturated', 
                                   'Intensity')]
  
  # Calculate relative intensity
  total_intensity <- sum(data_subset$Intensity, na.rm = TRUE)
  data_subset$relative_intensity <- (data_subset$Intensity / total_intensity)
  
  # Calculate NOSC for each row
  data_subset$NOSC <- 4 - ((4 * data_subset$C + data_subset$H - 3 * data_subset$N - 2 * data_subset$O) / data_subset$C)
  
  # Multiply and create new columns
  cols_to_multiply <- c('mz', 'H.C', 'O.C', 'AI', 'AI.mod', 'DBE', 'Aromatic', 'Highly.unsaturated', 'Unsaturated','NOSC')
  for(col in cols_to_multiply) {
    new_col_name <- paste(col, 'weighted', sep = '_')
    data_subset[[new_col_name]] <- data_subset[[col]] * data_subset$relative_intensity
  }
  
  # Store the modified dataframe in the list
  modified_dataframes[[col_name]] <- data_subset
  
  # Optionally, reset the column name for the next iteration
  names(A12_BH_orbitrap)[names(A12_BH_orbitrap) == 'Intensity'] <- col_name
}

#rename modified dataframes for revised code
modified_dataframes_A12BH <-   modified_dataframes
# Now, modified_dataframes list contains a dataframe for each 'A12..' column

# Iterate over each dataframe in the list
for(i in seq_along(modified_dataframes_A12BH)) {
  # Extract the name of the dataframe
  df_name <- names(modified_dataframes_A12BH)[i]
  
  # Extract the dataframe
  df <- modified_dataframes_A12BH[[i]]
  
  # Create a filename using the dataframe name
  filename <- paste0(df_name, ".csv")
  
  # Save the dataframe as a CSV file
  write.csv(df, filename, row.names = FALSE)
}

#looking at all spectrum

# Loop through each dataframe in the list
for(i in seq_along(modified_dataframes)) {
  # Get the name of the current sample
  sample_name <- names(modified_dataframes)[i]
  
  # Get the current dataframe
  current_dataframe <- modified_dataframes[[sample_name]]
  
  # Create the ggplot; relative intensity of each sample
  p <- ggplot(current_dataframe, aes(x = mz, y = relative_intensity)) +
    geom_line() +
    theme_minimal() +
    labs(title = paste("Mass Spectrum -", sample_name),
         x = "m/z", y = "Relative Intensity")
  
  # Convert to an interactive plotly object
  interactive_plot <- ggplotly(p)
  
  # Print the interactive plot
  print(interactive_plot)
  
  # Optionally, save the interactive plot to an HTML file
  html_file <- paste0("spectrum_", sample_name, ".html")
  htmlwidgets::saveWidget(interactive_plot, html_file)
}


#now looking at A12 bottom
setwd("~/Desktop/R_Anita_clean_Aug2023")

A12_all_bot <- read_excel("A12_A13_allregsampling.xlsx",sheet="A12_all_bottom")

# Initialize an empty dataframe to store the combined results
full_combined_dataframe <- data.frame()

# Loop through each dataframe in modified_dataframes
for(i in seq_along(modified_dataframes)) {
  # Get the current dataframe
  current_dataframe <- modified_dataframes[[i]]
  
  # Perform the left join with A12_all_top
  joined_df <- left_join(current_dataframe, A12_all_bot, by = "Sample_ID")
  
  # Combine the joined dataframe with the full_combined_dataframe
  full_combined_dataframe <- bind_rows(full_combined_dataframe, joined_df)
}

# Now, full_combined_dataframe contains all joined dataframes combined into one and filtered out NAs
#all peaks combining all A12 bot geochem data with orbitrap data
full_combined_dataframe_A12bot <-   full_combined_dataframe %>% filter(!is.na(Intensity))

write.csv(full_combined_dataframe_A12bot, "full_combined_dataframe_A12bot.csv")

#reviewing combined data


ggplot(full_combined_dataframe_A12botPeak, aes(y = relative_intensity, x = H.C)) +
  geom_point() +
  labs(title = "Relative Intensity vs H.C",
       y = "Relative Intensity",
       x = "H.C")
ggplot(full_combined_dataframe_A12botPeak, aes(y = relative_intensity, x = date)) +
  geom_point() +
  labs(title = "",
       x = "Date",
       y = "Relative Intensity")

ggplot(full_combined_dataframe_A12bot, aes(x = Sample_ID, y = relative_intensity)) +
  geom_boxplot() +
  #ylim(0,0.1)+
  labs(title = "Distribution of Relative Intensities per Sample",
       x = "Sample ID",
       y = "Relative Intensity")

##investigating relationships
# Select columns of interest
selected_columns <- full_combined_dataframe_A12botPeak[, c("relative_intensity", "pH", "EC", "doc_mgL", "dic_mgL", "Sr", "Suva254", 
                                                           "Fe_mgL", "Pb_mgL", "Zn_mgL", "Mn_mgL", "Ni_mgL","Al_mgL", "As_mgL", "Cu_mgL", "Cd_mgL",  "AI.mod_weighted", "DBE_weighted")]

# Create the GGally plot
ggpairs_plot <- ggpairs(selected_columns)

# Print the plot
print(ggpairs_plot)

###creating van krevelen plots#####

#just change the data name to look at A12 top or A12 bottom
# Create a new column 'Class' based on the class columns
data <- full_combined_dataframe_A12top%>%
  mutate(Class = case_when(
    Aromatic == 1 ~ "aromatic",
    Highly.unsaturated == 1 ~ "highly unsaturated",
    Unsaturated == 1 ~ "unsaturated",
    Saturated == 1 ~ "saturated",
    TRUE ~ "other"  # For rows that don't fall into any category
  ))

filtered_data <- data %>%
  filter(Class != "other")

van_krevelen_plot <- ggplot(filtered_data, aes(x = O.C, y = H.C, color = Class)) +
  geom_point(alpha = 0.6) +
  theme_minimal() +
  labs(x = "O/C Ratio", y = "H/C Ratio", title= "A12 Top") +
  scale_color_manual(values = c("aromatic" = "red", 
                                "highly unsaturated" = "blue",
                                "unsaturated" = "green",
                                "saturated" = "orange")) +
  theme_minimal()+theme_bw() +theme(panel.grid.major = element_blank(),
                                    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),  
                                    axis.text.x = element_text(size=14), axis.text.y = element_text(size=14),
                                    axis.title.x = element_text(size=14), axis.title.y = element_text(size=14),
                                    legend.text = element_text(size=14), legend.title = element_text(size=14))


# Print the plot
print(van_krevelen_plot)
# Convert to an interactive plotly object
interactive_plot <- ggplotly(van_krevelen_plot)

# Print the interactive plot
print(interactive_plot)

# Optionally, save the interactive plot to an HTML file
html_file <- paste0("Van Krevelen_A12bot", ".html")
htmlwidgets::saveWidget(interactive_plot, html_file)
ggsave("van_krevelen_plot_A12TH.png", van_krevelen_plot, width = 10, height = 8, units = "in")


full_combined_dataframe_A12topPeak$date <- as.Date(full_combined_dataframe_A12topPeak$date)
full_combined_dataframe_A12botPeak$date <- as.Date(full_combined_dataframe_A12botPeak$date)

###combining A12 top and A12 bottom combined datasets####
#looking at peak data
ggplot()+
  geom_point(data=full_combined_dataframe_A12topPeak, aes(x=date, y=relative_intensity, color="A12 Top"), size=2)+
  geom_line(data=full_combined_dataframe_A12topPeak, aes(x=date, y=relative_intensity, color="A12 Top"))+
  
  geom_point(data=full_combined_dataframe_A12botPeak, aes(x=date, y=relative_intensity, color="A12 Bottom"), size=2)+
  geom_line(data=full_combined_dataframe_A12botPeak, aes(x=date, y=relative_intensity, color="A12 Bottom"))+
  
  scale_color_manual(values=c("A12 Top"="blue", "A12 Bottom"="red")) +
  labs(x = "Date", y = "Relative intensity at peak")+
  scale_x_date(date_breaks="1 month", labels= label_date_short(format = c( "%Y", "%b")),expand = c(0.005,0.005))+
  theme_minimal()+theme_bw() +theme(panel.grid.major = element_blank(),
                                    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),  
                                    axis.text.x = element_text(size=14), axis.text.y = element_text(size=14),
                                    axis.title.x = element_text(size=14), axis.title.y = element_text(size=14),
                                    legend.text = element_text(size=14), legend.title = element_blank(), legend.position = c(.9,.9))

ggsave("relative_intensity_peak.png", width = 10, height = 7.5)

full_combined_dataframe_A12top$date <- as.Date(full_combined_dataframe_A12top$date)
full_combined_dataframe_A12bot$date <- as.Date(full_combined_dataframe_A12bot$date)

 #looking at all data
# Assuming 'full_combined_dataframe_A12top' and 'full_combined_dataframe_A12bot' 
# already have a 'date' and 'relative_intensity' columns

# Add a 'group' column to both dataframes
full_combined_dataframe_A12top$group <- 'A12 Top'
full_combined_dataframe_A12bot$group <- 'A12 Bottom'

# Now plot using both dataframes
ggplot() +
  geom_point(data=full_combined_dataframe_A12top, aes(x=date, y=relative_intensity, color=group, shape=group), size=2) +
  geom_point(data=full_combined_dataframe_A12bot, aes(x=date, y=relative_intensity, color=group, shape=group), size=2) +
  scale_color_manual(values=c("A12 Top"="blue", "A12 Bottom"="red")) +
  scale_shape_manual(values=c("A12 Top"=17, "A12 Bottom"=25)) + # Use triangle for both but with different colors
  labs(x = "Date", y = "Relative intensity at peak", color = "Group", shape = "Group") +
  scale_x_date(date_breaks="1 month", date_labels="%b %Y", expand = c(0.005,0.005)) +
  theme_minimal() + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),  
        axis.text.x = element_text(size=14), 
        axis.text.y = element_text(size=14),
        axis.title.x = element_text(size=14), 
        axis.title.y = element_text(size=14),
        legend.text = element_text(size=14), 
        legend.title = element_blank(),
        legend.position = c(.9,.9)) +
  guides(color=guide_legend(), shape=guide_legend()) # Ensure that legends are merged

ggsave("relative_intensity_peak.png", width = 20, height = 7.5)


ggplot()+
  geom_point(data=full_combined_dataframe_A12top, aes(x=date, y=average_relative_intensity, color="A12 Top"), size=2)+
  geom_line(data=full_combined_dataframe_A12top, aes(x=date, y=average_relative_intensity, color="A12 Top"))+
  
  geom_point(data=full_combined_dataframe_A12bot, aes(x=date, y=average_relative_intensity, color="A12 Bottom"), size=2)+
  geom_line(data=full_combined_dataframe_A12bot, aes(x=date, y=average_relative_intensity, color="A12 Bottom"))+
  
  scale_color_manual(values=c("A12 Top"="blue", "A12 Bottom"="red")) +
  labs(x = "Date", y = "Average Relative intensity")+
  scale_x_date(date_breaks="1 month", labels= label_date_short(format = c( "%Y", "%b")),expand = c(0.005,0.005))+
  theme_minimal()+theme_bw() +theme(panel.grid.major = element_blank(),
                                    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),  
                                    axis.text.x = element_text(size=14), axis.text.y = element_text(size=14),
                                    axis.title.x = element_text(size=14), axis.title.y = element_text(size=14),
                                    legend.text = element_text(size=14), legend.title = element_blank(),legend.position = c(.9,.9))

ggsave("Average_relative_intensity.png", width = 10, height = 7.5)


###looking at distributions with boxplots

# Assuming full_combined_dataframe_A12top$Sample_ID is a character vector
full_combined_dataframe_A12top$Sample_ID <- gsub("A12_", "", full_combined_dataframe_A12top$Sample_ID)
full_combined_dataframe_A12top$Sample_ID <- gsub("TH1_", "", full_combined_dataframe_A12top$Sample_ID)

ggplot(full_combined_dataframe_A12top, aes(x = Sample_ID, y = O.C_weighted)) +
  geom_boxplot() +
  #ylim(0,0.1)+
  labs(title = "A12 Top",
       x = "Sample ID",
       y = "O.C weighted")+
  theme_minimal()+theme_bw() +theme(panel.grid.major = element_blank(),
                                    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),  
                                    axis.text.x = element_text(size=14), axis.text.y = element_text(size=14),
                                    axis.title.x = element_text(size=14), axis.title.y = element_text(size=14),
                                    legend.text = element_text(size=14), legend.title = element_blank(),legend.position = c(.9,.9))

ggsave("O.C_weighted_A12top.png", width = 10, height = 7.5)

# Assuming full_combined_dataframe_A12bot$Sample_ID is a character vector
full_combined_dataframe_A12bot$Sample_ID <- gsub("A12_", "", full_combined_dataframe_A12bot$Sample_ID)
full_combined_dataframe_A12bot$Sample_ID <- gsub("BH1_", "", full_combined_dataframe_A12bot$Sample_ID)

ggplot(full_combined_dataframe_A12bot, aes(x = Sample_ID, y = H.C_weighted)) +
  geom_boxplot() +
  #ylim(0,0.1)+
  labs(title = "A12 Bottom",
       x = "Sample ID",
       y = "H.C weighted")+
  theme_minimal()+theme_bw() +theme(panel.grid.major = element_blank(),
                                    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),  
                                    axis.text.x = element_text(size=14), axis.text.y = element_text(size=14),
                                    axis.title.x = element_text(size=14), axis.title.y = element_text(size=14),
                                    legend.text = element_text(size=14), legend.title = element_blank(),legend.position = c(.9,.9))

ggsave("H.C_weighted_A12bot.png", width = 10, height = 7.5)



######saxony mine sampling ################################################################################
#formula attribution for Saxony from icbm
setwd("~/Desktop/orbit trap saxony/processed_data")
Saxonytour_orbitrap <- read.csv("Formula_attribution_Saxonysampling.csv")

names(Saxonytour_orbitrap)[names(Saxonytour_orbitrap) == "Sample1_DS1.csv"] <- "S_DS"
names(Saxonytour_orbitrap)[names(Saxonytour_orbitrap) == "Sample2_FAS1.csv"] <- "S_FAS"
names(Saxonytour_orbitrap)[names(Saxonytour_orbitrap) == "Sample3_FS1.csv"] <- "S_FS"
names(Saxonytour_orbitrap)[names(Saxonytour_orbitrap) == "Sample4_GAS1.csv"] <- "S_GAS"
names(Saxonytour_orbitrap)[names(Saxonytour_orbitrap) == "Sample5_HU1.csv"] <- "S_HU"
names(Saxonytour_orbitrap)[names(Saxonytour_orbitrap) == "Sample6_KWTE1.csv"] <- "S_KWTE"
names(Saxonytour_orbitrap)[names(Saxonytour_orbitrap) == "Sample7_MLS1.csv"] <- "S_MLS"
names(Saxonytour_orbitrap)[names(Saxonytour_orbitrap) == "Sample8_MRS1.csv"] <- "S_MRS"
names(Saxonytour_orbitrap)[names(Saxonytour_orbitrap) == "Sample9_MS1.csv"] <- "S_MS"
names(Saxonytour_orbitrap)[names(Saxonytour_orbitrap) == "Sample10_Mu.hS1.csv"] <- "S_MuhS"
names(Saxonytour_orbitrap)[names(Saxonytour_orbitrap) == "Sample11_NB1.csv"] <- "S_NB"
names(Saxonytour_orbitrap)[names(Saxonytour_orbitrap) == "Sample12_OF1.csv"] <- "S_OF"
names(Saxonytour_orbitrap)[names(Saxonytour_orbitrap) == "Sample13_RS1.csv"] <- "S_RS"
names(Saxonytour_orbitrap)[names(Saxonytour_orbitrap) == "Sample14_S1461.csv"] <- "S_S146"
names(Saxonytour_orbitrap)[names(Saxonytour_orbitrap) == "Sample15_SGE1.csv"] <- "S_SGE"
names(Saxonytour_orbitrap)[names(Saxonytour_orbitrap) == "Sample16_TF1.csv"] <- "S_TF"
names(Saxonytour_orbitrap)[names(Saxonytour_orbitrap) == "Sample17_THDS1.csv"] <- "S_THDS"
names(Saxonytour_orbitrap)[names(Saxonytour_orbitrap) == "Sample18_THGS1.csv"] <- "S_THGS"
names(Saxonytour_orbitrap)[names(Saxonytour_orbitrap) == "Sample19_TKS1.csv"] <- "S_TKS"
names(Saxonytour_orbitrap)[names(Saxonytour_orbitrap) == "Sample20_TMS1.csv"] <- "S_TMS"
names(Saxonytour_orbitrap)[names(Saxonytour_orbitrap) == "Sample21_TREF1.csv"] <- "S_TReF"
names(Saxonytour_orbitrap)[names(Saxonytour_orbitrap) == "Sample22_TRS1.csv"] <- "S_TRS"
names(Saxonytour_orbitrap)[names(Saxonytour_orbitrap) == "Sample23_TSBL1.csv"] <- "S_TSBL"
names(Saxonytour_orbitrap)[names(Saxonytour_orbitrap) == "Sample24_TSBS1.csv"] <- "S_TSBS"
names(Saxonytour_orbitrap)[names(Saxonytour_orbitrap) == "Sample25_TSS1.csv"] <- "S_TSS"
names(Saxonytour_orbitrap)[names(Saxonytour_orbitrap) == "Sample26_VGS1.csv"] <- "S_VGS"
names(Saxonytour_orbitrap)[names(Saxonytour_orbitrap) == "Sample27_Z1.csv"] <- "S_Z"



# Find all columns that start with 'A12'
sample_cols <- grep("^S_", names(Saxonytour_orbitrap), value = TRUE)

# Initialize a list to store each modified dataframe
modified_dataframes <- list()

# Loop through each 'Sample..' column
for(col_name in sample_cols) {
  # Rename the 'Sample..' column to 'Intensity' and create 'Sample_ID'
  names(Saxonytour_orbitrap)[names(Saxonytour_orbitrap) == col_name] <- 'Intensity'
  Saxonytour_orbitrap$Sample_ID <- col_name
  
  # Select the relevant columns
  data_subset <- Saxonytour_orbitrap[c('Sample_ID','id', 'mz','diff','reference', 'formula', 'H.C', 'O.C','C','H','O','N','S','P','AI', 'AI.mod', 
                                   'DBE','Aromatic','Highly.unsaturated','Unsaturated', 'Saturated', 
                                   'Intensity')]
  
  # Calculate relative intensity
  total_intensity <- sum(data_subset$Intensity, na.rm = TRUE)
  data_subset$relative_intensity <- (data_subset$Intensity / total_intensity)
  
  # Calculate NOSC for each row
  data_subset$NOSC <- 4 - ((4 * data_subset$C + data_subset$H - 3 * data_subset$N - 2 * data_subset$O) / data_subset$C)
  
  # Multiply and create new columns
  cols_to_multiply <- c('mz', 'H.C', 'O.C', 'AI', 'AI.mod', 'DBE', 'Aromatic', 'Highly.unsaturated', 'Unsaturated','NOSC')
  for(col in cols_to_multiply) {
    new_col_name <- paste(col, 'weighted', sep = '_')
    data_subset[[new_col_name]] <- data_subset[[col]] * data_subset$relative_intensity
  }
  
  # Store the modified dataframe in the list
  modified_dataframes[[col_name]] <- data_subset
  
  # Optionally, reset the column name for the next iteration
  names(Saxonytour_orbitrap)[names(Saxonytour_orbitrap) == 'Intensity'] <- col_name
}

#rename modified dataframes for revised code
modified_dataframes_Saxonytour <-   modified_dataframes
# Now, modified_dataframes list contains a dataframe for each 'A12..' column

# Iterate over each dataframe in the list
for(i in seq_along(modified_dataframes_Saxonytour)) {
  # Extract the name of the dataframe
  df_name <- names(modified_dataframes_Saxonytour)[i]
  
  # Extract the dataframe
  df <- modified_dataframes_Saxonytour[[i]]
  
  # Create a filename using the dataframe name
  filename <- paste0(df_name, ".csv")
  
  # Save the dataframe as a CSV file
  write.csv(df, filename, row.names = FALSE)
}

#looking at all spectrum

# Loop through each dataframe in the list
for(i in seq_along(modified_dataframes)) {
  # Get the name of the current sample
  sample_name <- names(modified_dataframes)[i]
  
  # Get the current dataframe
  current_dataframe <- modified_dataframes[[sample_name]]
  
  # Create the ggplot; relative intensity of each sample
  p <- ggplot(current_dataframe, aes(x = mz, y = relative_intensity)) +
    geom_line() +
    theme_minimal() +
    labs(title = paste("Mass Spectrum -", sample_name),
         x = "m/z", y = "Relative Intensity")
  
  # Convert to an interactive plotly object
  interactive_plot <- ggplotly(p)
  
  # Print the interactive plot
  print(interactive_plot)
  
  # Optionally, save the interactive plot to an HTML file
  html_file <- paste0("spectrum_", sample_name, ".html")
  htmlwidgets::saveWidget(interactive_plot, html_file)
}


#now looking at A12 bottom
setwd("~/Desktop/Mine project/summer_mine_tour")

Saxony_tour <- read_excel("Saxony_mines_data.xlsx",sheet="data_drainages")

# Initialize an empty dataframe to store the combined results
full_combined_dataframe <- data.frame()

# Loop through each dataframe in modified_dataframes
for(i in seq_along(modified_dataframes)) {
  # Get the current dataframe
  current_dataframe <- modified_dataframes[[i]]
  
  # Perform the left join with A12_all_top
  joined_df <- left_join(current_dataframe, Saxony_tour, by = "Sample_ID")
  
  # Combine the joined dataframe with the full_combined_dataframe
  full_combined_dataframe <- bind_rows(full_combined_dataframe, joined_df)
}

# Now, full_combined_dataframe contains all joined dataframes combined into one and filtered out NAs

full_combined_dataframe_Saxony_tour <-   full_combined_dataframe %>% filter(!is.na(Intensity))

write.csv(full_combined_dataframe_Saxony_tour, "full_combined_dataframe_Saxony_tour.csv")

#reviewing combined data


ggplot(full_combined_dataframe_Saxony_tour, aes(y = relative_intensity, x = H.C)) +
  geom_point() +
  labs(title = "Relative Intensity vs H.C",
       y = "Relative Intensity",
       x = "H.C")
ggplot(full_combined_dataframe_Saxony_tour, aes(y = relative_intensity, x = date)) +
  geom_point() +
  labs(title = "",
       x = "Date",
       y = "Relative Intensity")

ggplot(full_combined_dataframe_Saxony_tour, aes(x = Sample_ID, y = relative_intensity)) +
  geom_boxplot() +
  #ylim(0,0.1)+
  labs(title = "Distribution of Relative Intensities per Sample",
       x = "Sample ID",
       y = "Relative Intensity")

##investigating relationships
# Select columns of interest
selected_columns <- full_combined_dataframe_Saxony_tour[, c("relative_intensity", "pH", "EC", "doc_mgL", "dic_mgL", "Sr", "Suva254", 
                                                           "Fe_mgL", "Pb_mgL", "Zn_mgL", "Mn_mgL", "Ni_mgL","Al_mgL", "As_mgL", "Cu_mgL", "Cd_mgL",  "AI.mod_weighted", "DBE_weighted")]

# Create the GGally plot
ggpairs_plot <- ggpairs(selected_columns)

# Print the plot
print(ggpairs_plot)

###creating van krevelen plots#####

#just change the data name to look at A12 top or A12 bottom
# Create a new column 'Class' based on the class columns
data <- full_combined_dataframe_Saxony_tour%>%
  mutate(Class = case_when(
    Aromatic == 1 ~ "aromatic",
    Highly.unsaturated == 1 ~ "highly unsaturated",
    Unsaturated == 1 ~ "unsaturated",
    Saturated == 1 ~ "saturated",
    TRUE ~ "other"  # For rows that don't fall into any category
  ))

filtered_data <- data %>%
  filter(Class != "other")

van_krevelen_plot <- ggplot(filtered_data, aes(x = O.C, y = H.C, color = Class)) +
  geom_point(alpha = 0.6) +
  theme_minimal() +
  labs(x = "O/C Ratio", y = "H/C Ratio", title= "All saxony") +
  scale_color_manual(values = c("aromatic" = "red", 
                                "highly unsaturated" = "blue",
                                "unsaturated" = "green",
                                "saturated" = "orange")) +
  theme_minimal()+theme_bw() +theme(panel.grid.major = element_blank(),
                                    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),  
                                    axis.text.x = element_text(size=14), axis.text.y = element_text(size=14),
                                    axis.title.x = element_text(size=14), axis.title.y = element_text(size=14),
                                    legend.text = element_text(size=14), legend.title = element_text(size=14))


# Print the plot
print(van_krevelen_plot)
# Convert to an interactive plotly object
interactive_plot <- ggplotly(van_krevelen_plot)

# Print the interactive plot
print(interactive_plot)

# Optionally, save the interactive plot to an HTML file
html_file <- paste0("Van Krevelen_Saxony_tour", ".html")
htmlwidgets::saveWidget(interactive_plot, html_file)
ggsave("van_krevelen_plot_Saxony.png", van_krevelen_plot, width = 10, height = 8, units = "in")


######saxony mine sampling ################################################################################
#formula attribution for Saxony from icbm
setwd("~/Desktop/orbitrap B3/processed_data")
B3_orbitrap <- read.csv("Formula_attribution_B3.csv")

names(B3_orbitrap)[names(B3_orbitrap) == "Sample1_B3.03.csv"] <- "B3_03"
names(B3_orbitrap)[names(B3_orbitrap) == "Sample2_B3.04.csv"] <- "B3_04"
names(B3_orbitrap)[names(B3_orbitrap) == "Sample3_B3.05.csv"] <- "B3_05"
names(B3_orbitrap)[names(B3_orbitrap) == "Sample4_B3.06.csv"] <- "B3_06"
names(B3_orbitrap)[names(B3_orbitrap) == "Sample5_B3.07.csv"] <- "B3_07"
names(B3_orbitrap)[names(B3_orbitrap) == "Sample6_B3.09.csv"] <- "B3_09"
names(B3_orbitrap)[names(B3_orbitrap) == "Sample7_B3.11.csv"] <- "B3_11"
names(B3_orbitrap)[names(B3_orbitrap) == "Sample8_B3.12.csv"] <- "B3_12"
names(B3_orbitrap)[names(B3_orbitrap) == "Sample9_B3.15.csv"] <- "B3_15"
names(B3_orbitrap)[names(B3_orbitrap) == "Sample10_B3.16.csv"] <- "B3_16"
names(B3_orbitrap)[names(B3_orbitrap) == "Sample11_B3.20.csv"] <- "B3_20"
names(B3_orbitrap)[names(B3_orbitrap) == "Sample12_B3.21.csv"] <- "B3_21"
names(B3_orbitrap)[names(B3_orbitrap) == "Sample13_B3.22.csv"] <- "B3_22"
names(B3_orbitrap)[names(B3_orbitrap) == "Sample14_B3.23.csv"] <- "B3_23"
names(B3_orbitrap)[names(B3_orbitrap) == "Sample15_B3.25.csv"] <- "B3_25"
names(B3_orbitrap)[names(B3_orbitrap) == "Sample16_B3.26.csv"] <- "B3_26"
names(B3_orbitrap)[names(B3_orbitrap) == "Sample17_B3.27.csv"] <- "B3_27"
names(B3_orbitrap)[names(B3_orbitrap) == "Sample18_B3.28.csv"] <- "B3_28"
names(B3_orbitrap)[names(B3_orbitrap) == "Sample19_B3.29.csv"] <- "B3_29"
names(B3_orbitrap)[names(B3_orbitrap) == "Sample20_B3.31.csv"] <- "B3_31"
names(B3_orbitrap)[names(B3_orbitrap) == "Sample21_B3.32.csv"] <- "B3_32"
names(B3_orbitrap)[names(B3_orbitrap) == "Sample22_B3.33.csv"] <- "B3_33"
names(B3_orbitrap)[names(B3_orbitrap) == "Sample23_B3.35.csv"] <- "B3_35"
names(B3_orbitrap)[names(B3_orbitrap) == "Sample24_B3.36.csv"] <- "B3_36"
names(B3_orbitrap)[names(B3_orbitrap) == "Sample25_B3.37.csv"] <- "B3_37"


# Find all columns that start with 'A12'
sample_cols <- grep("^B3", names(B3_orbitrap), value = TRUE)

# Initialize a list to store each modified dataframe
modified_dataframes <- list()

# Loop through each 'Sample..' column
for(col_name in sample_cols) {
  # Rename the 'Sample..' column to 'Intensity' and create 'Sample_ID'
  names(B3_orbitrap)[names(B3_orbitrap) == col_name] <- 'Intensity'
  B3_orbitrap$Sample_ID <- col_name
  
  # Select the relevant columns
  data_subset <- B3_orbitrap[c('Sample_ID','id', 'mz','diff','reference', 'formula', 'H.C', 'O.C','C','H','O','N','S','P','AI', 'AI.mod', 
                                       'DBE','Aromatic','Highly.unsaturated','Unsaturated', 'Saturated', 
                                       'Intensity')]
  
  # Calculate relative intensity
  total_intensity <- sum(data_subset$Intensity, na.rm = TRUE)
  data_subset$relative_intensity <- (data_subset$Intensity / total_intensity)
  
  # Calculate NOSC for each row
  data_subset$NOSC <- 4 - ((4 * data_subset$C + data_subset$H - 3 * data_subset$N - 2 * data_subset$O) / data_subset$C)
  
  # Multiply and create new columns
  cols_to_multiply <- c('mz', 'H.C', 'O.C', 'AI', 'AI.mod', 'DBE', 'Aromatic', 'Highly.unsaturated', 'Unsaturated','NOSC')
  for(col in cols_to_multiply) {
    new_col_name <- paste(col, 'weighted', sep = '_')
    data_subset[[new_col_name]] <- data_subset[[col]] * data_subset$relative_intensity
  }
  
  # Store the modified dataframe in the list
  modified_dataframes[[col_name]] <- data_subset
  
  # Optionally, reset the column name for the next iteration
  names(B3_orbitrap)[names(B3_orbitrap) == 'Intensity'] <- col_name
}

#rename modified dataframes for revised code
modified_dataframes_B3 <-   modified_dataframes

# Iterate over each dataframe in the list
for(i in seq_along(modified_dataframes_B3)) {
  # Extract the name of the dataframe
  df_name <- names(modified_dataframes_B3)[i]
  
  # Extract the dataframe
  df <- modified_dataframes_B3[[i]]
  
  # Create a filename using the dataframe name
  filename <- paste0(df_name, ".csv")
  
  # Save the dataframe as a CSV file
  write.csv(df, filename, row.names = FALSE)
}


#looking at all spectrum

# Loop through each dataframe in the list
for(i in seq_along(modified_dataframes)) {
  # Get the name of the current sample
  sample_name <- names(modified_dataframes)[i]
  
  # Get the current dataframe
  current_dataframe <- modified_dataframes[[sample_name]]
  
  # Create the ggplot; relative intensity of each sample
  p <- ggplot(current_dataframe, aes(x = mz, y = relative_intensity)) +
    geom_line() +
    theme_minimal() +
    labs(title = paste("Mass Spectrum -", sample_name),
         x = "m/z", y = "Relative Intensity")
  
  # Convert to an interactive plotly object
  interactive_plot <- ggplotly(p)
  
  # Print the interactive plot
  print(interactive_plot)
  
  # Optionally, save the interactive plot to an HTML file
  html_file <- paste0("spectrum_", sample_name, ".html")
  htmlwidgets::saveWidget(interactive_plot, html_file)
}


#now looking at A12 bottom
setwd("~/Desktop/R_Anita_clean_Aug2023")

B3_data <- read.csv("fulldat_address_Jan172024.csv")
B3_data <- B3_data %>% filter(site_id =="B3")
B3_data <- B3_data[,-1]

# Initialize an empty dataframe to store the combined results
full_combined_dataframe <- data.frame()

# Loop through each dataframe in modified_dataframes
for(i in seq_along(modified_dataframes)) {
  # Get the current dataframe
  current_dataframe <- modified_dataframes[[i]]
  
  # Perform the left join with A12_all_top
  joined_df <- left_join(current_dataframe, B3_data, by = "Sample_ID")
  
  # Combine the joined dataframe with the full_combined_dataframe
  full_combined_dataframe <- bind_rows(full_combined_dataframe, joined_df)
}

# Now, full_combined_dataframe contains all joined dataframes combined into one and filtered out NAs

full_combined_dataframe_B3 <-   full_combined_dataframe %>% filter(!is.na(Intensity))

write.csv(full_combined_dataframe_B3, "full_combined_dataframe_B3.csv")

#reviewing combined data


ggplot(full_combined_dataframe_B3, aes(y = relative_intensity, x = H.C)) +
  geom_point() +
  labs(title = "Relative Intensity vs H.C",
       y = "Relative Intensity",
       x = "H.C")
ggplot(full_combined_dataframe_B3, aes(y = relative_intensity, x = date)) +
  geom_point() +
  labs(title = "",
       x = "Date",
       y = "Relative Intensity")

ggplot(full_combined_dataframe_B3, aes(x = Sample_ID, y = relative_intensity)) +
  geom_boxplot() +
  #ylim(0,0.1)+
  labs(title = "Distribution of Relative Intensities per Sample",
       x = "Sample ID",
       y = "Relative Intensity")

##investigating relationships
# Select columns of interest
selected_columns <- full_combined_dataframe_B3[, c("relative_intensity", "pH", "EC", "doc_mgL", "dic_mgL", "Sr", "Suva254", 
                                                            "Fe_mgL", "Pb_mgL", "Zn_mgL", "Mn_mgL", "Ni_mgL","Al_mgL", "As_mgL", "Cu_mgL", "Cd_mgL",  "AI.mod_weighted", "DBE_weighted")]

# Create the GGally plot
ggpairs_plot <- ggpairs(selected_columns)

# Print the plot
print(ggpairs_plot)

###creating van krevelen plots#####

#just change the data name to look at A12 top or A12 bottom
# Create a new column 'Class' based on the class columns
data <- full_combined_dataframe_B3%>%
  mutate(Class = case_when(
    Aromatic == 1 ~ "aromatic",
    Highly.unsaturated == 1 ~ "highly unsaturated",
    Unsaturated == 1 ~ "unsaturated",
    Saturated == 1 ~ "saturated",
    TRUE ~ "other"  # For rows that don't fall into any category
  ))

filtered_data <- data %>%
  filter(Class != "other")

van_krevelen_plot <- ggplot(filtered_data, aes(x = O.C, y = H.C, color = Class)) +
  geom_point(alpha = 0.6) +
  theme_minimal() +
  labs(x = "O/C Ratio", y = "H/C Ratio", title= "B3") +
  scale_color_manual(values = c("aromatic" = "red", 
                                "highly unsaturated" = "blue",
                                "unsaturated" = "green",
                                "saturated" = "orange")) +
  theme_minimal()+theme_bw() +theme(panel.grid.major = element_blank(),
                                    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),  
                                    axis.text.x = element_text(size=14), axis.text.y = element_text(size=14),
                                    axis.title.x = element_text(size=14), axis.title.y = element_text(size=14),
                                    legend.text = element_text(size=14), legend.title = element_text(size=14))


# Print the plot
print(van_krevelen_plot)
# Convert to an interactive plotly object
interactive_plot <- ggplotly(van_krevelen_plot)

# Print the interactive plot
print(interactive_plot)

# Optionally, save the interactive plot to an HTML file
html_file <- paste0("Van Krevelen_B3", ".html")
htmlwidgets::saveWidget(interactive_plot, html_file)
ggsave("van_krevelen_plot_B3.png", van_krevelen_plot, width = 10, height = 8, units = "in")

