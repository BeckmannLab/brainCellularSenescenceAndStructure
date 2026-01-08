import os
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('outputs.csv')

merged_values = pd.concat([df['left_probe'].fillna(''), df['right_probe'].fillna('')], ignore_index=True)
    clean_values = merged_values[merged_values != '']
    value_counts = clean_values.value_counts().sort_values(ascending=True)  

# Sort for horizontal bars
plt.barh(value_counts.index, value_counts.values, color='skyblue', edgecolor='black', height=0.7)  # Adjust height for spacing
plt.xlabel("Frequency")
plt.margins(y=0.2) 
plt.subplots_adjust(left=0.3, right=0.7) 
plt.savefig(os.path.join(results_dir, "labels_hist.jpg"), format="jpg", bbox_inches="tight")
plt.close()  


library(data.table)
library(ggplot2)
df = fread('outputs.csv')


# Merge values from 'left_probe' and 'right_probe', replacing NA with empty string
merged_values <- c(ifelse(is.na(df$left_probe), "", df$left_probe), 
                   ifelse(is.na(df$right_probe), "", df$right_probe))

# Remove empty values
clean_values <- merged_values[merged_values != ""]

# Count occurrences and sort in ascending order
value_counts <- sort(table(clean_values), decreasing = FALSE)
value_counts2 = as.data.frame(value_counts)

value_counts2$updated_names = value_counts2$clean_values
value_counts2$updated_names = gsub("ctx_rh_G_front_sup", "Right Hemisphere\nSuperior Frontal Gyrus", value_counts2$updated_names)
value_counts2$updated_names = gsub("ctx_lh_G_front_sup", "Left Hemisphere\nSuperior Frontal Gyrus", value_counts2$updated_names)
value_counts2$updated_names = gsub("ctx_rh_G_front_middle", "Right Hemisphere\nMiddle Frontal Gyrus", value_counts2$updated_names)
value_counts2$updated_names = gsub("ctx_lh_G_front_middle", "Left Hemisphere\nMiddle Frontal Gyrus", value_counts2$updated_names)


g <- ggplot(value_counts2, aes(x = Freq, y = reorder(updated_names, Freq))) +
  geom_bar(stat = "identity", fill = "grey40", color = "black", width = 0.4) +  # Reduce width for thinner bars
  labs(x = "Number of Samples", y = "Biopsy Location") +
  scale_x_continuous(
    breaks = seq(0, max(value_counts2$Freq + 10), by = 20),  # Adjust tick frequency
    expand = expansion(mult = c(0.02, 0.1))  # Add more space on the right
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 15),  # Increase x-axis text size
    axis.text.y = element_text(size = 15),  # Increase y-axis text size
    axis.title = element_text(size = 20) 
  )
ggsave("ct_biopsy_location.pdf", g, width = 8, height = 5)
