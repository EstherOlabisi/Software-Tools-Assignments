##Esther Olabisi-Adeniyi 
##Student ID: #1238826
##Software Tools (BINF6210) – Assignment #1 Code ----
##10/7/22


##Run tidyverse
library(tidyverse)

##I will be studying the bird taxon Apodidae (common name: swifts) for this assignment. The first step is obtaining data from BOLD using the BOLD API. I have obtained the base link to which I added the taxon i.e. "taxon=Apodidae", and tsv format "format=tsv."
Swifts <- read_tsv("http://www.boldsystems.org/index.php/API_Public/specimen?taxon=Apodidae&format=tsv")

##Set working directory through Session > Set Working Directory > Choose Directory. Check:
getwd()

##Write the data to hard disk and pass to new object:
write_tsv(Swifts, "Swift_BOLD_data.tsv")

Swifts <- read_tsv("Swift_BOLD_data.tsv")

Swifts


##My goal for this project is to see the abundance of Apodidae across different climates based on latitude. Therefore, I will be extracting the latitude column from the Apodidae tibble to create a histogram, allowing me to analyze number of samples across latitudes.

##Using tidyverse:
Lat_Swift <- Swifts %>%
  filter(!is.na(lat)) %>%
  pull(lat)

##First line involves passing the Swift data to the object name that we want for our latitude data eventually. Next, I filtered out NA values from the latitude column. The third line isolates the latitude column into the object Lat_Swift. 

#Checking the new vector:
Lat_Swift
class(Lat_Swift)

#Making sure there are no NA values; '0' is the output, good.
sum(is.na(Lat_Swift))

#Now, I will create a histogram for the latitude values in the vector Swifts_lat
hist(Lat_Swift, breaks = 5, col = "pink", main = paste("Histogram of Latitude for Apodidae (n = 71)"), xlab = "Latitude (degrees)", ylab = "Number of Apodidae Members")

##In the hist() function, I chose 5 breaks based on the length of latitudes. I also passed titles and axes labels as arguments.

##The histogram came out left skewed suggesting that there may be outliers. Therefore, I decided to do a box plot to check for outliers:
ggplot(data = as.data.frame(Lat_Swift)) +
  geom_boxplot(mapping = aes(x = Lat_Swift), outlier.color = "red", outlier.fill = "yellow", outlier.size = 2) +
  labs(title = "Box plot of Latitude (n = 71)", x = "Latitude (degrees)", y = "Number of Apodidae Members") +
  coord_flip()

##In the first line, I used the as.data.frame() function to change Lat_Swift from vector to data frame since ggplot() uses the latter. In the second line, I passed the arguments for any outliers in the box plot. In the third line, I created my title and axes labels. 

##The box plot happens to show no outliers; therefore, I will work with the histogram as is. I will run summary() and standard deviation for Lat_Swift to obtain values that may help with analyses. 
summary(Lat_Swift)
sd(Lat_Swift)

