#Esther Olabisi-A
#Assignment 2


####DATA EXPLORATION ----
##Load libraries
library(rentrez)
library(tidyverse)
library(Biostrings)
library(randomForest)
library(ggplot2)
library(colorspace)

#View databases
entrez_dbs()

#Search nuccore for asteroidea or starfish COI sequences. The max hit is 1,978; less than 3,000 but still okay. Setting use_history as TRUE since there are many sequences.
Asteroidea <- entrez_search(db = "nuccore", term = "Asteroidea[ORGN] OR Starfish[ORGN] AND COI[Gene] AND 100:600[SLEN]", retmax = 3000, use_history = T)

#The large retmax in the entrez_search above does not allow me to view a summary since the summary output is limited to 500 hits. So, I would have to reduce the retmax to 500 to obtain a summary and take a peek at the data before fetching. Here, I am using the webhistory from above.  
Ast_summ <- entrez_summary(db = "nuccore", web_history = Asteroidea$web_history, retmax = 500)

#Check for the organisms to make sure they are what I need. A quick Google search on some of the names yields Asteroida organisms as expected.
unique(extract_from_esummary(Ast_summ, "organism"))

#Remove summary because it is not needed downstream
rm(Ast_summ)

#Obtain large (web_history) fasta records for the sequences using the code in Entrez_Functions.R from courselink.  
source("Entrez_Functions.R") #first load the functions so that I can use them here.

#I used 1000 sequences per file just to avoid ohaving many separate files
Ast_fetch <- FetchFastaFiles(searchTerm = "Asteroidea[ORGN] OR Starfish[ORGN] AND COI[Gene] AND 100:600[SLEN]", seqsPerFile = 1000, fastaFileName = "Asteroidea_COI")

#Using another code from Entrez_Functions.R, merge the fetched files:
Asteroidea_COI <- MergeFastaFiles(filePattern = "Asteroidea_COI")

#Add species column which will contain the species name coined from the 2nd (2L) and 3rd(3L) words of the title. Then, rearrange columns
Asteroidea_COI$Species <- word(Asteroidea_COI$Title, 2L, 3L)
Asteroidea_COI <- Asteroidea_COI[, c("Title", "Species", "Sequence")]


#Search for tunicates. Max hit is 2,140; close to Asteroida's count and good for class balance
Tunicata <- entrez_search(db = "nuccore", term = "Tunicata[ORGN] OR Tunicate[ORGN] AND COI[Gene] AND 100:600[SLEN]" , retmax = 3000, use_history = TRUE)

#Fetch fasta records and merge files into a dataframe
Tun_fetch <- FetchFastaFiles(searchTerm = "Tunicata[ORGN] OR Tunicate[ORGN] AND COI[Gene] AND 100:600[SLEN]", seqsPerFile = 1000, fastaFileName = "Tunicata_COI")
Tunicata_COI <- MergeFastaFiles(filePattern = "Tunicata_COI")

#Add species column, then rearrange columns
Tunicata_COI$Species <- word(Tunicata_COI$Title, 2L, 3L)
Tunicata_COI <- Tunicata_COI[, c("Title", "Species", "Sequence")]

#Next, I would like to merge both both COI dataframe rows into one so that we only work with one dataframe onward. However, I would first add another column indicating the taxonomy for each group so that things don't get mixed up
Asteroidea_COI$Taxon <- "Asteroidea"
Tunicata_COI$Taxon <- "Tunicata"

#Add dataframes together by rows
dfAsteroidea_Tunicata <- rbind(Asteroidea_COI, Tunicata_COI)

#Rearrange the columns so that "taxon" is next to species"
dfAsteroidea_Tunicata <- dfAsteroidea_Tunicata[, c("Title", "Species", "Taxon", "Sequence")]

#Checking how many of each taxa we have just to be sure it worked without any mix up
table(dfAsteroidea_Tunicata$Taxon)



####QUALITY CONTROL ---- 
#checking the sequences for non-nucleotide letters using uniqueLetters(). I could have used alphabetFrequency() but I think uniqueLetters() is more straightforward in that I do not have to look for counts; the letters are either there or not.
#First convert to a DNA String set. 
dfAsteroidea_Tunicata$Sequence <- DNAStringSet(dfAsteroidea_Tunicata$Sequence)

#The results show that both data sets contain various non-nucleotide letters (including Ns).
uniqueLetters(dfAsteroidea_Tunicata$Sequence)

#Therefore, I will filter out sequences containing those alphabets (representing low quality data) while only leaving behind those containing N for further filtering. I am not taking out all sequences with Ns; only those with a certain percentage of N, so that step will occur later.
Filtered_dfCOI <- dfAsteroidea_Tunicata[!grepl('M|R|W|S|Y|K|V|H|D|B', dfAsteroidea_Tunicata$Sequence), ]

#Check that the filtering worked
uniqueLetters(Filtered_dfCOI$Sequence)

#Next, I will further filter the sequences to remove any NAs, bordering non-nucleotide strings, and rows containing Ns that make up over 5% of their sequences. I think 5% is a safe limit since we're working at higher taxonomic levels. I did not exactly see any dashes or bordering Ns upon visual inspection but I think it is better to carry out these steps to be on the safe side

Filtered_dfCOI$Sequence <- as.character(Filtered_dfCOI$Sequence) #Convert the sequences to character first so that I can use tidyverse functions 

Filtered_dfCOI2 <- Filtered_dfCOI %>%
    filter(!is.na(Sequence)) %>%
    mutate(Sequence2 = Sequence) %>%
    mutate(Sequence2 = str_remove(Sequence2, "^[-N]+")) %>%
    mutate(Sequence2 = str_remove(Sequence2, "[-N]+$")) %>%
    mutate(Sequence2 = str_remove_all(Sequence2, "-+")) %>%
    filter(str_count(Sequence2, "N") <= (0.05 * str_count(Sequence2)))
  
#Check for sequence length variations among both taxa. 
#Asteroidea's summary shows wide-ranging sequence lengths indicated by sequences as short at 125 and as long as 600. This variation suggests that constraining the sequence lengths would be better for consistency in classifying the taxa. That is because, as learned in previous scripts, extremely short and extremely long sequences can have different properties, potentially causing issues with identification. So, constraining the lengths within a more balanced range would be better for consistency.
  Ast_summ <- summary(nchar(Filtered_dfCOI2$Sequence2[Filtered_dfCOI2$Taxon == "Asteroidea"]))
  Ast_summ
  
  #The histogram shows the same thing as we see bars demonstrating various sequence lengths
  hist(x = nchar(Filtered_dfCOI2$Sequence2[Filtered_dfCOI2$Taxon == "Asteroidea"]), main = "Histogram of Asteroidea Sequence Lengths", xlab = "Sequence Length (Number of Nucleotides)")  

#Repeat for Tunicata_COI. We see a similar situation here as in Asteroidea's sequence lengths although most of the sequences here are on the longer end. 
  Tun_summ <- summary(nchar(Filtered_dfCOI2$Sequence2[Filtered_dfCOI2$Taxon == "Tunicata"]))
  Tun_summ
  
  hist(x = nchar(Filtered_dfCOI2$Sequence2[Filtered_dfCOI2$Taxon == "Tunicata"]), main = "Histogram of Tunicata Sequence Lengths", xlab = "Sequence Length (Number of Nucleotides)")  
#To constrain sequence lengths using the 1st and 3rd quartiles as boundaries, define 1st and 3rd quartiles
Ast_Q1 <- quantile(nchar(Filtered_dfCOI2$Sequence2[Filtered_dfCOI2$Taxon == "Asteroidea"]), probs = 0.25, na.rm = TRUE)
Ast_Q3 <-quantile(nchar(Filtered_dfCOI2$Sequence2[Filtered_dfCOI2$Taxon == "Asteroidea"]), probs = 0.75, na.rm = TRUE)
 
Tun_Q1 <- quantile(nchar(Filtered_dfCOI2$Sequence2[Filtered_dfCOI2$Taxon == "Tunicata"]), probs = 0.25, na.rm = TRUE)
Tun_Q3 <- quantile(nchar(Filtered_dfCOI2$Sequence2[Filtered_dfCOI2$Taxon == "Tunicata"]), probs = 0.75, na.rm = TRUE)

Ast_Q1
Ast_Q3
Tun_Q1
Tun_Q3

#Filter Asteroida and Tunicata sequences to keep only those within the 1st and 3rd quartile boundaries where the former acts as the minimum and the latter as the maximum length. 
Filtered_dfCOI2$Sequence2 <- as.character(Filtered_dfCOI2$Sequence2)

Asteroidea_Tunicata_COI <- Filtered_dfCOI2 %>%
  filter(str_count(Sequence2) >= Ast_Q1 & str_count(Sequence2) <= Ast_Q3 & Taxon == "Asteroidea" | str_count(Sequence2) >= Tun_Q1 & str_count(Sequence2) <= Tun_Q3 & Taxon == "Tunicata")  
                    
#Check if the above worked and plot new histograms
summary(nchar(Asteroidea_Tunicata_COI$Sequence2[Asteroidea_Tunicata_COI$Taxon == "Asteroidea" ]))
summary(nchar(Asteroidea_Tunicata_COI$Sequence2[Asteroidea_Tunicata_COI$Taxon == "Tunicata" ]))


#Using ggplot for customizable histograms. In the first line, the data to be used includes all rows of taxon Asteroidea. For x in the second line, we input all sequence lengths (nchar) belonging to the taxon Asteroidea. Subsequent lines are for visual changes to axes titles, plot titles, fonts, etc.   
ggplot(data = Asteroidea_Tunicata_COI[Asteroidea_Tunicata_COI$Taxon == "Asteroidea",]) +
  geom_histogram(mapping = aes(x = nchar(Sequence2[Taxon == "Asteroidea"])),
                 fill = "coral3", col = "brown", alpha = 0.8, size = 1, binwidth = 20) +
  labs(title = "Histogram of Asteroidea Sequence Lengths", x = "COI Sequence Length (Number of Nucleotides)", y = "Frequency" ) +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 12),
    axis.title.y = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 12),
    axis.line = element_line(size = 1, colour = "darkred")
  )

#Repeat for Tunicata's histogram.
ggplot(data = Asteroidea_Tunicata_COI[Asteroidea_Tunicata_COI$Taxon == "Tunicata",]) +
  geom_histogram(mapping = aes(x = nchar(Sequence2[Taxon == "Tunicata"])), fill = "coral3", col = "brown", alpha = 0.8, size = 1, binwidth = 10) +
  labs(title = "Histogram of Tunicata Sequence Lengths", x = "COI Sequence Length (Number of Nucleotides)", y = "Frequency" ) +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 12),
    axis.title.y = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 12),
    axis.line = element_line(size = 1, colour = "darkred")
  )
#The new histograms show that the sequences now have lengths more uniformly distributed over a specific range (i.e. within the quartiles)


#Remove datasets that will not be needed downstream 
rm(Ast_fetch, Tun_fetch, Asteroidea, Tunicata, Asteroidea_COI, Tunicata_COI, Filtered_dfCOI)

#Now, we are left with 1,497 sequences for Asteroidea and 1,312 for Tunicata. The latter would be our maximum sample size for each group. Also, as previously learned, the sample sizes do not have to be equal but I am choosing to keep it this way to avoid potential issues with class balance.


####Obtaining sequence features---- 
#We will have two taxonomic classifiers where one uses nucleotide proportion and the other, trinucleotide (3mers) frequencies, as distinguishing features.

#change sequences to DNA String set to use Biostrings functions
Asteroidea_Tunicata_COI$Sequence2 <- DNAStringSet(Asteroidea_Tunicata_COI$Sequence2)

#Compute the frequencies of all four nucleotides  
Asteroidea_Tunicata_COI <- cbind(Asteroidea_Tunicata_COI, as.data.frame(letterFrequency(Asteroidea_Tunicata_COI$Sequence2, letters = c("A", "G", "T", "C"))))

#Computing the proportion of all nucleotides since we want to have frequency in relation to each sequence's length. I am calculating for all four nucleotides so that they can all be visually depicted downstream
Asteroidea_Tunicata_COI$Aprop <- (Asteroidea_Tunicata_COI$A) / (Asteroidea_Tunicata_COI$A + Asteroidea_Tunicata_COI$G + Asteroidea_Tunicata_COI$T + Asteroidea_Tunicata_COI$C)

Asteroidea_Tunicata_COI$Gprop <- (Asteroidea_Tunicata_COI$G) / (Asteroidea_Tunicata_COI$A + Asteroidea_Tunicata_COI$G + Asteroidea_Tunicata_COI$T + Asteroidea_Tunicata_COI$C)

Asteroidea_Tunicata_COI$Tprop <- (Asteroidea_Tunicata_COI$T) / (Asteroidea_Tunicata_COI$A + Asteroidea_Tunicata_COI$G + Asteroidea_Tunicata_COI$T + Asteroidea_Tunicata_COI$C)

Asteroidea_Tunicata_COI$Cprop <- (Asteroidea_Tunicata_COI$C) / (Asteroidea_Tunicata_COI$A + Asteroidea_Tunicata_COI$G + Asteroidea_Tunicata_COI$T + Asteroidea_Tunicata_COI$C)

#Next, we will compute trinucloetide frequencies (k-mer with 3 nucleotides)
Asteroidea_Tunicata_COI$Sequence2 <- DNAStringSet(Asteroidea_Tunicata_COI$Sequence2)
Asteroidea_Tunicata_COI <- cbind(Asteroidea_Tunicata_COI, as.data.frame(trinucleotideFrequency(Asteroidea_Tunicata_COI$Sequence2, as.prob = TRUE)))

#Set aside 393 (~30% of 1,312) randomly selected sequences to validate the trained classifiers later
#convert back to character first
Asteroidea_Tunicata_COI$Sequence2 <- as.character(Asteroidea_Tunicata_COI$Sequence2) 

#Set.seed() so that we will have reproducible results. Group_by() allows the subsequent function run categorically (the taxa are the groups in this case).
set.seed(150)
dfValidation <- Asteroidea_Tunicata_COI %>% 
  group_by(Taxon) %>%
  sample_n(393)

#Filter out sequences already included in the validation dataset so that the classifiers will not include them in training. %in% helps with that. Then, randomly select 919 sequences (70% of the sample size) for randomForest training. 
set.seed(100)
dfTraining <- Asteroidea_Tunicata_COI %>%
  filter(!Title %in% dfValidation$Title) %>%
  group_by(Taxon) %>%
  sample_n(919)
 
#Check that dfValidation and dfTraining contain equal amounts of both taxa. All set!
table(dfValidation$Taxon)
table(dfTraining$Taxon)



####MAIN ANALYSIS: TRAINING THE CLASSIFIERS ---- 
#Taxonomy group is the response variable (Asteroidea or Tunicata). The distingiushing features are nucleotide and trinucleotide (3-mer) proportions. 

#Creating the (SINGLE) NUCLEOTIDE BASED CLASSIFIER: After trial and error with different ntree values, ~500 appears to produce the best predictions. Importance is set to TRUE so that randomForest calculates the level of importance for each distinguishing variable. I will need this information later. Columns 10:13 contain single-nucleotide proportions.
set.seed(100)
Nucleotide_Based_Classifier <- randomForest::randomForest(x = dfTraining[, c(10:13)], y = as.factor(dfTraining$Taxon), ntree = 500, importance = TRUE)

#Creating the KMER BASED CLASSIFIER. ~1000 appears to be ideal for consistent results
set.seed(200)
Kmer_Based_Classifier <- randomForest::randomForest(x = dfTraining[, c(14:77)], y = as.factor(dfTraining$Taxon), ntree = 1000, importance = TRUE)

#Test both!
Nucleotide_Based_Classifier     #All Asteroidea correctly identified. 915/919 Tunicata correctly identified. Nice performance with a 0.22% error rate!

Kmer_Based_Classifier           #918/919 Asteroidea correctly identified and 914/919 Tunicata correctly identified. A 0.33% error rate, also decent.



####Testing the classifiers with unseen data ----
Nucleotide_predict <- predict(Nucleotide_Based_Classifier, dfValidation[, c(10:13)])
Kmer_predict <- predict(Kmer_Based_Classifier, dfValidation[, c(13:77)])

Nucleotide_predict
Kmer_predict

#Create confusion matrices for both predictors 
cf_Nucleotide_predict <- table(observed = dfValidation$Taxon, predicted = Nucleotide_predict)
cf_Kmer_predict <- table(observed = dfValidation$Taxon, predicted = Kmer_predict)

#Consistent performances in both cases!
cf_Nucleotide_predict #100% performance
cf_Kmer_predict #almost perfect

#The importance of each distinguishing feature/variable based on Mean Decrease Accuracy (MDA) and Mean Decrease Gini (MDG) but I will be considering MDA for my project.
importance(Nucleotide_Based_Classifier)
importance(Kmer_Based_Classifier)

#Let's plot MDA for both classifiers!




####Figures depicting variable importance ---- 
#With help from https://stackoverflow.com/a/56306139 and https://stackoverflow.com/a/16968999

#For the NUCLEOTIDE BASED CLASSIFIER: First convert the "importance" data to a dataframe. We need to do this so that we can create custom plots with ggplot. Next, the row names (variables) will be moved to a column so that I can use the content downstream.  
Impt_Nucleotide <- as.data.frame(importance(Nucleotide_Based_Classifier)) %>%
  rownames_to_column() 
Impt_Nucleotide

#Plot MeanDecreaseAccuracy. The reorder() function under aes() sorts the MDA values numerically while the variables (i.e "rowname" column) follow. geom_col() plots MDA to Nucleotide Proportion (the variable). Subsequent lines apply various visual customizations to the plot. Coord_flip() allows for a better visual representation of the decreasing MDA.

#The MDA shows Cprop to be the most importance feature in distinguishing between both taxonomies based on COI gene sequences.
ggplot(data = Impt_Nucleotide) +
  geom_col(mapping = aes(x = reorder(rowname, MeanDecreaseAccuracy), y = MeanDecreaseAccuracy, fill = MeanDecreaseAccuracy)) +
  scale_fill_gradientn(colors = sequential_hcl (4, palette = "Burg", power = 0.1, rev = TRUE, alpha = 1)) +
  labs(title = "Variable Importance: Nucleotide-Based Classifier", fill = "Mean Decrease   Accuracy (MDA)") +
  xlab("Variable (Nucleotide Proportion)") +
  ylab("Mean Decrease Accuracy (MDA)") +
  theme(
    plot.title = element_text(size = 25, face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 15),
    axis.title.y = element_text(face = "bold", size = 15),
    axis.text = element_text(size = 12),
    axis.line = element_line(size = 1.5, colour = "gray"),
  ) +
  coord_flip()


#Next is for the KMER BASED CLASSIFIER. Convert "importance" data to dataframe
Impt_Kmer <- as.data.frame(importance(Kmer_Based_Classifier)) %>%
  rownames_to_column()

#Plot MeanDecreaseAccuracy. 
ggplot(data = Impt_Kmer) + 
  geom_col(mapping = aes(x = reorder(rowname, MeanDecreaseAccuracy), y = MeanDecreaseAccuracy, fill = MeanDecreaseAccuracy)) +
  scale_fill_gradientn(colors = sequential_hcl (n = length(Impt_Kmer$MeanDecreaseAccuracy), palette = "Emrld", rev = TRUE, alpha = 0.8)) +
  labs(title = "Variable Importance: Kmer-Based Classifier", fill = "Mean Decrease Accuracy (MDA)") +
  xlab("Variable (3-mer Proportion)") +
  ylab("Mean Decrease Accuracy (MDA)") +
  theme(
    plot.title = element_text(size = 25, face = "bold", hjust = 0.5), 
    axis.title.x = element_text(face = "bold", size = 15),
    axis.title.y = element_text(face = "bold", size = 15),
    axis.text = element_text(size = 8),
    axis.line = element_line(size = 1.5, colour = "gray"),
  ) +
  coord_flip()


#All done!

