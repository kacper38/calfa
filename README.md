# calfa
reconstruction of ss structure give alpha carbon coordinates

Prediction of secondary structure given only alpha carbons coordinates, 
Coordinates saved in PCA format
Secondary structure verification done using DSSP

Prediction based on neural network


currently trying different types of file formats to feed into the network
:to do:

functional mask for padding, as can't affors noise 

try predictions for segments of vector == human like reading of protein structures
 
reading vector as whole seems better for machine, as sheet may end at 3/4 of subset lenght, which may conclude in local mismatches, and add up to overfitting certain parts, more, as aminoacids are not equally frequent, -> then the structure of coils and sheets, are not equally distributed


maybe add a probability of full sequence instead of sequuence of probs for each part of seq
-> ...
