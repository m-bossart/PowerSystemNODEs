#PSEUDO CODE OUTLINING THE VARIOUS STEPS FOR TRAINING UODE.

#Load a library of fault data in the form of Fourier coefficients
    #Divide into training and test sets.

#Build Full System + IB system in PSID. This is the source of ground truth data.
#Randomly assign parameter values over a certain range.

#Case A: Largest

    #Evaluate loss over the full test set.



#Write a surrogate model with GFM structure plus two additional states for the NODE current source at output.

#Two stage training process.
function initialize_surrogate()
