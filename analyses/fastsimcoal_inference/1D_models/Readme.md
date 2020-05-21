For the single population models, all populations use the same template files. Only sample size varies between populations, and sample sizes are pulled from the file projectionValues.txt.

1D.1Epoch refers to a single epoch model, where a constant net population size is inferred

1D.2Epoch refers to a two epoch model, where a current size is inferred, and then allowed to change size at some point T generations in the past. The time of the size change is also inferred.

1D.3Epoch is a three epoch model, with a current population size which changes at some time T generations in the past. The population is allowed to stay at this size for 20 generations, and T+20 generations in the past the size is allowed to change again.

 
