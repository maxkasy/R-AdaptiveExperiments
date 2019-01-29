### Upload csv data in the following format:
**Previous data**: Columns labeled
* *outcome*:  binary outcome variable
* *treatment1, treatment2, ...*, (for treatment dummies)  
or *treatment* (for categorical variable)
* *covar1, covar2, ...*, for covariates (strata will be created from full interactions of controls)  
The app allows for the case that no covariates are provided.

**New wave**: Columns labeled 
* *covar1, covar2, ...* (same as previous data)
* Any other variables, which will be ignored by the app, but might serve as identifiers in treatment dataset.  

**Replicate draws**: A larger choice will yield less random treatment counts, but results in slower computation.

