# rgmirt
Project on fitting regularised mirt models to sparse data. The mirt directory contains an adaptation of the R-package [mirt](https://github.com/philchalmers/mirt), that allows for l1 and l2 regularisation of the likelihood. This is useful for sparse data applications. Most code is simply copied from the original repo, with the additional option of adding a penalty to the likelihood. 
The project directory contains code for running simulation studies. It is written to run on a computing cluster. 
