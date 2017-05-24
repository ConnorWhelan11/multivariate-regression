# multivariate-regression

### Building the code
To build the code do as follows
```sh
git clone https://github.com/unh-hpc/project-1-0800fc577294c34e0b28ad2839435945.git
cd project-1-0800fc577294c34e0b28ad2839435945
autoreconf -iv
./configure
```
### Running the program
```
[cpf37@fishercat mpi]$ ./main
y = (0.080604)*x0 + (0.042525)*x1 + (0.018824)*x2 + error

```

### The Data
The data used in this problem came from the UCI Machine Learning repository at https://archive.ics.uci.edu/ml/datasets/Auto+MPG. The following information has been copied from their website:

Source:

This dataset was taken from the StatLib library which is maintained at Carnegie Mellon University. The dataset was used in the 1983 American Statistical Association Exposition.


Data Set Information:

This dataset is a slightly modified version of the dataset provided in the StatLib library. In line with the use by Ross Quinlan (1993) in predicting the attribute "mpg", 8 of the original instances were removed because they had unknown values for the "mpg" attribute. The original dataset is available in the file "auto-mpg.data-original". 

"The data concerns city-cycle fuel consumption in miles per gallon, to be predicted in terms of 3 multivalued discrete and 5 continuous attributes." (Quinlan, 1993)


Attribute Information:

1. mpg: continuous 
2. cylinders: multi-valued discrete 
3. displacement: continuous 
4. horsepower: continuous 
5. weight: continuous 
6. acceleration: continuous 
7. model year: multi-valued discrete 
8. origin: multi-valued discrete 
9. car name: string (unique for each instance)


### Multivariate Regressin to predict Auto MPG
The goal is to solve the Ax=b problem for the vector of coefficients x. My implentation does so through QR factorization of matrix A into the multiplication of an orthonormal matrix Q and upper triangular matrix R.
/begin{equation*}
x^2
/end{equation*}
