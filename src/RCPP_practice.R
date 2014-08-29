rm(list=c(ls()))

library(Rcpp)
library(inline)
src <- 'int n = as<int>(ns); double x = as<double>(xs);
  for (int i=0; i<n; i++) x=1/(1+x);
  return wrap(x);'
l <- cxxfunction(signature(ns="integer",
                           xs="numeric"),
                 body=src, plugin="Rcpp")


#This next code does not use inline...
cppFunction('
  int add(int x, int y, int z) {
    int sum = x + y + z;
    return sum;
  }'
)
add(1,2,3)

cppFunction('
  NumericVector zip(NumericVector x) {
    int n = x.size();
    for (int i=0; i < n; ++i) {
      x[i] = x[i] * x[i];
    }
    return(x);
  }'
)
v <- c(1,2,3)
zip(v)

cppFunction('
int one() {
int x(10);
return x;
}
            ')


cppFunction('
LogicalMatrix one() {
LogicalMatrix x(5,5);
int n = x.nrow();
for (int i=0;i<n;i++) {
  x(i,0) = TRUE;
}
return x;
}
')

cppFunction('
NumericMatrix rvector() {
NumericMatrix X(5,5);
for (int i=0;i<5;i++) {
NumericVector x = runif(5);
for (int j=0;j<5;j++) {
X(i,j) = x(j);
}
}
return X;
}
')

cppFunction('
  NumericVector rowSumsC(NumericMatrix x) {
    int nrow = x.nrow(), ncol = x.ncol();
    NumericVector out(nrow);

    for (int i = 0; i < nrow; i++) {
      double total = 0;
      for (int j = 0; j < ncol; j++) {
        total += x(i, j);
      }
      out[i] = total;
    }
    return out;
  }
')
x <- matrix(sample(100),10)
rowSumsC(x)


#Read in and extract elements of a list
#Use this approach to read in nearest neighbors for all nodes!
cppFunction('
  int readlist(List x_1) {
    int n = x_1.size();
    NumericVector out(n);
    for (int i = 0; i < n; i++) {
      NumericVector sublist = as<NumericVector>(x_1[i]);
      out[i] = sum(sublist);
    }
    return out(1:2);
  }            
')
x <- list(c(1,1,1),c(2,2,2),c(3,3,3,3))
readlist(x)


cppFunction('
int whileloop() {
  int count = 0;
  while (count < 100) {
count = count + 1;
  }
return count;
}
            ')

cppFunction('
int findmax(IntegerVector X) {
int m = which_max(X);
int out = X(m);
return out;
}
            ')

cppFunction('
            int findmin(IntegerVector X) {
            int m = which_min(X);
            int out = m;
            return out;
            }
            ')

cppFunction('
            int findmin(IntegerVector X) {
            int m = min(X);
            int out = m;
            return out;
            }
            ')


cppFunction('
bool testvector(IntegerVector X) {
if ((X(0) == 1) && (X(1) == 1)) {
return true;} 
else {return false;}
}
')

cppFunction('
IntegerVector randomnum(int x, int min_num, int max_num) {
IntegerVector r(x);
for (int i=0;i<x;i++){
r(i) = rand() % (max_num + 1 - min_num) + min_num;
}
return r;
}              
            ')

cppFunction('
int randx() {
return rand();
}
')

cppFunction('
NumericVector mult(IntegerVector X) {
NumericVector new_X(X.size()) = X % 2;
return new_X;
}            
')



cppFunction('
int func(int x) {
int new_x;
new_x = x + 5;
return new_x;
}
            ')

cppFunction('
NumericVector func_norm(double M, double SD) {
NumericVector x = rnorm(20,M,SD);
return x;
}
            ')

cppFunction('
double randx() {
NumericVector x = runif(1);
double X = as<double>(x);
return X;
}
')

cppFunction('
int sizelist(List x) {
NumericVector n= as<NumericVector>(x[0]);
int output = n.size();
return output;
}
            ')




