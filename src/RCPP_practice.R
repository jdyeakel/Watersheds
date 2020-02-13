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

#How to backwards iterate
cppFunction('
NumericVector backit(NumericVector x) {
int l = x.size();
NumericVector y(l);
int tic = 0;
for (int i=l; i --> 0;) {
  y(tic) = x(i);
  tic = tic + 1;
}
NumericVector output = y;
return output;
}
            ')


#Grab integer component
cppFunction('
double intpart(double x) {
double y = floor(x);
return y;
}
            ')

#Apply Changes to a List

cppFunction('
List changelist(List x) {
int l = x.size();
for (int i=0;i<l;i++) {
NumericMatrix m = x[i];
m(1,1) = m(1,1) + 1;
x[i] = m;
}
return x;
}            
            ')


cppFunction('
IntegerVector matrixsize(NumericMatrix M) {
IntegerVector x(2);
x(0) = M.nrow();
x(1) = M.ncol();
return x;
}            
            ')

cppFunction('
int randomint(int num_res) {
NumericVector rdraw_v = runif(1);
double rdraw = as<double>(rdraw_v);
int draw = (int) floor(rdraw*(num_res));
return draw;
}            
            ')

cppFunction('
int rf(double x) {
int y = round(x);
            return y;
            }            
            ')

cppFunction('
double RandomFloat(float min, float max)
{
    float r = (float)rand() / (float)RAND_MAX;
    return min + r * (max - min);
}
            ')

cppFunction('
int mod(int x, int y)
{
    int z = x%y;
    return z;
}
            ')



cppFunction('
int nn(int x, int L) {
      x = x+1;
      int new_x;
      int check = 0;
      int mod = (x)%((L+2));
      //Bottom Row
      if (x <= (L+2)) {new_x = x + L*(L+2); check = 1;}
      //Top Row
      if (x >= (pow((L+2),2) - (L+1))) {new_x = x - L*(L+2); check = 1;}
      //Left Row
      if (mod == 1) {new_x = x + L; check = 1;}
      //Right Row
      if (mod == 0) {new_x = x - L; check = 1;}
      //Bottom Left
      if ((x <= (L+2)) && (mod == 1)) {new_x = x + L*(L+2) + L; check = 1;}
      //Bottom Right
      if ((x <= (L+2)) && (mod == 0)) {new_x = x + L*(L+2) - L; check = 1;}
      //Top Left
      if ((x >= (pow((L+2),2) - (L+1))) && (mod == 1)) {new_x = x - L*(L+2) + L; check = 1;}
      //Top Right
      if ((x >= (pow((L+2),2) - (L+1))) && (mod == 0)) {new_x = x - L*(L+2) - L; check = 1;}
      //Middle
      if(check == 0) {new_x = x;}
      return new_x - 1;
}            
            
            ')

cppFunction('
double mmin(NumericVector x) {
NumericVector::iterator it = std::min_element(x.begin(), x.end());
double y = *it;
return y;
} 
            ')

cppFunction('
IntegerVector erase(IntegerVector vec, int i) {
vec.erase(vec.begin() + i);
return vec;
}         
            ')


cppFunction('
List listed(List x) {
List y(2);
List z(1);
y(0) = x;
y(1) = x;
z(0) = y;
return z;
}
            ')

cppFunction('
IntegerMatrix mm(int x) {
IntegerMatrix m(x,x);
m(x,x) = 1;
return m;
}
            ')

cppFunction('
int findmax(IntegerVector x) {
int y = which_max(x);
return y;
}
            ')

cppFunction('
int dr(IntegerVector nn) {
int num_draw = nn.size();
NumericVector rdraw_v = runif(1);
double rdraw = as<double>(rdraw_v);
int draw = (int) floor(rdraw*(num_draw));
int new_loc = nn(draw);
return new_loc;
}
            ')

cppFunction('
NumericVector drawexp(double num, double x) {
double lam = 1/x;
NumericVector draw = round(rexp(num,lam),0);
return draw;
}
            ')

cppFunction('
NumericVector addend(NumericVector v, double x) {
v.push_back(x);
return v;
}
            ')

cppFunction('
int randint(double min, double max) {
int id = min + (rand() % (int)(max - min + 1.L));
return id;
}
            ')

cppFunction('
double randnum() {
double r = ((double) rand() / (RAND_MAX));
return r;
}
            ')


cppFunction('
List growlist(NumericVector tostart, NumericVector toadd, int n) {
List x(1);
x(0) = tostart;
for (int i=0;i<n;i++) {
x.push_back(toadd);
}
return x;
}
            ')
