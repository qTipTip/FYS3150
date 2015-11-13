#include <iostream>
// Computes the Legendre polynomial of degree N
double Legendre( int n, double x) {
  double r, s, t;
  int m;
  
  r = 0; 
  s = 1;
  // Using recurrence relation 
  for (m = 0; m < n; ++m) {
    t = r;
    r = s;
    s = (2 * m + 1)*x*r - m*t;
    s = s / (m+1);
  }
  return s;
}

int main(int argc, char *argv[])
{
  int n = atoi(argv[1]);
  double x = atof(argv[2]);

  std::cout << "The Legendre polynomial of order " << n << " at "
            << x << " is " << Legendre(n, x) << std::endl;  
  return 0;
}
