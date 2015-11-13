double cartesian_function(double x1, double y1, double z1, double x2, double y2, double z2);
double spherical_function(double r1, double r2, double theta1, double theta2, double phi1, double phi2);

void legendre(double x1, double x2, double x[], double w[], int n);
void laguerre(double x[], double w[], int n, double alph);
void bruteforce_MC(double lower_lim, double upper_lim, int N);
void importance_sampl_MC(int N);

double gammln(double xx);

double legendre_integral(double *x, double *w, int N, double (*func)(double, double, double, double, double, double));
double laguerre_integral(double *xr, double *wr, double *x_theta, double *w_theta, double *x_phi, double *w_phi, int N, double (*func)(double, double, double, double, double, double));

