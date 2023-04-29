#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define nx 4
#define ny 4

/*************************************************************/

double d(double x, double y);

double f(double x, double y);

void Boundary(double x[], double y[], int n, double **a, double *rhs);

double *heat_equation(double x[], double y[], double d(double x, double y), double f(double x, double y));

void interior(double x[], double y[], double d(double x, double y), double f(double x, double y),  double **a, double *rhs);

double choleskydc(double **A, int n);

void choleskylin(double **A, int n, double **L);

double CGm(double **M, double *vetor, double *d, int m, int n, int p);

/************************************************************/
double d(double x, double y){
  double value = 1.0;

  return value;
}

double f(double x, double y){
  
  double value=1.0;
  
  return value;
}

void Boundary(double x[], double y[], int n, double **a, double *rhs){

  int i, j, kc;

  /* Left boundary. */
  j = 0;
  for ( i = 1; i < ny - 1; i++ )
  {
    kc = i * nx + j;
    a[kc][kc] += 1.0;
    rhs[kc] = 10.0;
  }

  /* Right boundary. */
  j = nx - 1;
  for ( i = 1; i < ny - 1; i++ )
  {
    kc = i * nx + j;
    a[kc][kc] += 1.0;
    rhs[kc] = 100.0;
  }

  /* Lower boundary. */
  i = 0;
  for ( j = 0; j < nx; j++ )
  {
    kc = i * nx + j;
    a[kc][kc] += 1.0;
    rhs[kc] = 0.0;
  }

  /* Upper boundary. */
  i = ny - 1;
  for ( j = 0; j < nx; j++ )
  {
    kc = i * nx + j;
    a[kc][kc] += 1.0;
    rhs[kc] = 0.0;
  }
}

void interior(double x[], double y[], double d(double x, double y), double f(double x, double y), double **a, double *rhs) {
    /* stes up the matrix and the right hand side part at all interior nodes */
    
    double dce, dcn, dcs, dcw, dx, dy;
    int ic, in, is, jc, je, jw, kc, ke, kn, ks, kw;

    dx = x[1] - x[0];
    dy = y[1] - y[0];

    for(ic = 1; ic < ny - 1; ic++) {
        for(jc = 1; jc < nx - 1; jc++) {
            in = ic + 1;
            is = ic - 1;
            je = jc + 1;
            jw = jc - 1;

            kc = ic * nx + jc;
            ke = kc + 1;
            kw = kc - 1;
            kn = kc + nx;
            ks = kc - nx;

            dce = d(0.5 * (x[jc] + x[je]), y[ic]);
            dcw = d(0.5 * (x[jc] + x[jw]), y[ic]);
            dcn = d(x[jc], 0.5 * (y[ic] + y[in]));
            dcs = d(x[jc], 0.5 * (y[ic] + y[is]));

            a[kc][kc] = (dce + dcw) / dx / dx + (dcn + dcs) / dy / dy;
            a[kc][ke] = -dce / dx / dx;
            a[kc][kw] = -dcw / dx / dx;
            a[kc][kn] = -dcn / dy / dy;
            a[kc][ks] = -dcs / dy / dy;

            rhs[kc] = f(x[jc], y[ic]);
        }
    }
}

double *subregion(int n, double a, double b){
  /*the interval is divided into N-1 even subintervals,
    and the endpoints of intervals are used as the points.*/
  
  int i;
  double *x;

  x = ( double * ) malloc ( n * sizeof ( double ) );

  if ( n == 1 )
  {
    x[0] = ( a + b ) / 2.0;
  }
  else
  {
    for ( i = 0; i < n; i++ )
    {
      x[i] = ( (double)(n-1-i)*a+(double)(i)*b )/(double)(n-1);
    }
  }
  return x;
}

void mesh2d(double xvec[], double yvec[], double xmat[], double ymat[]){

  /*creates a  2D-mesh from X and Y vectors*/
  int i;
  int j;

  for ( j = 0; j < ny; j++ )
  {
    for ( i = 0; i < nx; i++ )
    {
      xmat[i+j*nx] = xvec[i];
    }
  }

  for ( j = 0; j < ny; j++ )
  {
    for ( i = 0; i < nx; i++ )
    {
      ymat[i+j*nx] = yvec[j];
    }
  }

 return;
}

void meshvectorprint(int n, double **a, char *title){
   int i,j;

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "%s\n", title );
  fprintf ( stdout, "\n" );
  for ( i = 0; i < n; i++ ){
    for (j=0; j<n;j++)
      fprintf ( stdout, "  %8d\t%8d : %14f\n", i,j, a[i][j] );
    }

  return;
}

double *heat_equation(double x[], double y[],double d(double x, double y), double f(double x, double y)){

  /*The generation's grid process:*/
  double **a, *u;
  int n,i,j;

  n = nx * ny; 

  /* the matrix and the right hand side:*/
  a = (double **)malloc(n * sizeof(double *));
  for(i = 0; i < n; i++) {
    a[i] = (double *)malloc(n * sizeof(double));
  }
  
  u = (double *)malloc(n * sizeof(double));
  
  /*matrix interior points*/
   
  interior(x, y, d, f, a, u);
  
  /*handle with the boundary conditions:*/
  Boundary(x, y, n, a, u);
 
  /*Solving the linear system*/
  
  choleskydc(a,n);
  choleskylin(a, n, &u); 
 
  for (int i = 0; i < n; i++) {
    free(a[i]);  // Free each row of the matrix
  }
  free(a);
  
  return u;
}

/************************************************************/

double **ReadMatrix(char *nome, int *m, int *n) {
  double **M;
  int i, j;
  FILE *f;

  f = fopen(nome, "r");

  fscanf(f, "%d %d", m, n);

  M = (double **)malloc(*m * sizeof(double *));

  for (i = 0; i < *m; i++)
    M[i] = (double *)malloc(*n * sizeof(double));

  for (i = 0; i < *m; i++)
    for (j = 0; j < *n; j++)
      fscanf(f, "%lf", &M[i][j]);

  return M;
}

void PrintMatrix(double **M, int m, int n) {
  int i, j;

  for (i = 0; i < m; i++) {
    printf("\n\t");
    for (j = 0; j < n; j++) {
      printf("%.2lf\t", M[i][j]);
    }
  }
}

/***************************************************************/

double VectorNorm(double *Vetor, int m, int p) {
  double max=0;
  int i;
  
  if(p == 0){
    max=fabs(Vetor[0]);
    for(i=0;i<m;i++){
      if(fabs(Vetor[i])>max){
        max=fabs(Vetor[i]);
      }
    }
  }
    else{
      for(i=0;i<m;i++){
        max+=pow(fabs(Vetor[i]),p);
      }
    max = pow(max,1./p);
    }
  return max;
}

double DotProd(double *vetor, double *vetor2, int m) {
  double prod=0;
  int i;
  
  for (i=0; i<m; i++) {
    prod += vetor[i]*vetor2[i];
  }
  return prod;
}

double *MatrixVector(double **Matriz, double *vetor, int m, int n, int mc) {
  int i,k;
  double *R;
  
  if (m != mc) {
    puts("It is not possible to continue the product.");
    return NULL;
  }
  
  R = malloc(m*sizeof(double));
  
    for(i=0; i<m; i++){
      
      R[i] = 0;
      
      for (k=0; k<mc; k++) {
        R[i] += (Matriz[i][k]*vetor[k]);
      }
    }
  return R;
}

double *residue(double **M, double *vetor, int m, int n){
  int i;
  double *Ax, *r;

  Ax = malloc(m*sizeof(double)); 
  r = malloc(m*sizeof(double));

  Ax = MatrixVector(M,vetor,m,n,m);

  for(i=0; i<m; i++){
    r[i] = M[i][n-1] - Ax[i];
  }
  
  return r;
}

/******************************************************************/

double CGm(double **M, double *vetor, double *d, int m, int n, int p){
  int i, it = 0;
  double *r, *Ad, Add, alpha, beta, *d0, rr;
  
  r = d0 = residue(M, vetor, m, n);
  rr = DotProd(r, r, m);
  Ad = MatrixVector(M, d, m, n, m); 
  Add = DotProd(Ad, d, m);
  alpha = rr/Add;
    
  for (i=0; i<m; i++) {
    vetor[i] += alpha*d[i]; 
    r[i] -= alpha*Ad[i];
  }  

  beta = DotProd(r,r, m)/rr;

  // d^{i+1} 
  for(i=0; i<m; i++){
    d[i] = r[i]+ (beta*d[i]);
  }
  
  return VectorNorm(r, m, p);
}

void showMatrix(double *A, int n) {
  for (int i = 0; i < n; i++) {
    printf("\n\t");
    for (int j = 0; j < n; j++) {
      printf("%.2lf\t", A[i * (n + j)]);
    }
    printf("\n");
  }
}

double choleskydc(double **A, int n) { // A=L*(L^T)
  int i, j, k;
  double sum;

  for (int i = 0; i < n; i++) {
    for (int j = i; j < n; j++) {
      for (sum = A[i][j], k = i - 1; k >= 1; k--)
        sum -= A[i][k] * A[j][k];
      if (i == j) {
        if (sum <= 0.0){
          return -1; // indicates failure if the matrix is not positive definite.
        }
        A[i][i] = sqrt(sum);
      } else
        A[j][i] = sum / A[j][j];
    }
  }
  return 0;
}

void choleskylin(double **A, int n, double **L) {
  int i, j, k;
  double sum;

  for (i = 0; i < n; i++) { // L*y=b, acumulating y em x:
    for(j = 0; j<i; j++){
      sum = A[i][j];
      for(k=0; k<j; k++)
        sum -= L[i][k] * L[j][k];
      if(i==j){
        if(sum<=0)
          printf("\nA is not Positive definite!");
        L[i][i]=sqrt(sum);
      }
      L[i][j] = sum / L[i][j];
    }
  }
}

void test()
{
  /* le gránd finale*/

  char command_filename[] = "test_commands.txt";
  FILE *command_unit;
  char data_filename[] = "test_data.txt";
  FILE *data_unit;
  int i;
  int j;
  double *umat;
  double *chmat;
  double cgmat;
  double u_mean;
  double *xmat;
  double *xvec;
  double *ymat;
  double *yvec;
/*
  Specify the spatial grid.
*/
  xvec = subregion( nx, 0.0, 2.0 );

  yvec = subregion( ny, 0.0, 1.0 );
  
  xmat = ( double * ) malloc ( nx * ny * sizeof ( double ) );
  ymat = ( double * ) malloc ( nx * ny * sizeof ( double ) );
  mesh2d(xvec, yvec, xmat, ymat);
/*
  Solve the finite difference approximation to the steady 2D heat equation.
*/
  umat =  heat_equation(xvec, yvec, d, f);
  /*
  cgmat = CGm(umat, xvec, yvec, nx, ny, 1);

  Create a data file.

  data_unit = fopen ( data_filename, "wt" );
  for ( j = 0; j < ny; j++ )
  {
    for ( i = 0; i < nx; i++ )
    {
      fprintf ( data_unit, "  %14g  %14g  %14g  %14g\n", 
        xmat[i+j*nx], ymat[i+j*nx], umat[i+j*nx], cgmat);
    }
    fprintf ( data_unit, "\n" );
  }
  fclose ( data_unit );
  printf ( "\n" );
  printf ( "  Created graphics data file '%s'\n", data_filename );
/*
  Create the command file.

  command_unit = fopen ( command_filename, "wt" );

  fprintf ( command_unit, "# %s\n", command_filename );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "# Usage:\n" );
  fprintf ( command_unit, "#  gnuplot < %s\n", command_filename );
  fprintf ( command_unit, "#\n" );
  fprintf ( command_unit, "set term png\n" );
  fprintf ( command_unit, "set output 'test.png'\n" );
  fprintf ( command_unit, "set xlabel '<---X--->'\n" );
  fprintf ( command_unit, "set ylabel '<---Y--->'\n" );
  fprintf ( command_unit, "set zlabel '<---U(X,Y)--->'\n" );
  fprintf ( command_unit, "set title 'Sample Solution'\n" );
  fprintf ( command_unit, "set contour\n" );
  fprintf ( command_unit, "set cntrparam levels 10\n" );
  fprintf ( command_unit, "set view 75, 75\n" );
  fprintf ( command_unit, "unset key\n" );
  fprintf ( command_unit, "splot '%s'\n", data_filename );

  fclose ( command_unit );

  printf ( "  Created graphics command file '%s'\n", command_filename );
/*
  Report the average value of U.

  u_mean = 0.0;
  for ( j = 0; j < ny; j++ )
  {
    for ( i = 0; i < nx; i++ )
    {
      u_mean = u_mean + umat[i+j*nx];
    }
  }
  u_mean = u_mean / ( double ) ( nx * ny );

  printf ( "\n" );
  printf ( "  Mean value of U is %g\n", u_mean );
/*
  Free memory.
*/
  free ( umat );
  free ( xmat );
  free ( xvec );
  free ( ymat );
  free ( yvec );

  return;
}

int main(int argc, char **argv){

  printf ( "\n" );
  printf ( "PTC - CHOLESKY AND CGm COMPARISON:\n" );
  printf ( "  HPC 1 - BASE CODE \n" );
  printf ( "    Test the heat steady PTC equation library.\n" );

  test( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "HEAT EQUATION SOLVING EXECUTION:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );

  return 0;
}
