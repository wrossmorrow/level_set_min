#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define SGN(x)   ( ((x)==0) ? 0 : ( ((x)>0) ? 1 : -1 ) )
#define MAX(x,y) ( ((x)>(y)) ? (x) : (y) )
#define MIN(x,y) ( ((x)<(y)) ? (x) : (y) )

/* 
 * take an upwind step to solving
 *      ut - Vf.Vu = 0
 * ("V" is grad operator, "." is dot-product)
 * 
 */
void lsostep(double *U, int N, double h, double k , double *F, double *Fx, double *Fy)
{
	int i , j;
	double  Dx , Dy; // numerical derivatives
	
	for( i=0; i<=N ; i++ ) {
		for( j=0 ; j<=N ; j++ ) {
			
			// define directional derivatives, conditional on their existence
			// upwind derivatives come from the gradient direction
			if (*(Fx+i*(N+1)+j)>=0) // non-decreasing in x
				Dx = ( (i==N) ? 0 : (*(U+(i+1)*(N+1)+j) - *(U+i*(N+1)+j))/h );
			else // decreasing in x
				Dx = ( (i==0) ? 0 : (*(U+i*(N+1)+j) - *(U+(i-1)*(N+1)+j))/h );
			
			if (*(Fy+i*(N+1)+j)>=0) // non-decreasing in y
				Dy = ( (i==N) ? 0 : (*(U+i*(N+1)+(j+1)) - *(U+i*(N+1)+j))/h ); 
			else // decreasing in y
				Dy = ( (i==0) ? 0 : (*(U+i*(N+1)+j) - *(U+i*(N+1)+(j-1)))/h );
			
			// update level set function
			*(U+i*(N+1)+j) += k*( (*(Fx+i*(N+1)+j))*Dx + (*(Fy+i*(N+1)+j))*Dy );
		}
	}
}

/* 
 * reinitialize U to be a signed distance function from the interface
 * uses numerical method of Russo & Smereka (2000)
 * 
 */
void reinitialize(double *U, double *V, int N, double h, int S)
{
	int s , i , j;
	double k , D , G , a , b , c , d , ap , am , bp , bm , cp , cm , dp , dm;
	
	// CFL condition
	k = 0.5*h;
	
	// initialize V as U
	for( i=0 ; i<=N; i++ ) for( j=0 ; j<=N; j++ ) 
		*(V+i*(N+1)+j) = *(U+i*(N+1)+j);
	
	// reinitialization iterations
	for( s=1 ; s<=S ; s++ ) {
	
		for( i=0 ; i<=N ; i++ ) {
			for( j=0 ; j<=N ; j++ ) {
				
				if( (i==0)||(i==N)||(j==0)||(j==N) ) { // boundary point
					
					// define directional differences, conditional on their existence
					a = ( (i==0) ? 0 : (*(V+i*(N+1)+j) - *(V+(i-1)*(N+1)+j))/h );
					b = ( (i==N) ? 0 : (*(V+(i+1)*(N+1)+j) - *(V+i*(N+1)+j))/h );
					c = ( (j==0) ? 0 : (*(V+i*(N+1)+j) - *(V+i*(N+1)+(j-1)))/h );
					d = ( (j==N) ? 0 : (*(V+i*(N+1)+(j+1)) - *(V+i*(N+1)+j))/h );
					
					if( *(U+i*(N+1)+j) >= 0 ) {
						ap = MAX(a,0); bm = MIN(b,0);
						cp = MAX(c,0); dm = MIN(d,0);
						G = sqrt( MAX(ap*ap,bm*bm) + MAX(cp*cp,dm*dm) ) - 1;
					} else {
						am = MIN(a,0); bp = MAX(b,0);
						cm = MIN(c,0); dp = MAX(d,0);
						G = sqrt( MAX(am*am,bp*bp) + MAX(cm*cm,dp*dp) ) - 1;
					}
					
					*(V+i*(N+1)+j) -= k*SGN(*(U+i*(N+1)+j))*G;
					
				} else { // interior point
					
					// these conditionals only defined on the interior
					if (  ( *(U+i*(N+1)+j)*(*(U+(i-1)*(N+1)+j)) < 0 ) 
					   || ( *(U+i*(N+1)+j)*(*(U+(i+1)*(N+1)+j)) < 0 ) 
					   || ( *(U+i*(N+1)+j)*(*(U+i*(N+1)+(j-1))) < 0 ) 
					   || ( *(U+i*(N+1)+j)*(*(U+i*(N+1)+(j+1))) < 0 ) ) {
							// within one cell of zero level-set
						
						D = (2*h*(*(U+i*(N+1)+j)))/sqrt( pow( *(U+(i+1)*(N+1)+j) - *(U+(i-1)*(N+1)+j) , 2 )
								+ pow( *(U+i*(N+1)+(j+1)) - *(U+i*(N+1)+(j-1)) , 2 ) );
							
						*(V+i*(N+1)+j) -= (k/h)*( SGN(*(U+i*(N+1)+j))*abs(*(V+i*(N+1)+j)) - D );
						
					} else {
					
						// directional differences
						a = (*(V+i*(N+1)+j) - *(V+(i-1)*(N+1)+j))/h;
						b = (*(V+(i+1)*(N+1)+j) - *(V+i*(N+1)+j))/h;
						c = (*(V+i*(N+1)+j) - *(V+i*(N+1)+(j-1)))/h;
						d = (*(V+i*(N+1)+(j+1)) - *(V+i*(N+1)+j))/h;
						
						if( *(U+i*(N+1)+j) >= 0 ) {
							ap = MAX(a,0); bm = MIN(b,0);
							cp = MAX(c,0); dm = MIN(d,0);
							G = sqrt( MAX(ap*ap,bm*bm) + MAX(cp*cp,dm*dm) ) - 1;
						} else {
							am = MIN(a,0); bp = MAX(b,0);
							cm = MIN(c,0); dp = MAX(d,0);
							G = sqrt( MAX(am*am,bp*bp) + MAX(cm*cm,dp*dp) ) - 1;
						}
						
						*(V+i*(N+1)+j) -= k*SGN(*(U+i*(N+1)+j))*G;
					}
				}
				
			}
		}
	}
}

/*
 * Primary routine
 *
 */
int main (void)
{
	// N is number of "grid cells"; meaning an (N+1)x(N+1) grid
	int i , j , t , N = 25 , T = 500 , s = 10; 
	double h , r , k , z , *x , *y , *F , *Fx , *Fy , *U , *V;
	
	FILE *obj , *lsf;
	
	// output data file for the objective function
	obj = fopen("obj.csv","w");
	// output data file for the level-set function
	lsf = fopen("lsf.csv","w");
	// check for an error with either output file
	if ((obj==NULL)||(lsf==NULL)) { fprintf(stderr,"opening output file failed.\n"); exit(EXIT_FAILURE); }
	
	//////////////////////////////////////////////////////////////////////////////////////////////////
	// ALLOCATIONS ///////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////
	
	// allocate memory for x: grid points (1-d array)
	x = (double *) calloc( sizeof(double) , (N+1) );
	if( x==NULL ) { fprintf(stderr,"calloc failed.\n"); exit(EXIT_FAILURE); }
	// allocate memory for y: grid points (1-d array)
	y = (double *) calloc( sizeof(double) , (N+1) );
	if( y==NULL ) { fprintf(stderr,"calloc failed.\n"); exit(EXIT_FAILURE); }
	
	// allocate memory for F: objective function values over the grid
	F = (double *) calloc( sizeof(double) , (N+1)*(N+1) );
	if( F==NULL ) { fprintf(stderr,"calloc failed.\n"); exit(EXIT_FAILURE); }
	// allocate memory for Fx: partial derivatives, in x-direction, over the grid
	Fx = (double *) calloc( sizeof(double) , (N+1)*(N+1) );
	if( Fx==NULL ) { fprintf(stderr,"calloc failed.\n"); exit(EXIT_FAILURE); }
	// allocate memory for Fy: partial derivatives, in y-direction, over the grid
	Fy = (double *) calloc( sizeof(double) , (N+1)*(N+1) );
	if( Fy==NULL ) { fprintf(stderr,"calloc failed.\n"); exit(EXIT_FAILURE); }
	
	// allocate memory for U: level-set function
	U = (double *) calloc( sizeof(double) , (N+1)*(N+1) );
	if( U==NULL ) { fprintf(stderr,"calloc failed.\n"); exit(EXIT_FAILURE); }
	// allocate memory for V: blank array for reinitializations
	V = (double *) calloc( sizeof(double) , (N+1)*(N+1) );
	if( V==NULL ) { fprintf(stderr,"calloc failed.\n"); exit(EXIT_FAILURE); }
	
	//////////////////////////////////////////////////////////////////////////////////////////////////
	// INITIALIZATION ////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////
	
	// (uniform) mesh size
	h = (double)6.0/N;
	
	// grid points (really just a convenience)
	for( i=0 ; i<=N ; i++ ) *(x+i) = -3 + i*h;
	for( j=0 ; j<=N ; j++ ) *(y+j) = -3 + j*h;
	
	// objective function and its partial derivatives
	for( i=0 ; i<=N ; i++ ) {
		for( j=0 ; j<=N ; j++ ) {
			// objective function
			*(F+i*(N+1)+j)  = pow(*(x+i)-1,2) + pow(*(y+j)-1,2);
			// partial derivative in x-direction
			*(Fx+i*(N+1)+j) = 2*(*(x+i)-1);
			// partial derivative in y-direction
			*(Fy+i*(N+1)+j) = 2*(*(y+j)-1);
		
			// update heuristic CFL denominator with the norm of the 
			// gradient of the objective function (see write-up)
			if ((i==0)&&(j==0)) z = sqrt( pow(*(Fx+i*(N+1)+j),2) + pow(*(Fy+i*(N+1)+j),2) );
			else z = MAX( sqrt( pow(*(Fx+i*(N+1)+j),2) + pow(*(Fy+i*(N+1)+j),2) ) , z );
		}
	}
	
	printf("%f\n",h);
	printf("%f\n",z);
	
	// output objective function values
	for( i=0 ; i<=N ; i++ ) {
		for( j=0 ; j<=N ; j++ ) {
			if (j==N) fprintf(obj,"%f",*(F+i*(N+1)+j));
			else fprintf(obj,"%f,",*(F+i*(N+1)+j));
		}
		fprintf(obj,"\n");
	}
	fprintf(obj,"\n");
	// partial derivative in x-direction
	for( i=0 ; i<=N ; i++ ) {
		for( j=0 ; j<=N ; j++ ) {
			if (j==N) fprintf(obj,"%f",*(Fx+i*(N+1)+j));
			else fprintf(obj,"%f,",*(Fx+i*(N+1)+j));
		}
		fprintf(obj,"\n");
	}
	fprintf(obj,"\n");
	// partial derivative in y-direction
	for( i=0 ; i<=N ; i++ ) {
		for( j=0 ; j<=N ; j++ ) {
			if (j==N) fprintf(obj,"%f",*(Fy+i*(N+1)+j));
			else fprintf(obj,"%f,",*(Fy+i*(N+1)+j));
		}
		fprintf(obj,"\n");
	}
	
	// time-step parameter via *heuristic* CFL condition
	k = 0.5*h/z;
	
	// radius of zero-level set circle
	r = 1.75;
	
	// initialize U to have the circle of radius r as a zero-level set 
	for( i=0 ; i<=N ; i++ ) for( j=0 ; j<=N ; j++ ) 
		*(U+i*(N+1)+j) = pow(*(x+i),2) + pow(*(y+j),2) - pow(r,2);
		
	// store initial level-set function
	for ( i=0 ; i<=N ; i++ ) {
		for ( j=0 ; j<=N ; j++ ) {
			if (j==N) fprintf(lsf,"%f",*(U+i*(N+1)+j));
			else fprintf(lsf,"%f,",*(U+i*(N+1)+j));
		}
		fprintf(lsf,"\n");
	}
	fprintf(lsf,"\n");
	
	//////////////////////////////////////////////////////////////////////////////////////////////////
	// ITERATIONS OF THE LEVEL-SET FUNCTION //////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////
	
	for( t=1 ; t<=T ; t++ ) {
		
		// take an upwind step
		lsostep(U,N,h,k,F,Fx,Fy);
		
		/*
		// reinitialize U to be a signed distance function from the zero-level set
		// (U is needed in its "original form" in this algorithm)
		reinitialize(U,V,N,h,100);
		
		// replace U with V, a reinitialized level-set function
		for ( i=0 ; i<=N ; i++ ) for ( j=0 ; j<=N ; j++ )
			*(U+i*(N+1)+j) = *(V+i*(N+1)+j);
		*/
		
		// store data at appropriate intervals
		if(fmod(t,s)==0) {
			for ( i=0 ; i<=N ; i++ ) {
				for ( j=0 ; j<=N ; j++ ) {
					if (j==N) fprintf(lsf,"%f",*(U+i*(N+1)+j));
					else fprintf(lsf,"%f,",*(U+i*(N+1)+j));
				}
				fprintf(lsf,"\n");
			}
			// include this return so that data from separate
			// iterations can be understood as such
			fprintf(lsf,"\n");
		}
		
	}
	
	//////////////////////////////////////////////////////////////////////////////////////////////////
	// CUSTODIAL /////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////
	
	// release dynamically allocated memory (V should point to the beginning of this block)
	free(x);  // grid points (1-d array)
	free(y);  // grid points (1-d array)
	free(F);  // objective function values over the grid
	free(Fx); // partial derivatives of the objective, in x-direction, over the grid
	free(Fy); // partial derivatives of the objective, in y-direction, over the grid
	free(U);  // level-set function
	free(V);  // blank array for reinitializations
	
	// close output files
	fclose(obj);
	fclose(lsf);
	
	return 0;
}