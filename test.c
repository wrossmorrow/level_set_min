#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define SQ(x)    x*x
#define SGN(x)   ( ((x)==0) ? 0 : ( ((x)>0) ? 1 : -1 ) )
#define MAX(x,y) ( ((x)>(y)) ? (x) : (y) )
#define MIN(x,y) ( ((x)<(y)) ? (x) : (y) )

void reinitialize(double *U, double *V, int N, double h, int S)
{
	int s , i , j;
	double k = 0.5*h , D , G , a , b , c , d;
	
	// initialize V
	for( i=0 ; i<=N; i++ ) for( j=0 ; j<=N; j++ ) 
		*(V+i*(N+1)+j) = *(U+i*(N+1)+j);
	
	// reinitialization iterations
	for( s=1 ; s<=S ; s++ ) {
	
		for( i=0 ; i<=N ; i++ ) {
			for( j=0 ; j<=N ; j++ ) {
				
				if( (i==0)||(i==N)||(j==0)||(j==N) ) { 
					// boundary point
					
					if (i>0) a = (*(V+i*(N+1)+j) - *(V+(i-1)*(N+1)+j))/h; else a = 0;
					if (i<N) b = (*(V+(i+1)*(N+1)+j) - *(V+i*(N+1)+j))/h; else b = 0;
					if (j>0) c = (*(V+i*(N+1)+j) - *(V+i*(N+1)+(j-1)))/h; else c = 0;
					if (j<N) d = (*(V+i*(N+1)+(j+1)) - *(V+i*(N+1)+j))/h; else d = 0;
					
					if( *(U+i*(N+1)+j) >= 0 )
						G = sqrt(MAX(SQ(MAX(a,0)),SQ(MIN(b,0)))+MAX(SQ(MAX(c,0)),SQ(MIN(d,0))));
					else
						G = sqrt(MAX(SQ(MIN(a,0)),SQ(MAX(b,0)))+MAX(SQ(MIN(c,0)),SQ(MAX(d,0))));
					
					*(V+i*(N+1)+j) -= k*SGN(*(U+i*(N+1)+j))*G;
					
				} else {
					// interior point
					
					if ( ( *(U+i*(N+1)+j)*(*(U+(i-1)*(N+1)+j)) < 0 ) 
							|| ( *(U+i*(N+1)+j)*(*(U+(i+1)*(N+1)+j)) < 0 ) 
							|| ( *(U+i*(N+1)+j)*(*(U+i*(N+1)+(j-1))) < 0 ) 
							|| ( *(U+i*(N+1)+j)*(*(U+i*(N+1)+(j+1))) < 0 ) ) {
						
						D = (2*h*(*(U+i*(N+1)+j)))/sqrt( SQ( *(U+(i+1)*(N+1)+j) - *(U+(i-1)*(N+1)+j) )
								+ SQ( *(U+i*(N+1)+(j+1)) - *(U+i*(N+1)+(j-1)) ) );
							
						*(V+i*(N+1)+j) -= (k/h)*( SGN(*(U+i*(N+1)+j))*abs(*(V+i*(N+1)+j)) - D );
						
					} else {
					
						a = (*(V+i*(N+1)+j) - *(V+(i-1)*(N+1)+j))/h;
						b = (*(V+(i+1)*(N+1)+j) - *(V+i*(N+1)+j))/h;
						c = (*(V+i*(N+1)+j) - *(V+i*(N+1)+(j-1)))/h;
						d = (*(V+i*(N+1)+(j+1)) - *(V+i*(N+1)+j))/h;
						
						if( *(U+i*(N+1)+j) >= 0 )
							G = sqrt(MAX(SQ(MAX(a,0)),SQ(MIN(b,0)))+MAX(SQ(MAX(c,0)),SQ(MIN(d,0))));
						else
							G = sqrt(MAX(SQ(MIN(a,0)),SQ(MAX(b,0)))+MAX(SQ(MIN(c,0)),SQ(MAX(d,0))));
						
						*(V+i*(N+1)+j) -= k*SGN(*(U+i*(N+1)+j))*G;
					}
				}
				
			}
		}
	}
}

int main (void)
{
	int i , j , N=6; // N is number of "grid cells"; meaning an (N+1)x(N+1) grid
	double *U , *V;
	
	FILE *output;	
	
	// allocate memory for U
	U = (double *) calloc( sizeof(double) , (N+1)*(N+1) );
	if( U==NULL ) { fprintf(stderr,"calloc failed.\n"); exit(EXIT_FAILURE); }
	
	// allocate memory for V
	V = (double *) calloc( sizeof(double) , (N+1)*(N+1) );
	if( V==NULL ) { fprintf(stderr,"calloc failed.\n"); exit(EXIT_FAILURE); }
	
	// output data file
	output = fopen("data.csv","w");
	if(output==NULL) { fprintf(stderr,"opening output file failed.\n"); exit(EXIT_FAILURE); }
	
	// initialize U
	for( i=0 ; i<=N ; i++ ) for( j=0 ; j<=N ; j++ ) 
		*(U+i*(N+1)+j) = (double)(i*j);
	
	// print
	for ( i=0 ; i<=N ; i++ ) { 
		for ( j=0 ; j<=N ; j++ ) { 
			printf(" %f ",*(U+i*(N+1)+j));
		} 
		printf("\n"); 
	}
	
	// reinitialize U to be a signed distance function from the interface
	// uses numerical method of Smereka & Russo [1]
	reinitialize(U,V,N,0.1,2);
	
	// replace U with V
	for ( i=0 ; i<=N ; i++ ) {
		for ( j=0 ; j<=N ; j++ ) {
			// printf(" %f ",*(V+i*(N+1)+j));
			if (j==N) fprintf(output,"%f",*(V+i*(N+1)+j));
			else fprintf(output,"%f,",*(V+i*(N+1)+j));
			// *(U+i*(N+1)+j) = *(V+i*(N+1)+j);
			// printf(" %f ",*(U+i*(N+1)+j));
		}
		// printf("\n");
		fprintf(output,"\n");
	}
	printf("\n");
	
	// release dynamically allocated memory (V should point to the beginning of this block)
	free(U);
	free(V);
	
	// close output file
	fclose(output);
	
	return 0;
}