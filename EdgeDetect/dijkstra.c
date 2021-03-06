/*
* 
* Implementation of Dijkstra's method using a Matlab sparse matrix
* as an adjacency matrix.  Zero entries represent non-existent edges.
* Uses linear search for simplicity 
*
 Inputs:
 -------

   A        - sparse adjacency matrix (N x N)
   s        - Source node index in [1,...,N]
   d        - Destination node index in [1,...,N]


   
 Outputs:
 -------

   path     - distance vector from Dijkstra (1 x m)
   pathcost - Cost of the path



Example 1
---------


A                            = sparse([0 1 0 0 0 0 ; 1 0 1 0 0 1 ; 0 1 0 1 0 0 ; 0 0 1 0 1 1 ; 0 0 0 1 0 1 ; 0 1 0 1 1 0]);
s                            = 1;
d                            = 5;
[path , pathcost]            = dijkstra(A , s , d);




Example 1
---------

close all

N                            = 2000;
L                            = 1000;
R                            = 2*L/sqrt(N);%200;
s                            = 1;
d                            = 10;
X                            = L*rand(2 , N);

A                            = (Radjacency(X  ,R));
%A(A~=0)                      = 1;
[path , pathcost]            = dijkstra(A , s , d);
hold on,h=plot(X(1 , :) , X(2 , :) , '+' , X(1 , path) , X(2 , path) , 'r-+', X(1 , s) , X(2 , s) , 'ko' , X(1 , d) , X(2 , d) , 'mo' , 'linewidth' , 3);,hold off

legend(h(2:3) , 'Dijkstra')





mex -g -output dijkstra.dll dijkstra.c

mex -f mexopts_intel10amd.bat -output dijkstra.dll dijkstra.c


 Author : Sébastien PARIS : sebastien.paris@lsis.org
 -------  Date : 11/01/2007

 Reference ""


 *
 *Change Log
 *2019/08/26 RML changed mxErrMsgTxt to mexErrMsgTxt (running in MATLAB 2019a)
*/

#include <math.h>
#include <mex.h>




void mexFunction(int nlhs, mxArray** plhs ,  int nrhs, const mxArray** prhs)
{
    
	double *A;
	
	double *path , *pathcost;
	
	double *pathtemp;
	
	int *jc, *ir;
	
	int i , j , n , s , d , m=0 , u , v , t  , p;
	
	
    double du , inf = mxGetInf() , dist , Auv;
	
	
	double *distance;
	
	int *parent, *visited;
	
    if (nrhs != 3)
	{
        
		mexErrMsgTxt("Too few arguments");
		
	}
    
	/* Input 1 */
	
	if (!mxIsSparse(prhs[0]) || mxGetM(prhs[0]) != mxGetN(prhs[0]))
	{
        mexErrMsgTxt("Graph must be a square matrix");
	}
	
    A        = mxGetPr(prhs[0]);
	
	jc       = mxGetJc(prhs[0]);
	
	ir       = mxGetIr(prhs[0]);
	
    n        = mxGetN(prhs[0]);
	
	/* Input 2 */
	
    s        = ((int) mxGetScalar(prhs[2])) - 1;
	
    if (s < 0 || s > n)
	{
        mexErrMsgTxt("Source identifier is out of range");
	}
	
	/* Input 3 */
	
    d         = ((int) mxGetScalar(prhs[1])) - 1;
	
    if (d < 0 || d > n)
	{
        mexErrMsgTxt("Destination identifier is out of range");
	}
	
	
	
    distance  = mxMalloc(n*sizeof(double));
	
	pathtemp  = mxMalloc(sizeof(double));
	
    parent    = mxMalloc(n*sizeof(int));
	
    visited   = mxMalloc(n*sizeof(int));
	
	
	
	/* Main Call */
	
	
    plhs[1]   = mxCreateDoubleMatrix(1, 1, mxREAL);
	
    pathcost  = mxGetPr(plhs[1]);
	
	
	
	
	for (i = 0 ; i < n ; i++)
	{
		
		distance[i] = inf;
		
		parent[i]   = 0;
		
		visited[i]  = 0;
			
	}
	
	
	distance[s]     = 0.0;
	
	for (i = 0 ; i < n - 1 ; i++)
	{
		
		
		u        = 0;
		
		du       = inf;
		
		
		for (j = 0 ; j < n ; j++)
		{
			
			if(visited[j] == 0)
			{
				
				dist     = distance[j];
				
				if(dist < du)
					
				{
					
					du       = dist;
					
					u        = j;
					
				}
			}
			
		}
		
		if (u == d)
		{
			
			break;
			
		}
		
		
		visited[u]  = 1;  
		
		for (j = jc[u] ; j < jc[u + 1] ; j++) 
		{
			v   = ir[j];
			
			Auv = A[j]; 
			
			if ( (distance[v] > du + Auv) ) // (Auv != 0.0) &&
			{
				distance[v] = du + Auv;
				
				parent[v]   = u;
			}
			
		}  	   
		
	}
	
	
	/* BackTracking */
	
	pathcost[0] = distance[d];
	
	if(parent[d] != 0)
		
	{
		
		t           = d;
		
		pathtemp[0] = (double)(d + 1);
		
		
		while (t != s)
		{
			
			
			pathtemp    = mxRealloc(pathtemp , (m+2)*sizeof(double)); 
			
			p           = parent[t];
			
			m++;
			
			pathtemp[m] = (double)(p + 1);
			
			t           = p;
			
			
		}
		
		m++;
	}
	
	
	
	/* Ouputs */
	
    plhs[0]   = mxCreateDoubleMatrix(m , 1 , mxREAL);
	
	path      = mxGetPr(plhs[0]);
	
	
	for (i = 0 ; i < m ; i++)
	{
		
		path[i] = pathtemp[i];
		
	}
	
	
	mxFree(distance);
	
	mxFree(parent);
	
	mxFree(visited);
	
	mxFree(pathtemp);
	
}


