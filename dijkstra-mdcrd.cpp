// A C / C++ program for Dijkstra's single source shortest
// path algorithm. The program is for adjacency matrix
// representation of the graph.
// taken from http://www.geeksforgeeks.org/printing-paths-dijkstras-shortest-path-algorithm/
#include <stdio.h>
#include <limits.h>
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <cmath>
using namespace std;


// Number of vertices in the graph

// A utility function to find the vertex with minimum distance
// value, from the set of vertices not yet included in shortest
// path tree
int minDistance(int dist[], bool sptSet[], int V)
{
	// Initialize min value
  	int min = INT_MAX, min_index;

	for (int v = 0; v < V; v++)
		if (sptSet[v] == false && dist[v] <= min)
			min = dist[v], min_index = v;

	return min_index;
}

// Function to print shortest path from source to j
// using parent array 
// prints atom numbers
void printPath(int parent[], int j, int numat[])
{
	// Base Case : If j is source
	if (parent[j]==-1)
		return;

	printPath(parent, parent[j],numat);

	printf("%d ", numat[j]);
}

// A utility function to print the constructed distance
// array
int printSolution(int dist[], int n, int parent[], int numat[])
{
	int src = 0;
	printf("Vertex\t Distance\tPath");
	for (int i = 1; i < n; i++)
	{
		printf("\n%d -> %d \t\t %d\t\t%d ", src, i, dist[i], src);
		printPath(parent, i, numat);
		printf("\n");
	}
}

// Funtion that implements Dijkstra's single source shortest path
// algorithm for a graph represented using adjacency matrix
// representation
void dijkstra(int src , int V, int *graph, int numat[])
{
	int dist[V]; // The output array. dist[i] will hold
				// the shortest distance from src to i

	// sptSet[i] will true if vertex i is included / in shortest
	// path tree or shortest distance from src to i is finalized
	bool sptSet[V];

	// Parent array to store shortest path tree
	int parent[V];

	// Initialize all distances as INFINITE and stpSet[] as false
	for (int i = 0; i < V; i++)
	{
		parent[0] = -1;
		dist[i] = INT_MAX;
		sptSet[i] = false;
	}

	// Distance of source vertex from itself is always 0
	dist[src] = 0;

	// Find shortest path for all vertices
	for (int count = 0; count < V-1; count++)
	  {
		// Pick the minimum distance vertex from the set of
		// vertices not yet processed. u is always equal to src
		// in first iteration.
	    int u = minDistance(dist, sptSet, V);

		// Mark the picked vertex as processed
	    sptSet[u] = true;

		// Update dist value of the adjacent vertices of the
		// picked vertex.
	    for (int v = 0; v < V; v++)
	      if (!sptSet[v] && graph[u+v*V] &&
		  dist[u] + graph[u+v*V] < dist[v])
		{
		  parent[v] = u;
		  dist[v] = dist[u] + graph[u+v*V];
		} 
	  }

	// print the constructed distance array
	printSolution(dist, 3, parent,numat);
}

double distance (double xi, double xf, double yi, double yf,double zi,double zf)
{ 
 // Function to calculate distance
  double d,dummy1,dummy2,dummy3,dummy4;
  dummy1= (xf-xi)*(xf-xi);
  dummy2=(yf-yi)*(yf-yi);
  dummy3= (zf-zi)*(zf-zi);
  dummy4=dummy1+dummy2+dummy3;
  d=sqrt(dummy4);
  return d;
}


// driver program to test above function
int main()
{
	/* Heavily edited by I. Sumner for water wire problem*/
	  string name1,name2,name3,name4;
	  int * numat, numres, natoms, * graph, nstart,nfinish1,nfinish2,numatdum[3],iline,nwater;
	  double * x,* y,* z,dum1,dum2,dist1,dist2,dist3,cutoff,xdum[3],ydum[3],zdum[3];
	  string filename;
	  //read in file name, expects mdcrd format
	  cin >> filename;
  ifstream inputfile;
  inputfile.open(filename.c_str());
  string line;
  int i=0;
  iline=0;
  //read in number of atoms
  cin >> natoms;
  x =  new double[3*natoms];
  y =  new double[natoms];
  z =  new double[natoms];
  numat =  new int[natoms];
  int numlines = (3*natoms)%10;
  int nat = numlines;
  if (numlines > 0){
    numlines = int(3*natoms/10.0)+1;
  }
  else{
    numlines =  int(3*natoms/10.0);
  }
  //read in endpoints for pathways
  cin >> nstart;
  cin >> nfinish1;
  cin >> nwater;
  cin >> cutoff;
  //  cout << "numlines " << numlines << "\n";
  int noxygens=(natoms-nwater-1)/3;
  getline(inputfile,line);
  while (getline(inputfile,line))
    {
      //only read in lines with cartesian coordinates
      if( i+1 != numlines  ){
	istringstream ss(line);
	ss >> x[i*10] >> x[i*10+1] >> x[i*10+2] >> x[i*10+3] >> x[i*10+4] >> x[i*10+5] >> x[i*10+6] >> x[i*10+7] >> x[i*10+8] >> x[i*10+9];
	i++;
      }
      else if(i+1 == numlines){
	iline++;
	istringstream ss(line);
	int j = 0;
	for(j=0;j<nat;j++){
	  ss>>x[i*10+j];
	}	
	xdum[1]=x[(nfinish1-1)*3];
	ydum[1]=x[(nfinish1-1)*3+1];
	zdum[1]=x[(nfinish1-1)*3+2];
	xdum[0]=x[(nstart-1)*3];
	ydum[0]=x[(nstart-1)*3+1];
	zdum[0]=x[(nstart-1)*3+2];
	for(j=0;j<=noxygens;j++){
	  x[3+j]=x[(nwater-1)*3+j*9];
	  y[3+j]=x[(nwater-1)*3+j*9+1];
	  z[3+j]=x[(nwater-1)*3+j*9+2];
	  numat[3+j]=(nwater-1)+j*3+1;
	}
	x[0]=xdum[0];
	y[0]=ydum[0];
	z[0]=zdum[0];
	x[1]=xdum[1];
	y[1]=ydum[1];
	z[1]=zdum[1];
	numat[0]=nstart;
	numat[1]=nfinish1;
	i++;
	natoms=noxygens+3;
	j =0;
	  //cut down on number of waters. Only consider waters within 15 ang of endpoints
	  for(i=0;i<natoms;i++){
	    dist1=distance(x[0],x[i],y[0],y[i],z[0],z[i]);
	    dist2=distance(x[1],x[i],y[1],y[i],z[1],z[i]);
	    if((dist1<15)||(dist2<15)){
	      x[j]=x[i];
	      y[j]=y[i];
	      z[j]=z[i];
	      numat[j]=numat[i];
	      j++;
	    }
	  }
	//only viable pathways are those with water oxygens within cutoff
	//If the distance is less than the cutoff, the distance is set to 1
        //If the distance is greater than the cutoff, the distance is set to 1000
	  natoms=j;
	  graph = new int[natoms*natoms];
	  cout << "geom " << iline << "\n";
	  //diagonal portion of distance graph
	  for(i=0;i<natoms;i++){
	    graph[i+i*natoms]=1;
	  }
	  //off diagonal 
	  for(i=0;i<natoms;i++){
	    for(j=0;j<i;j++){
	      dist1=distance(x[j],x[i],y[j],y[i],z[j],z[i]);
	      graph[i+j*natoms]=1;
	      graph[j+i*natoms]=1;
	      if(dist1>cutoff){
		graph[i+j*natoms]=1000;
		graph[j+i*natoms]=1000;
	      }
	    }
	  }
       //  The code could be sped up greatly if only one of the off-diagonal blocks is considered
       // This is a symmetric matrix
	  dijkstra(0,natoms,graph,numat);
	  getline(inputfile,line);
	  i=0;
      }
    }
  inputfile.close();
	//	return 0;
}
