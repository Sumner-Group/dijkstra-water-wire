// A C / C++ program for Dijkstra's single source shortest
// path algorithm. The program is for adjacency matrix
// representation of the graph.
// taken from http://www.geeksforgeeks.org/printing-paths-dijkstras-shortest-path-algorithm/
// Many of these subroutines were edited by I. Sumner for the water wire problem.
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
	int i = 1;
	// Even though all paths are computed, we only print 
	// out the distance from donor to acceptor,
	// since that is all we care about. 
	if(dist[i]<1000){
	  printf("\t Start -> End\t Number of Waters\tPath");
	  printf("\n\t %d -> %d \t\t %d\t\t%d ", numat[0], numat[1], dist[i]-1, numat[0]);
	  printPath(parent, i, numat);
	  printf("\n");
	}
	else {
	  printf("No Viable Path\n");
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
	printSolution(dist, 2, parent,numat);
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


int mdcrd(string filename, double x[], double y[], double z[],double xdum[], double ydum[], double zdum[], int numat[], int numlines, int iline, int nat, int nfinish1, int nstart, int noxygens, int nwater, int natoms, double cutoff)
{
  double dist1,dist2;
  int * graph ;
  ifstream inputfile;
  inputfile.open(filename.c_str());
  string line;
  getline(inputfile,line);
  int i = 0;
  while (getline(inputfile,line))
    {
      //read in coordinates
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
	//saving water oxygen coordinates. Assuming that the atom order for each water molecule is:
	// oxygen, hydrogen, hydronge
	for(j=0;j<=noxygens;j++){
	  x[2+j]=x[(nwater-1)*3+j*9];
	  y[2+j]=x[(nwater-1)*3+j*9+1];
	  z[2+j]=x[(nwater-1)*3+j*9+2];
	  //numat array stores the atom numbers of the donor, acceptor and water oxygens
	  numat[2+j]=(nwater-1)+j*3+1;
	}
	//save donor/acceptor coordinates and atom numbers
	x[0]=xdum[0];
	y[0]=ydum[0];
	z[0]=zdum[0];
	x[1]=xdum[1];
	y[1]=ydum[1];
	z[1]=zdum[1];
	numat[0]=nstart;
	numat[1]=nfinish1;
	i++;
	natoms=noxygens+2;
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
	  printf("geom %d ",iline);
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
	  //This algorithm will calculate the path from each atom stored
	  //to every other atom. However, it only prints out donor/acceptor path
	  dijkstra(0,natoms,graph,numat);
	  getline(inputfile,line);
	  i=0;
      }
    }
  inputfile.close();
  return 0;
}

int pdb(string filename, double x[], double y[], double z[],double xdum[], double ydum[], double zdum[], int numat[], int iline, int nfinish1, int nstart, int nwater, int natoms, double cutoff)
{
  double dist1,dist2, dum1, dum2;
  int * graph, numres, numatdum[3];
  string name1, name2, name3, name4;

  ifstream inputfile;
  inputfile.open(filename.c_str());
  string line;
  getline(inputfile,line);
  int i = 0;
  while (!inputfile.eof())
    {
      getline(inputfile,line);

      //only read in lines with cartesian coordinates
      //This is sloppy, but works in my test cases.
      if(line.length() == 80 ){
	istringstream ss(line);
	ss >> name1 >> numat[i] >> name2 >> name3 >> numres >> x[i] >> y[i] >> z[i] >> dum1 >> dum2 >> name4;  
	//only store water oxygens and endpoints.
	//Clearly, the code assums that waters are labeld as "WAT" and the oxygen as "O"
	if(name3=="WAT"&&name2=="O"){
	    i++;
	  }
	  else if(numat[i]==nfinish1){
	    xdum[1]=x[i];
	    ydum[1]=y[i];
	    zdum[1]=z[i];
	    numatdum[1]=numat[i];
	    i++;
	  }
	  else if(numat[i]==nstart){
	    xdum[0]=x[i];
	    ydum[0]=y[i];
	    zdum[0]=z[i];
	    numatdum[0]=numat[i];
	    i++;
	  }
      }
      else if(line.length() < 10 ){
  	istringstream ss(line);
  	ss >> name1;
	//Assuming that the pdb frames end with "ENDMDL" This is the case for using AmberTools to convert
	//a trajectory file to a pdb file.
  	if( name1 == "ENDMDL"){
	  iline++;
  	  name1 = "dummy";
	  x[0]=xdum[0];
	  y[0]=ydum[0];
	  z[0]=zdum[0];
	  x[1]=xdum[1];
	  y[1]=ydum[1];
	  z[1]=zdum[1];
	  numat[0]=numatdum[0];
	  numat[1]=numatdum[1];
	  natoms=i;
	  int j =0;
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
	  natoms=j;
	  graph = new int[natoms*natoms];
	  printf("geom %d ",iline);
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
	  //This algorithm will calculate the path from each atom stored
	  //to every other atom. However, it only prints out donor/acceptor path
	  dijkstra(0,natoms,graph,numat);
	  i=0;
	}
      }
    }
  inputfile.close();
  return 0;
}
// driver program to test above function
int main()
{
	/* Heavily edited by I. Sumner for water wire problem*/
	  string name1,name2,name3,name4;
	  int * numat, numres, natoms, * graph, nstart,nfinish1,nfinish2,numatdum[3],iline,nwater;
	  double * x,* y,* z,dum1,dum2,dist1,dist2,dist3,cutoff,xdum[3],ydum[3],zdum[3];
	  string filename,type;
	  cin >> type;
	  //read in file name, expects mdcrd format
	  cin >> filename;
  int i=0;
  iline=0;
  //read in number of atoms
  cin >> natoms;
  x =  new double[3*natoms];
  y =  new double[natoms];
  z =  new double[natoms];
  numat =  new int[natoms];
  //assumes standard mdcrd file form with 10 entries per line
  if(type=="mdcrd"){
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
    int noxygens=(natoms-nwater-1)/3;
    int done = mdcrd(filename, x, y, z, xdum, ydum, zdum, numat, numlines, iline, nat, nfinish1, nstart, noxygens, nwater, natoms, cutoff);
  }
  else if(type=="pdb"){
    cin >> nstart;
    cin >> nfinish1;
    cin >> cutoff;
    int done = pdb(filename, x, y, z, xdum, ydum, zdum, numat, iline, nfinish1, nstart, nwater, natoms, cutoff);
  }
  return 0;
}
