#include<mpi.h>
#include<iostream>
#include<random>
#include<math.h>
#include<bits/stdc++.h>
using namespace std;



int main(int argc,char* argv[])
{
    MPI_Init(&argc, &argv);

    int p, rank;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_size(comm,&p);
    MPI_Comm_rank(comm,&rank);
    cout<<"Hello world!";
    MPI_Finalize();

}
