#include<mpi.h>
#include<iostream>
#include<random>
#include<math.h>
#include<bits/stdc++.h>
using namespace std;
int p, rank_process;

#define PI 3.14159265

//Helper function to check if point within circle
bool isWithinCircle(double x, double y, double R)
{
    return (abs(x)<= R/sqrt(2) && abs(y)<=R/sqrt(2))?true:false;

}

//Helper function to calculate the sum of the vector
double calculate_sum(vector<double> &buf)
{
    double temp_sum = 0;
    for(int i=0;i<buf.size();i++)
    {
            temp_sum += buf[i];
    }

    return temp_sum;

}


int dboard(int N)
{
   int darts_per_processor, remainder, rank_iterator, counter;
   double R;
   vector<int> buf;
   remainder = N%p;
   if(rank_process < remainder){ darts_per_processor = N/p + 1;}
   else if(rank_process >= remainder){darts_per_processor = N/p;}

   //Darting simulation :
   // Radius of the circle is R=10
   //Assuming r=(0,10) and theta=(0,360)[In degree]
   counter = 0;
   R = 10.00;


   for(int i=0;i<darts_per_processor;i++)
   {
        //Caveat: The distriution for r is not uniform.
        double r = R * sqrt((rand()%101)/100.00);
        double theta = rand()%361;

        double x = r*cos(theta * PI/180.00);
        double y = r*sin(theta * PI/180.00);

        if(isWithinCircle(x,y,R)){counter += 1;}

   }
    return counter;
}


int main(int argc,char* argv[])
{
    MPI_Init(&argc, &argv);
    double t_start = MPI_Wtime();
    MPI_Comm comm = MPI_COMM_WORLD;

    MPI_Comm_size(comm,&p);
    MPI_Comm_rank(comm,&rank_process);

    //Setting values and declarations
    if(argc < 3)
    {
        cout<<"Please insert the values of N and R";
        exit(1);
    }
    //Setting random seed
    srand(rank_process);

    int N,R,darts_in_square,sum_processes,size_of_subarray;
    double sum_simulations;
    vector<double> results;

    vector<double>::iterator ite_begin, ite_end;

    if(rank_process == 0){N = stoi(argv[1]); R = stoi(argv[2]);}

    //Broadcasting N and R
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&R, 1, MPI_INT, 0, MPI_COMM_WORLD);

    //Dart simulation for R iterations
    //Vector results has the values across different iterations.
    for(int j=0;j<R;j++)
    {
        darts_in_square = dboard(N);
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Reduce(&darts_in_square, &sum_processes, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

        if(rank_process==0)
        {
            double calculated_pi = 2 * (double)N/sum_processes;
            results.push_back(calculated_pi);
        }

    }

    vector<double> buf(R/p);

    //Scattering and reducing 'results' vector for faster computation
    MPI_Scatter(&results[0], R/p, MPI_DOUBLE, &buf[0], R/p, MPI_DOUBLE, 0,MPI_COMM_WORLD);
    double temp_sum = calculate_sum(buf);
    double total_sum;
    MPI_Reduce(&temp_sum, &total_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    double t_end = MPI_Wtime();
    if(rank_process == 0)
    {
        cout<<"N = "<<N<<"  R = "<<R<<"  P = "<<p<<"  PI = "<<(double)total_sum/R;
        cout<<"\nTime taken "<<t_end-t_start<<"s";
    }

    MPI_Finalize();
    return 0;

}
