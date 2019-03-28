#include <mpi.h>
#include <stdio.h>
#include <math.h>

int main(int argc,char **argv)
{
    int n = 0, myid, numprocs, i, root = 0;
    double startwtime, endwtime;
    int  namelen;
    char processor_name[MPI_MAX_PROCESSOR_NAME];

    MPI_Init(&argc,&argv);
/*  Se inicia el trabajo con MPI */
           
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
/*  MPI almacena en numprocs el numero total de procesos que se pusieron a correr */
 
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);
/*  MPI almacena en myid la identificacion del proceso actual */

    MPI_Get_processor_name(processor_name,&namelen);
/*  MPI almacena en processor_name el nombre de la computadora en la que corre el
    proceso actual, y en namelen la longitud de este */

    MPI_Barrier(MPI_COMM_WORLD);

	if (myid == root){
        
	} else {
		
	}

           
    MPI_Finalize();
    return 0;
}