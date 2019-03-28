#include <mpi.h>
#include <stdio.h>
#include <math.h>

int main(int argc,char **argv)
{
    int n = 0, myid, numprocs, i, root = 0, tp;
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
	
	if(myid==root){
		//Crear matriz A y B, con dimensiÃ³n nxn dada por el usuario.
		//el tamaño de las matrices debe ser multiplo de la cantidad de procesos, siendo k 
		// la cantidad de procesos, se le da n/k filas a cada proceso
		//ej. n = 1000, k = 100 => cada proceso maneja 10 filas.
		//A con valores random entre 0 y 5, B con valores entre 0 y 2.
	}
    MPI_Barrier(MPI_COMM_WORLD);

	/*
		Calcular de forma distribuida entre todos los procesos:
			- Matriz M = AxB.
			- tp = Numero total de valores primos en M
			- Vector P tal que P[i] = a la cantidad de primos en la columna i de M.
			- Crear matriz C tal que:
				C[i,j] = M[i,j] + M[i+1,j] + M[i-1,j] + M[i,j+1] + M[i,j-1]
					(Considerar filas y columnas extremas de la matriz.)
	*/

	if (myid == root){
        
	} else {
		
	}

           
    MPI_Finalize();
    return 0;
}