#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


/*
	Matrices como vector:
	-M[k][n] = M[k*n]
	-M[fila][columna] = M[fila*(dimension de la fila)+columna]
*/

/*
	Llena la matriz recibida por vector con numeros aleatorios en el rango de los
	enteros recibidos por par�metros.
	
	Para el proposito del programa, el valor m�nimo es cero siempre.
*/
void LlenarMatriz(int matriz[],int dimensionMatriz, int valorMax){
	time_t t;
	srand(time( &t ));
	for(int i= 0; i < dimensionMatriz; i++){
		matriz[i] = (int)rand() % (valorMax+1);
	}
}


int main(int argc,char **argv)
{
    int n = 0, myid, numprocs, i, root = 0, tp;
    double startwtime, endwtime;
    int  namelen;
	int *A, *B, *M, *C, *P;//Todas las matrices que se van a manejar.
	//P contiene la cantidad de primos por fila.
    char processor_name[MPI_MAX_PROCESSOR_NAME];
		//Solicitar n para las dimensiones de las matrices.
		//Crear matriz A y B, con dimensión nxn dada por el usuario.
		//el tamano de las matrices debe ser multiplo de la cantidad de procesos, siendo k 
		// la cantidad de procesos, se le da n/k filas a cada proceso
		//ej. n = 1000, k = 100 => cada proceso maneja 10 filas.
		//A con valores random entre 0 y 5, B con valores entre 0 y 2.
		
		
	
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
	
	if(myid == root){
		printf("Digite el valor para las dimensiones de las matrices: \n");
		scanf("%d",&n);
		printf("\n");
		A = (int*)malloc(sizeof(int)*(n*n));
		B = (int*)malloc(sizeof(int)*(n*n));
		C = (int*)malloc(sizeof(int)*(n*n));
		M = (int*)malloc(sizeof(int)*(n*n));
		P = (int*)malloc(sizeof(int)*(n));
		LlenarMatriz(A, n*n, 5);
		LlenarMatriz(B, n*n, 2);
	}	

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
	/*
		imprimir todo en batch si la cantidad de procesos es < 100
		generar un archivo para cada matriz en caso contrario.
	*/

           
    MPI_Finalize();
	//Limpiamos memoria
	free(A);
	free(B);
	free(C);
	free(M);
	free(P);
    return 0;
}