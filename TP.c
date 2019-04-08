#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


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

/*
	Evalua si el numero recibido por parametro es primo haciendo uso de la raiz cuadrada.
*/
bool PruebaDePrimalidad(int valorAEvaluar){
	bool rst = true;
	int valSQRT = (int)sqrt(valorAEvaluar)+1;
	for(int valorAProbar = 2; valorAProbar < valSQRT && rst; valorAProbar++){
		if(valorAEvaluar%valorAProbar == 0) rst = false;
	}
	if(valSQRT <=1) rst = false;
	return rst;
}

/*
	Imprime el arreglo recibido por parametro.
*/
void ImprimirArreglo(int arregloAImprimir[], int numElem, int cntFilas, FILE* output){
	int cntRecorrida = 0;
	for(int index  = 0; index < numElem; index++){
		fprintf(output, "%d ", arregloAImprimir[index]);
		cntRecorrida++;
		if(cntRecorrida == cntFilas){
			cntRecorrida = 0;
			fprintf(output,"\n");
		}
	}
	fprintf(output,"\n");	
}

/*
	Hace la multiplicacion de matrices y va contando cuantos primos tiene M.
*/
int MultMatriz(int parteA[], int B[], int numElem, int numFilas, int parteM[]){
	int posM = 0, posA = 0, posB = 0, carry = 0, cntPrimos = 0; //valores
	for(int i = 0; i < numFilas; i++){//Me muevo por la fila de M
		for(int j = 0; j < numElem; j++){ //Me muevo por la columna de M con n y B con j;
			posM = (i*numElem)+j;
			carry = 0;
			for(int k = 0; k < numElem; k++){
				posA = (i*numElem)+k;
				posB = (k*numElem)+j;
				carry += parteA[posA]*B[posB];
			}
			if(PruebaDePrimalidad(carry)){ 
				cntPrimos++;
			};
			parteM[posM] = carry;
		}
	}
	return cntPrimos;
}

/*
	- Vector P tal que P[i] = a la cantidad de primos en la columna i de M.
	- Crear matriz C tal que:
		C[i,j] = M[i,j] + M[i+1,j] + M[i-1,j] + M[i,j+1] + M[i,j-1]
			(Considerar filas y columnas extremas de la matriz.)
*/
void CalcularCyP(int parteSuperior[], int parteInferior[], int parteC[], int parteM[], int myId, int localP[], int numElem, int cntFilas, int cntProcs){
	bool usarParteSuperior = true, usarParteInferior = false;
	int indexC = 0, indexCol;
	for(int i = 0; i < cntFilas; i++){
		indexCol = 0;
		if(i == cntFilas -1) usarParteInferior = true;//Determina cuando se debe usar la parte inferior.
		for(int j = 0; j < numElem; j++){

			indexC = (i*numElem) + j;//Indice de C que esta siendo calculado.

			parteC[indexC] += parteM[(i*numElem) + j];//C[i][j] += M[i][j].
			
			if(myId != cntProcs - 1){//C[i][j] += M[i+1][j]. Para todos menos para el hilo que maneja el ultimo pedazo de M.
				if(!usarParteInferior){
					parteC[indexC] += parteM[((i+1)*numElem) + j];
				} else {
					parteC[indexC] += parteInferior[indexCol];
				}
			}

			if(myId != 0){//C[i][j] += M[i-1][j]. Para todos menos para el hilo root.
				if(!usarParteSuperior){
					parteC[indexC] += parteM[((i-1)*numElem + j)];
				} else {
					parteC[indexC] += parteSuperior[indexCol];
				}
			}


			if( indexCol == 0 || indexCol != cntFilas-1){//C[i][j] += M[i][j+1]. No se calcula cuando j es el extremo derecho.
				parteC[indexC] += parteM[(i*numElem)+(j+1)];
			}

			if(indexCol != 0){//C[i][j] += M[i][j-1]. No se calcula cuando j es el extremo izquierdo.
				parteC[indexC] += parteM[(i*numElem) + (j-1)];
			}

			if(PruebaDePrimalidad(parteC[indexC])){//Aumenta el contador de primos en la columna j si el valor de C es primo.
				localP[j]++;
			}
			indexCol++;
		}
		usarParteSuperior = false;//ya se uso la parte superior
	}
}


int main(int argc,char **argv)
{
    int n , myid, numprocs, i, root = 0, tp = 0, cantFilasPorProc;
    double startwtime, endwtime;
	int *A, *B, *M, *C, *P;//Todas las matrices que se van a manejar.
	int *localA, *localM;
	FILE* output = stdout;
	//P contiene la cantidad de primos por fila.
    char processor_name[MPI_MAX_PROCESSOR_NAME];
	
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);
    //MPI_Get_processor_name(processor_name,&namelen);
    MPI_Barrier(MPI_COMM_WORLD);
	
	MPI_Status estado;
	//Se solicita al usuario que ingrese la cantidad de filas y columnas que van a tener las matrices.
	if(myid == root){
		printf("Digite el valor para las dimensiones de las matrices, debe ser multiplo de %d: ", numprocs);
		scanf("%d",&n);
		printf("\n");
		while(n%numprocs != 0){
			printf("El valor ingresado no es multiplo de %d, por favor, vuelva a ingresar otro: ", numprocs);
			scanf("%d",&n);
			printf("\n");
		}
		//Reservamos memoria para todos los arreglos globales.
		A = (int*)calloc(sizeof(int),(n*n));
		B = (int*)calloc(sizeof(int),(n*n));
		C = (int*)calloc(sizeof(int),(n*n));
		M = (int*)calloc(sizeof(int),(n*n));
		P = (int*)calloc(sizeof(int),(n));

		//Se llenan las matrices
		LlenarMatriz(A, n*n, 5);
		LlenarMatriz(B, n*n, 2);
	}

	//Se envia el valor de n a todos los procesos.
	MPI_Bcast(&n, 1, MPI_INT, root, MPI_COMM_WORLD);
	if(myid != root)
		B = (int*)calloc(sizeof(int),(n*n));

	//Envio la matriz B a todos los procesos.
	MPI_Bcast(B,n*n, MPI_INT, root, MPI_COMM_WORLD);
	cantFilasPorProc = n/numprocs;//Cuantas filas recibe cada proceso.
	localA = (int*)calloc(n*cantFilasPorProc,sizeof(int));//se reserva espacio de la cantidad de filas que le toca a cada proc.
	localM = (int*)calloc(n*cantFilasPorProc,sizeof(int));//Inicializo la parte de M correspondiente a cada proceso.

	//Envio a cada hilo la cantidad de elementos de A que le corresponden.
	MPI_Scatter(A, n*cantFilasPorProc, MPI_INT, localA, n*cantFilasPorProc, MPI_INT, root, MPI_COMM_WORLD);
	
	//Calculamos la multiplicacion local de la matriz y cuantos primos tiene la misma de forma local.
	int localTP = MultMatriz(localA, B, n/*#elem de cada fila * cantidad de filas*/, cantFilasPorProc, localM);

	//Armamos la matriz M a partir de los resultados locales de los procesos.
	MPI_Gather(localM, n*cantFilasPorProc, MPI_INT, M, n*cantFilasPorProc, MPI_INT, root, MPI_COMM_WORLD);

	//Hacemos la recoleccion de los resultados de cada proceso con la cantidad de numeros primos.
	MPI_Reduce(&localTP,&tp,1,MPI_INT,MPI_SUM,root,MPI_COMM_WORLD);

	//Todos los hilos ya tienen su parte correspondiente de M, les falta la fila arriba y la fila abajo para poder hacer el calculo de C.
	int* parteSuperior = (int*)calloc(n,sizeof(int));
	int* parteInferior = (int*)calloc(n,sizeof(int));
	int* localC = (int*)calloc(n*cantFilasPorProc,sizeof(int));
	int* localP = (int*)calloc(n,sizeof(int));
	
	int desplazamientoParteSuperior = (n*cantFilasPorProc)-n;//Direccion en el arreglo donde empieza la ultima fila.
	/*
	n -> n+1 -> ultima fila
	n -> n-1 -> primer fila
	n <- n+1 -> primer fila
	n <- n-1 -> ultima fila
	*/
	if(myid != root) MPI_Send(localM, n, MPI_INT, myid - 1, 2, MPI_COMM_WORLD);//Envio primer fila si no soy cero.
	if(myid != numprocs-1) MPI_Send(localM+desplazamientoParteSuperior, n, MPI_INT, myid+1, 1, MPI_COMM_WORLD);//Envio ultima fila si no soy el ultimo
	if(myid != root) MPI_Recv(parteSuperior, n, MPI_INT, myid-1, 1, MPI_COMM_WORLD, &estado);//Recibo la fila de arriba si no soy root.
	if(myid != numprocs - 1) MPI_Recv(parteInferior, n, MPI_INT, myid+1, 2, MPI_COMM_WORLD, &estado);//Recibo la fila de abajo si no soy root.
	
	//Calculamos matriz C y P.
	CalcularCyP(parteSuperior, parteInferior, localC, localM, myid, localP, n, cantFilasPorProc, numprocs);

	MPI_Gather(localC, n*cantFilasPorProc, MPI_INT, C, n*cantFilasPorProc, MPI_INT, root, MPI_COMM_WORLD);

	MPI_Reduce(localP, P, n, MPI_INT, MPI_SUM, root, MPI_COMM_WORLD);
	//Aqui se hechan todas las impresiones.
	if(myid == root){
		printf("A: \n");
		ImprimirArreglo(A, n*n, n, output);
		printf("B: \n");
		ImprimirArreglo(B, n*n, n, output);
		printf("\n");
		printf("M: \n");
		ImprimirArreglo(M, n*n, n, output);
		printf("Total de primos en M: %d\n\n", tp);
		printf("C: \n");
		ImprimirArreglo(C, n*n, n, output);
		printf("P: \n");
		ImprimirArreglo(P, n, n, output);
	}

	if (myid == root){
        free(A);
		free(M);
		free(C);
		free(P);
	}
	/*
		imprimir todo en batch si la cantidad de procesos es < 100
		generar un archivo para cada matriz en caso contrario.
	*/
	//Limpiamos memoria
	free(B);
	free(localA);
	free(localM);
	free(parteSuperior);
	free(parteInferior);
	free(localP);
	free(localC);
    MPI_Finalize();
    return 0;
}