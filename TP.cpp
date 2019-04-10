#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <string>
using namespace std;

/**
 *	Tarea programada hecha por Adrian Jose Alvarez Rodriguez 
 *  	Carnet: B40340
 * 	Arquitectura de computadoras II semestre 2019.
 */



/*
	Llena la matriz recibida por vector con numeros aleatorios cuyo valor maximo es el recibido por parametro.
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
	if(valSQRT <= 1) rst = false;
	return rst;
}

/*
	Imprime el arreglo en el archuvo (puede ser consola o un archivo.)
*/
void ImprimirArreglo(int arregloAImprimir[], int numElem, int cntFilas, string nombreArchivo, ostream &archivo){
	archivo << "Matriz ";
	archivo << nombreArchivo;
	archivo << "\n";
	int cntRecorrida = 0;
	for(int index = 0; index < numElem; index++){
		archivo << arregloAImprimir[index];
		archivo << " ";
		cntRecorrida++;
		if(cntRecorrida == cntFilas){
			cntRecorrida = 0;
			archivo << "\n";
		}
	}
	archivo << "\n";
}

/*
	Imprime en el archivo (puede ser consola o un archivo) el arreglo P.
*/
void ImprimirArregloConPrimos(int arregloAImprimir[], int numElem, ostream &archivo){
	for(int index = 0; index < numElem; index++){
		archivo << "P[";
		archivo << index;
		archivo<< "] = ";
		archivo << arregloAImprimir[index];
		archivo << " \n";
	}
}

/*
	Hace la multiplicacion de matrices y va contando cuantos primos tiene M en cada columna.
	Retorna la cantidad de primos que fueron calculados para el proceso correspondiente.
*/
int MultMatriz(int parteA[], int B[], int numElem, int numFilas, int parteM[], int localP[]){
	int posM = 0, posA = 0, posB = 0, carry = 0, cntPrimos = 0; //valores
	for(int i = 0; i < numFilas; i++){
		for(int j = 0; j < numElem; j++){
			posM = (i*numElem)+j;//Calculo posision en M
			parteM[posM] = 0;
			for(int k = 0; k < numElem; k++){
				posA = (i*numElem)+k;//Calculo posicion en A
				posB = (k*numElem)+j;//Calculo posicion en B
				parteM[posM] += parteA[posA]*B[posB];
			}
			if(PruebaDePrimalidad(parteM[posM])){ // Si el valor calculado de M es primo, aumenta los contadores correspondientes.
				cntPrimos++;
				localP[j]++;
			};
		}
	}
	return cntPrimos;
}

/*
	- Calcular matriz C tal que:
		- C[i,j] = M[i,j] + M[i+1,j] + M[i-1,j] + M[i,j+1] + M[i,j-1]
		-parteSuperior[] y parteInferior[] son los arreglos que le faltan al proceso para poder realizar los calculos en
		C[i][j+-1].
		-parteC[] y parteM[] son los arreglos locales para realizar los calculos.
		-myId es el identificador del proceso actual que va a realizar los calculos para saber si debe usar los arreglos superior e inferior.
*/
void CalcularC(int parteSuperior[], int parteInferior[], int parteC[], int parteM[], int myId, int numElem, int cntFilas, int cntProcs){
	bool usarParteSuperior = true, usarParteInferior = false;
	int indexC, indexCol = 0;
	for(int i = 0; i < cntFilas; i++){
		if(i == cntFilas -1) usarParteInferior = true;//Determina cuando se debe usar la parte inferior.
		for(int j = 0; j < numElem; j++){
			indexCol = j%numElem;//Calculo el valor de la columna real actual.

			indexC = (i*numElem) + j;//Indice de C que esta siendo calculado.

			parteC[indexC] += parteM[(i*numElem) + j];//C[i][j] += M[i][j].
		
			//C[i][j] += M[i+1][j]. 
			if(usarParteInferior){//Si se requiere usar el arreglo con la fila faltante de abajo.
				if( myId != cntProcs - 1 ){//El ultimo proceso no ocupa hacer este calculo.
					parteC[indexC] += parteInferior[indexCol];
				}
			} else {//No estamos en la ultima fila de localM del proceso, por lo que no se ocupa el arreglo con la parte faltante.
				parteC[indexC] += parteM[((i+1)*numElem) + j];
			}

			//C[i][j] += M[i-1][j].
			if(usarParteSuperior){//Si se requiere usar el arreglo que tiene la parte superior faltante en el proc.
				if(myId != 0){//El proc raiz no utiliza esta fila.
					parteC[indexC] += parteSuperior[indexCol];
				}
			} else {//Ya no se ocupa hacer la fila con la parte faltante, todos los procs participan.
				parteC[indexC] += parteM[((i-1)*numElem + j)];
			}


			
			if(indexCol != numElem-1){//C[i][j] += M[i][j+1]. No se calcula cuando j es el extremo derecho.
				parteC[indexC] += parteM[(i*numElem)+(j+1)];
			}

			if(indexCol != 0){//C[i][j] += M[i][j-1]. No se calcula cuando j es el extremo izquierdo.
				parteC[indexC] += parteM[(i*numElem) + (j-1)];
			}
		}
		usarParteSuperior = false;//ya se uso la parte superior, por lo que se descarta su uso.
	}
}


int main(int argc,char **argv)
{
    int n , myid, numprocs, i, root = 0, tp = 0, cantFilasPorProc;
    double tiempoInicio, tiempoFin, tiempoIniCalculos, tiempoFinCalculos;
	int *A, *B, *M, *C, *P;//Todas las matrices que se van a manejar.
	int *localA, *localM;
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
		tiempoInicio = MPI_Wtime();
		printf("Digite el valor para las dimensiones de las matrices, debe ser multiplo de %d: ", numprocs);
		scanf("%d",&n);
		printf("\n");
		while(n%numprocs != 0){
			printf("El valor ingresado no es multiplo de %d, por favor, vuelva a ingresar otro: ", numprocs);
			scanf("%d",&n);
			printf("\n");
		}
		tiempoIniCalculos = MPI_Wtime();
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
	int* localP = (int*)calloc(n,sizeof(int));

	//Envio a cada hilo la cantidad de elementos de A que le corresponden.
	MPI_Scatter(A, n*cantFilasPorProc, MPI_INT, localA, n*cantFilasPorProc, MPI_INT, root, MPI_COMM_WORLD);
	
	//Calculamos la multiplicacion local de la matriz, arreglo P y cuantos primos tiene la misma de forma local.
	int localTP = MultMatriz(localA, B, n/*#elem de cada fila * cantidad de filas*/, cantFilasPorProc, localM, localP);

	//Armamos la matriz M a partir de los resultados locales de los procesos.
	MPI_Gather(localM, n*cantFilasPorProc, MPI_INT, M, n*cantFilasPorProc, MPI_INT, root, MPI_COMM_WORLD);

	//Hacemos la recoleccion de los resultados de cada proceso con la cantidad de numeros primos.
	MPI_Reduce(&localTP,&tp,1,MPI_INT,MPI_SUM,root,MPI_COMM_WORLD);

	//Se recogen y suman los resultados de cada proceso de la cuenta de primos en cada columna para saber el total.
	MPI_Reduce(localP, P, n, MPI_INT, MPI_SUM, root, MPI_COMM_WORLD);

	//Todos los hilos ya tienen su parte correspondiente de M, les falta la fila arriba y la fila abajo para poder hacer el calculo de C
	// en todas sus entradas.
	int* parteSuperior = (int*)calloc(n,sizeof(int));
	int* parteInferior = (int*)calloc(n,sizeof(int));

	int* localC = (int*)calloc(n*cantFilasPorProc,sizeof(int));//Cada proceso ocupa almacenar la parte de C que va a calcular.
	
	int desplazamientoParteSuperior = (n*cantFilasPorProc)-n;//Offset del arreglo local de M donde empieza la ultima fila.
	
	/*
		Reparto a todos los procesos la parte de M que les falta para calcular C, es decir, si soy el proceso n,
		al proceso n+1 le doy mi ultima fila y al proceso n-1 le doy mi primer fila.
	*/
	if(myid != root) MPI_Send(localM, n, MPI_INT, myid - 1, 2, MPI_COMM_WORLD);//Envio primer fila si no soy root.
	if(myid != numprocs-1) MPI_Send(localM+desplazamientoParteSuperior, n, MPI_INT, myid+1, 1, MPI_COMM_WORLD);//Envio ultima fila si no soy el proceso.
	if(myid != root) MPI_Recv(parteSuperior, n, MPI_INT, myid-1, 1, MPI_COMM_WORLD, &estado);//Recibo la fila de arriba si no soy root.
	if(myid != numprocs - 1) MPI_Recv(parteInferior, n, MPI_INT, myid+1, 2, MPI_COMM_WORLD, &estado);//Recibo la fila de abajo si no soy el ultimo proceso.
	
	//Calculamos matriz C.
	CalcularC(parteSuperior, parteInferior, localC, localM, myid, n, cantFilasPorProc, numprocs);

	//Se recogen los resultados de cada proceso para armar la matriz C en el proceso raiz.
	MPI_Gather(localC, n*cantFilasPorProc, MPI_INT, C, n*cantFilasPorProc, MPI_INT, root, MPI_COMM_WORLD);

	tiempoFinCalculos = MPI_Wtime();
	//Aqui se hechan todas las impresiones.
	if(myid == root){
		if(n <= 30 ){//Imprimo en consola
			ImprimirArreglo(A, n*n, n, "A", cout);
			ImprimirArreglo(B, n*n, n, "B", cout);
			ImprimirArreglo(M, n*n, n, "M", cout);
			ImprimirArreglo(C, n*n, n, "C", cout);
			ImprimirArregloConPrimos(P, n, cout);
		} else {//Imprimo cada matriz en el archivo correspondiente.
			ofstream archivoA ("A.txt");
			ImprimirArreglo(A, n*n, n, "A", archivoA);
			ofstream archivoB ("B.txt");
			ImprimirArreglo(B, n*n, n, "B", archivoB);
			ofstream archivoM ("M.txt");
			ImprimirArreglo(M, n*n, n, "M", archivoM);
			ofstream archivoC ("C.txt");
			ImprimirArreglo(C, n*n, n, "C", archivoC);
			ofstream archivoP ("P.txt");
			ImprimirArregloConPrimos(P, n, archivoP);
		}

		//Se imprimen los valores generales de los calculos y duracion de los mismos.
		printf("\n");
		printf("Valor de n: %d\n",n);
		printf("Cantidad de Procesos: %d\n",numprocs);
		printf("Total de primos en M: %d\n", tp);
		tiempoFin = MPI_Wtime();
		printf("Tiempo total que el programa tardo ejecutandose: %lf s\n", tiempoFin - tiempoInicio);
		printf("Tiempo que el programa duro realizando los calculos: %lf s\n", tiempoFinCalculos - tiempoIniCalculos);
		printf("\n");
		//Libero memoria de las matrices que solo existen en el proc raiz.
		free(A);
		free(M);
		free(C);
		free(P);
	}

	//Limpiamos memoria local de cada proceso.
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
