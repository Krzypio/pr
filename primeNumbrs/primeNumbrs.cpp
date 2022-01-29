#include <stdio.h>
#include <time.h>
#include <omp.h>
#include <cmath>

using namespace std;

int isPrime(int number);
int countPrimes(bool*& arrayBoolean, int arraySize);
void initializeArrayBoolean(bool*& arrayBoolean, int arraySize, int min, int max);
void computeDividers(bool*& dividers, int dividersSize);

int metodaDzielacaSekwencyjna(bool*& arrayBoolean, int arraySize, int min);
int metodaDzielacaRownolegla(bool*& arrayBoolean, int arraySize, int min);
int metodaSitaSekwencyjna(bool*& arrayBoolean, int arraySize, int min, int max);
int metodaSitaRownoleglaDomenowa(bool*& arrayBoolean, int arraySize, int min, int max);
int metodaSitaRownoleglaFunkcyjna(bool*& arrayBoolean, int arraySize, int min, int max);

void printArrayBoolean(bool* arrayBoolean, int arraySize);
int* convertBoolToIntPrimeArray(bool* arrayBoolean, int arraySize, int primeCount, int min);
void printArrayInt(int* arrayInt, int arraySize);

int main(int argc, char* argv[])
{
	const int min = 0;
	const int max = 10'000'000;
	const int arraySize = max - min + 1;
	bool* arrayBoolean = new bool[arraySize];

	//TIME COUNT
	clock_t duration = clock();

	//int primeArraySize = metodaDzielacaSekwencyjna(arrayBoolean, arraySize, min);
	//int primeArraySize = metodaDzielacaRownolegla(arrayBoolean, arraySize, min);
	//int primeArraySize = metodaSitaSekwencyjna(arrayBoolean, arraySize, min, max);
	//int primeArraySize = metodaSitaRownoleglaDomenowa(arrayBoolean, arraySize, min, max);
	int primeArraySize = metodaSitaRownoleglaFunkcyjna(arrayBoolean, arraySize, min, max);

	duration = clock() - duration;
	printf("It took me %d clicks (%f seconds).\n", duration, ((float)duration) / CLOCKS_PER_SEC);
	//TIME COUNT

	int* primeArray = convertBoolToIntPrimeArray(arrayBoolean, arraySize, primeArraySize, min);
	printf("Count of primary numbers in range <%d;%d>: %d\n", min, max, primeArraySize);
	//printArrayInt(primeArray, primeArraySize);
}

int isPrime(int number) {
	if (number <= 1)
		return false;
	for (int i = 2; i <= sqrt(number); i++)
		if (number % i == 0)
			return false;
	return true;
}
int countPrimes(bool*& arrayBoolean, int arraySize) {
	int primeCounter = 0;
#pragma omp parallel
	{
		int localPrimeCounter = 0;
#pragma omp for nowait
		for (int i = 0; i < arraySize; i++)
		{
			if (true == arrayBoolean[i]) {
				localPrimeCounter++;
			}//if
		}//for
#pragma omp atomic
		primeCounter += localPrimeCounter;
	}//pragma omp parallel

	return primeCounter;
}
void initializeArrayBoolean(bool*& arrayBoolean, int arraySize, int min, int max) {
#pragma omp for
	for (int i = 0; i < arraySize; i++) {
		arrayBoolean[i] = true;
	}//for i

	if (0 == min) {
		arrayBoolean[0] = false;
		if (max >= 1) {
			arrayBoolean[1] = false;
		}//if2
	}//if1
	if (1 == min) {
		arrayBoolean[0] = false;
	}//if
}
void computeDividers(bool*& dividers, int dividersSize) {
	for (int i = 0; i < dividersSize; i++) {
		dividers[i] = true;
	}

	if (dividersSize > 0)
		dividers[0] = false;
	if (dividersSize > 1)
		dividers[1] = false;

	// Filter non-prime from dividers
	for (int i = 2; i < dividersSize; i++) {
		for (int j = i * i; j < dividersSize; j = j + i) {
			if (true == dividers[j]) {
				dividers[j] = false;
			}//if
		}//for j
	}//for i
}

int metodaDzielacaSekwencyjna(bool*& arrayBoolean, int arraySize, int min) {
	int primeCounter = 0;
	for (int i = 0; i < arraySize; i++) {
		if (isPrime(min + i)) {
			arrayBoolean[i] = true;
			primeCounter++;
		}//if
	}//for
	return primeCounter;
}//metodaDzielacaSekwencyjna()
int metodaDzielacaRownolegla(bool*& arrayBoolean, int arraySize, int min) {
	int primeCounter = 0;
#pragma omp parallel
	{
		int localPrimeCounter = 0;
#pragma omp for schedule (dynamic, 100)
		for (int i = 0; i < arraySize; i++) {
			if (isPrime(min + i)) {
				arrayBoolean[i] = true;
				localPrimeCounter++;
			}//if
			else {
				arrayBoolean[i] = false;
			}//else
		}//for
#pragma omp atomic
		primeCounter += localPrimeCounter;
	}//pragma omp parallel
	return primeCounter;
}//metodaDzielacaRownolegla()
int metodaSitaSekwencyjna(bool*& arrayBoolean, int arraySize, int min, int max) {
	initializeArrayBoolean(arrayBoolean, arraySize, min, max);
	int dividersSize = (int)sqrt(max) + 1;

	bool* dividers = new bool[dividersSize];
	computeDividers(dividers, dividersSize);

	int nonPrimeCounter = 0;
	if (min <= 0)
		nonPrimeCounter++;
	if (min <= 1)
		nonPrimeCounter++;
		
	// Filter non-prime from arrayBoolean
	for (int i = 0; i < dividersSize; i++) {
		if (dividers[i]) {
			int j = i * floor((min - 1) / i);
			if (j < i) {
				j += i;
			}//j
			while (j + i <= max) {
				j += i;
				int index = j - min;
				if (arrayBoolean[index]) {
					arrayBoolean[index] = false;
					nonPrimeCounter++;

				}//if
			}//while
		}//if dividers
	}//for i
	return max - min + 1 - nonPrimeCounter;
}
int metodaSitaRownoleglaDomenowa(bool*& arrayBoolean, int arraySize, int min, int max) {
	initializeArrayBoolean(arrayBoolean, arraySize, min, max);
	int dividersSize = (int)sqrt(max) + 1;

	int nonPrimeCounter = 0;
	if (min <= 0)
		nonPrimeCounter++;
	if (min <= 1)
		nonPrimeCounter++;

#pragma omp parallel
	{
		int localNonPrimeCounter = 0;
		//initialize dividers
		bool* dividers = new bool[dividersSize];
		computeDividers(dividers, dividersSize);

		int threadMax = omp_get_max_threads();//
		int threadId = omp_get_thread_num();//
		int minIndex = min + threadId * arraySize / threadMax;//
		int maxIndex = min + (threadId + 1) * arraySize / threadMax - 1;//

		//printf("threadMax=%d threadId=%d minIndex=%d maxIndex=%d \n", threadMax, threadId, minIndex, maxIndex);

		// Filter non-prime from arrayBoolean
		for (int i = 2; i < dividersSize; i++) {
			if (dividers[i]) {
				int j = i * floor((minIndex - 1) / i);//
				if (j < i) {
					j += i;
				}//j
				while (j + i <= maxIndex) {//
					j += i;
					int index = j - min;
					if (arrayBoolean[index]) {
						arrayBoolean[index] = false;
						localNonPrimeCounter++;
					}//if
				}//while
			}//if dividers
		}//for i
#pragma omp atomic
		nonPrimeCounter += localNonPrimeCounter;
	}//pragma omp parallel

	return max - min + 1 - nonPrimeCounter;
}
int metodaSitaRownoleglaFunkcyjna(bool*& arrayBoolean, int arraySize, int min, int max) {
	initializeArrayBoolean(arrayBoolean, arraySize, min, max);
	int dividersSize = (int)sqrt(max) + 1;

	bool* dividers = new bool[dividersSize];
	computeDividers(dividers, dividersSize);

	int primeCounter = 0;

	// Filter non-prime from arrayBoolean
#pragma omp parallel
	{
#pragma omp for
		for (int i = 0; i < dividersSize; i++) {
			if (dividers[i]) {
				int j = i * floor((min - 1) / i);
				if (j < i) {
					j += i;
				}//j
				while (j + i <= max) {
					j += i;
					int index = j - min;
					if (arrayBoolean[index]) {
						arrayBoolean[index] = false;
					}//if
				}//while
			}//if dividers
		}//for i

		//count primes 
		int localPrimeCounter = 0;
	#pragma omp for nowait
		for (int i = 0; i < arraySize; i++){
			if (true == arrayBoolean[i]) {
				localPrimeCounter++;
			}//if
		}//for
	#pragma omp atomic
		primeCounter += localPrimeCounter;
	}//pragma omp parallel

	return primeCounter;
}

void printArrayBoolean(bool* arrayBoolean, int arraySize)
{
	for (int i = 0; i < arraySize; i++) {
		//Next line after each 10 values
		if (i % 10 == 0)
			printf("\n");
		//Print values with ", " excluding endling for last value
		printf("%d", arrayBoolean[i]);
		if (i < arraySize - 1)
			printf(",\t");
	}
	printf("\n");
}
void printArrayInt(int* arrayInt, int arraySize)
{
	for (int i = 0; i < arraySize; i++) {
		//Next line after each 10 values
		if (i % 10 == 0)
			printf("\n");
		//Print values with ", " excluding endling for last value
		printf("%d", arrayInt[i]);
		if (i < arraySize - 1)
			printf(",\t");
	}
	printf("\n");
}
int* convertBoolToIntPrimeArray(bool* arrayBoolean, int arraySize, int primeCount, int min) {
	int* resultArray = new int[primeCount];
	int primeIndex = 0;
	for (int i = 0; i < arraySize; i++) {
		if (true == arrayBoolean[i]) {
			resultArray[primeIndex] = i + min;
			primeIndex++;
		}//if
	}//for
	return resultArray;
}