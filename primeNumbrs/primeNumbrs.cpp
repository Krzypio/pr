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
int metodaDzielacaRownolegla(bool*& arrayBoolean, int arraySize, int min, int numberOfThreads);
int metodaSitaSekwencyjna(bool*& arrayBoolean, int arraySize, int min, int max);
int metodaSitaRownoleglaDomenowa(bool*& arrayBoolean, int arraySize, int min, int max, int numberOfThreads);
int metodaSitaRownoleglaFunkcyjna(bool*& arrayBoolean, int arraySize, int min, int max, int numberOfThreads);

void printArrayBoolean(bool* arrayBoolean, int arraySize);
int* convertBoolToIntPrimeArray(bool* arrayBoolean, int arraySize, int primeCount, int min);
void printArrayInt(int* arrayInt, int arraySize);

void wypiszNaglowekPelny();
void wypiszWynik(const char methodName[], int lowerBound, int upperBound, int threadNumber, int primeCounter, clock_t methodClock, clock_t referenceClock);

clock_t zbadajMetodaDzielacaSekwencyjna(int min, int max);
void zbadajMetodaDzielacaRownolegla(int min, int max, int numberOfThreads, clock_t referenceClock);
clock_t zbadajMetodaSitaSekwencyjna(int min, int max);
void zbadajMetodaSitaRownoleglaDomenowa(int min, int max, int numberOfThreads, clock_t referenceClock);
void zbadajMetodaSitaRownoleglaFunkcyjna(int min, int max, int numberOfThreads, clock_t referenceClock);
void wykonajAnalize(int lowerBound, int upperBound);
//----------------------------

int main(int argc, char* argv[])
{
	clock_t clockCliks = clock();

	wykonajAnalize(2, 10'000'000);

	clockCliks = clock() - clockCliks;
	printf("Analize przeprowadzono w czasie %f [sec];\n", ((float)clockCliks) / CLOCKS_PER_SEC);

	/*int min = 0;
	int max = 100'000;
	int numberOfThreads = 8;
	 
	clock_t clockCliks = clock();

	const int arraySize = max - min + 1;
	bool* arrayBoolean = new bool[arraySize];

	//int primeArraySize = metodaDzielacaSekwencyjna(arrayBoolean, arraySize, min);
	//int primeArraySize = metodaDzielacaRownolegla(arrayBoolean, arraySize, min, numberOfThreads);
	//int primeArraySize = metodaSitaSekwencyjna(arrayBoolean, arraySize, min, max);
	//int primeArraySize = metodaSitaRownoleglaDomenowa(arrayBoolean, arraySize, min, max, numberOfThreads);
	int primeArraySize = metodaSitaRownoleglaFunkcyjna(arrayBoolean, arraySize, min, max, numberOfThreads);
	
	printf("\tIt took me %d clicks (%f seconds).\n", clockCliks,((float)clockCliks) / CLOCKS_PER_SEC);

	clockCliks = clock() - clockCliks;

	int* primeArray = convertBoolToIntPrimeArray(arrayBoolean, arraySize, primeArraySize, min);
	printf("\tCount of primary numbers in range <%d;%d>: %d\n", min, max, primeArraySize);
	printArrayInt(primeArray, primeArraySize);*/
}

void wykonajAnalize(int lowerBound, int upperBound){

	int numberOfThreads = 0;
	int min = 0;
	int max = 0;

	wypiszNaglowekPelny();

	//Wariant A
	min = lowerBound;
	max = upperBound;
	clock_t referenceClockDzielenie = zbadajMetodaDzielacaSekwencyjna(min, max);

	zbadajMetodaDzielacaRownolegla(min, max, numberOfThreads = 1, referenceClockDzielenie);
	zbadajMetodaDzielacaRownolegla(min, max, numberOfThreads = 2, referenceClockDzielenie);
	zbadajMetodaDzielacaRownolegla(min, max, numberOfThreads = 4, referenceClockDzielenie);
	zbadajMetodaDzielacaRownolegla(min, max, numberOfThreads = 8, referenceClockDzielenie);

	clock_t referenceClockSito = zbadajMetodaSitaSekwencyjna(min, max);
	zbadajMetodaSitaRownoleglaDomenowa(min, max, numberOfThreads = 1, referenceClockSito);
	zbadajMetodaSitaRownoleglaDomenowa(min, max, numberOfThreads = 2, referenceClockSito);
	zbadajMetodaSitaRownoleglaDomenowa(min, max, numberOfThreads = 4, referenceClockSito);
	zbadajMetodaSitaRownoleglaDomenowa(min, max, numberOfThreads = 8, referenceClockSito);

	zbadajMetodaSitaRownoleglaFunkcyjna(min, max, numberOfThreads = 1, referenceClockSito);
	zbadajMetodaSitaRownoleglaFunkcyjna(min, max, numberOfThreads = 2, referenceClockSito);
	zbadajMetodaSitaRownoleglaFunkcyjna(min, max, numberOfThreads = 4, referenceClockSito);
	zbadajMetodaSitaRownoleglaFunkcyjna(min, max, numberOfThreads = 8, referenceClockSito);
	
	//Wariant B
	min = lowerBound;
	max = upperBound/2;
	referenceClockDzielenie = zbadajMetodaDzielacaSekwencyjna(min, max);
	referenceClockSito = zbadajMetodaSitaSekwencyjna(min, max);

	zbadajMetodaDzielacaRownolegla(min, max, numberOfThreads = 1, referenceClockDzielenie);
	zbadajMetodaDzielacaRownolegla(min, max, numberOfThreads = 2, referenceClockDzielenie);
	zbadajMetodaDzielacaRownolegla(min, max, numberOfThreads = 4, referenceClockDzielenie);
	zbadajMetodaDzielacaRownolegla(min, max, numberOfThreads = 8, referenceClockDzielenie);

	zbadajMetodaSitaRownoleglaDomenowa(min, max, numberOfThreads = 1, referenceClockSito);
	zbadajMetodaSitaRownoleglaDomenowa(min, max, numberOfThreads = 2, referenceClockSito);
	zbadajMetodaSitaRownoleglaDomenowa(min, max, numberOfThreads = 4, referenceClockSito);
	zbadajMetodaSitaRownoleglaDomenowa(min, max, numberOfThreads = 8, referenceClockSito);

	zbadajMetodaSitaRownoleglaFunkcyjna(min, max, numberOfThreads = 1, referenceClockSito);
	zbadajMetodaSitaRownoleglaFunkcyjna(min, max, numberOfThreads = 2, referenceClockSito);
	zbadajMetodaSitaRownoleglaFunkcyjna(min, max, numberOfThreads = 4, referenceClockSito);
	zbadajMetodaSitaRownoleglaFunkcyjna(min, max, numberOfThreads = 8, referenceClockSito);

	//Wariant C
	min = upperBound/2;
	max = upperBound;
	referenceClockDzielenie = zbadajMetodaDzielacaSekwencyjna(min, max);
	referenceClockSito = zbadajMetodaSitaSekwencyjna(min, max);

	zbadajMetodaDzielacaRownolegla(min, max, numberOfThreads = 1, referenceClockDzielenie);
	zbadajMetodaDzielacaRownolegla(min, max, numberOfThreads = 2, referenceClockDzielenie);
	zbadajMetodaDzielacaRownolegla(min, max, numberOfThreads = 4, referenceClockDzielenie);
	zbadajMetodaDzielacaRownolegla(min, max, numberOfThreads = 8, referenceClockDzielenie);

	zbadajMetodaSitaRownoleglaDomenowa(min, max, numberOfThreads = 1, referenceClockSito);
	zbadajMetodaSitaRownoleglaDomenowa(min, max, numberOfThreads = 2, referenceClockSito);
	zbadajMetodaSitaRownoleglaDomenowa(min, max, numberOfThreads = 4, referenceClockSito);
	zbadajMetodaSitaRownoleglaDomenowa(min, max, numberOfThreads = 8, referenceClockSito);

	zbadajMetodaSitaRownoleglaFunkcyjna(min, max, numberOfThreads = 1, referenceClockSito);
	zbadajMetodaSitaRownoleglaFunkcyjna(min, max, numberOfThreads = 2, referenceClockSito);
	zbadajMetodaSitaRownoleglaFunkcyjna(min, max, numberOfThreads = 4, referenceClockSito);
	zbadajMetodaSitaRownoleglaFunkcyjna(min, max, numberOfThreads = 8, referenceClockSito);
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
int metodaDzielacaRownolegla(bool*& arrayBoolean, int arraySize, int min, int numberOfThreads) {
	int primeCounter = 0;
	omp_set_num_threads(numberOfThreads);
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
int metodaSitaRownoleglaDomenowa(bool*& arrayBoolean, int arraySize, int min, int max, int numberOfThreads) {
	initializeArrayBoolean(arrayBoolean, arraySize, min, max);
	int dividersSize = (int)sqrt(max) + 1;

	int nonPrimeCounter = 0;
	if (min <= 0)
		nonPrimeCounter++;
	if (min <= 1)
		nonPrimeCounter++;

	omp_set_num_threads(numberOfThreads);
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
int metodaSitaRownoleglaFunkcyjna(bool*& arrayBoolean, int arraySize, 
	int min, int max, int numberOfThreads) {
	initializeArrayBoolean(arrayBoolean, arraySize, min, max);
	int dividersSize = (int)sqrt(max) + 1;

	bool* dividers = new bool[dividersSize];
	computeDividers(dividers, dividersSize);

	int primeCounter = 0;

	// Filter non-prime from arrayBoolean
	omp_set_num_threads(numberOfThreads);
#pragma omp parallel
	{
#pragma omp for schedule(dynamic)
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

clock_t zbadajMetodaDzielacaSekwencyjna(int min, int max) {
	const int arraySize = max - min + 1;
	bool* arrayBoolean = new bool[arraySize];

	clock_t clockCliks = clock();
	int primeArraySize = metodaDzielacaSekwencyjna(arrayBoolean, arraySize, min);
	clockCliks = clock() - clockCliks;
	wypiszWynik("metodaDzielacaSekwencyjna", min, max, 1, primeArraySize, clockCliks, 0);
	return clockCliks;
}
void zbadajMetodaDzielacaRownolegla(int min, int max, int numberOfThreads, clock_t referenceClock) {
	const int arraySize = max - min + 1;
	bool* arrayBoolean = new bool[arraySize];

	clock_t clockCliks = clock();
	int primeArraySize = metodaDzielacaRownolegla(arrayBoolean, arraySize, min, numberOfThreads);
	clockCliks = clock() - clockCliks;
	wypiszWynik("metodaDzielacaRownolegla", min, max, numberOfThreads, primeArraySize, clockCliks, referenceClock);
}
clock_t zbadajMetodaSitaSekwencyjna(int min, int max) {
	const int arraySize = max - min + 1;
	bool* arrayBoolean = new bool[arraySize];

	clock_t clockCliks = clock();
	int primeArraySize = metodaSitaSekwencyjna(arrayBoolean, arraySize, min, max);
	clockCliks = clock() - clockCliks;
	wypiszWynik("metodaSitaSekwencyjna", min, max, 1, primeArraySize, clockCliks, 0);
	return clockCliks;
}
void zbadajMetodaSitaRownoleglaDomenowa(int min, int max, int numberOfThreads, clock_t referenceClock) {
	const int arraySize = max - min + 1;
	bool* arrayBoolean = new bool[arraySize];

	clock_t clockCliks = clock();
	int primeArraySize = metodaSitaRownoleglaDomenowa(arrayBoolean, arraySize, min, max, numberOfThreads);
	clockCliks = clock() - clockCliks;
	wypiszWynik("metodaSitaRownoleglaDomenowa", min, max, numberOfThreads, primeArraySize, clockCliks, referenceClock);
}
void zbadajMetodaSitaRownoleglaFunkcyjna(int min, int max, int numberOfThreads, clock_t referenceClock) {
	const int arraySize = max - min + 1;
	bool* arrayBoolean = new bool[arraySize];

	clock_t clockCliks = clock();
	int primeArraySize = metodaSitaRownoleglaFunkcyjna(arrayBoolean, arraySize, min, max, numberOfThreads);
	clockCliks = clock() - clockCliks;
	wypiszWynik("metodaSitaRownoleglaFunkcyjna", min, max, numberOfThreads, primeArraySize, clockCliks, referenceClock);
}

void wypiszNaglowekPelny() {
	printf("methodName;<lowerBound,upperBound>;numberOfThreads;primeCounter;");
	printf("czasPrzetwarzania[ms];przyspieszeniePR;predkoscPrztwarzania[szt/ms];efektywnoscPR;\n");
}

void wypiszWynik(const char methodName[], int lowerBound, int upperBound, int threadNumber, int primeCounter, clock_t methodClock, clock_t referenceClock) {
	long czasPrzetwarzania = (long)((methodClock *1000) / CLOCKS_PER_SEC);
	float przyspieszeniePR = (float)referenceClock/methodClock;
	long predkoscPrztwarzania = (upperBound - lowerBound + 1) / czasPrzetwarzania;
	float efektywnoscPR = przyspieszeniePR / threadNumber;
	printf("%s;<%d,%d>;%d;%d;", methodName, lowerBound, upperBound, threadNumber, primeCounter);
	printf("%d;%f;%d;%f;\n", czasPrzetwarzania, przyspieszeniePR, predkoscPrztwarzania, efektywnoscPR);
}