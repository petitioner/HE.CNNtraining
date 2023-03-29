
#include "Tools.h"
#include "Scheme.h"
#include <assert.h>

int CNNtraining(double **traindata, int *trainlabel, long trainfactorDim,
		long trainsampleDim, double **testdata, int *testlabel,
		long testfactorDim, long testsampleDim);


void MyBootStrap(Scheme &scheme, SecretKey &secretKey,
		Ciphertext &encWData, long labelnum, long slots, long batch,
		long fdimNums, long logp, long logq, long logQ, long logT = 3, long logI = 4);

