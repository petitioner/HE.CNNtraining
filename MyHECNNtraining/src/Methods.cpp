
#include "Tools.h"
#include "Methods.h"

#include "Ciphertext.h"
#include "NTL/ZZX.h"
#include "Scheme.h"
#include "TestScheme.h"
#include "SecretKey.h"
#include "TimeUtils.h"
#include <cmath>

#include "SerializationUtils.h"  // error: ‘SerializationUtils’ has not been declared

#include "Tools.h"
#include <EvaluatorUtils.h>
#include <NTL/BasicThreadPool.h>
#include <NTL/RR.h>
#include <NTL/ZZ.h>

#include <iomanip>
#include <set>
#include <map>
#include <typeinfo>
#include <vector>
#include <algorithm>    // std::shuffle
#include <vector>        // std::vector
#include <random>       // std::default_random_engine
#include <chrono>       // std::chrono::system_clock

#include <unistd.h>

/*
 int CNNinference(double** testdata, double* testlabel, long factorDim, long sampleDim, double **CNNdate, long cnnWeightsLen, long *cnnWeightsDims)
 {
 long fdimBits = (long)ceil(log2(factorDim)); //ceil(x) : Rounds x upward, returning the smallest integral value that is not less than x.
 long sdimBits = (long)ceil(log2(sampleDim)); //log2(x) : Returns the binary (base-2) logarithm of x.

 long wBits = 45;                             // Δ (delta)
 long pBits = 20;

 long logN = Tools::suggestLogN(80, logQ);    // it should be the Security Parameter λ
 long slots = 1 << (logN - 1);                // slots := 
 if ( slots > (1 << fdimBits) * (1 << sdimBits) ) slots = (1 << fdimBits) * (1 << sdimBits);
 long sBits = (long)ceil(log2(slots));    

 long batch = long( slots / (1 << fdimBits) );// batch is the Number of several sample dimensions.
 long bBits = (long)ceil(log2(batch));        

 long rnum = (long)ceil((double)sampleDim / batch); // To Divide the whole Test Data into Several Batches

 cout << endl << endl;
 cout << "factorDim = " << factorDim << endl << "sampleDim = " << sampleDim << endl;
 cout << "batch = " << batch << ", slots = " << slots << ", rnum = " << rnum << endl;
 cout << "logQ = " << logQ << ", logN = " << logN << ", sdimBits = " << sdimBits << ", fdimBits = " << fdimBits << endl;

 TimeUtils timeutils;
 timeutils.start("Scheme generating...");
 Ring ring;
 SecretKey secretKey(ring);
 Scheme scheme(secretKey, ring); 
 //std::set<int> leftset, rightset;
 //std::set<int>::iterator it;
 int pos[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,28,30,32,34,36,38,40,42,44,46,48,50,54,56,60,64,66,72,78,84,90,96,102,108,112,114,120,126,128,140,168,256,512,1024,2048,3072,4096,5120,6144,7168,8192,9216,10240,11264,12288,13312,14336,15360,16384,17408,18432,19456,20480,21504,22528,23552,24576,25600,26624,27648,28672,29696,30720,31744};
 for(long i=0; i < 90; ++i) {
 //leftset.insert(pos[i] );
 scheme.addLeftRotKey(secretKey, pos[i] );
 }
 int ppos[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,28,30,32,34,36,38,40,42,44,46,48,50,54,56,60,64,66,72,78,84,90,96,102,108,112,114,120,126,128,140,168,256,369,420,457,485,497,499,505,508,511,512,513,517,523,528,529,535,536,541,547,548,550,553,556,559,561,562,564,565,568,569,571,574,575,577,579,580,581,583,585,586,587,589,591,592,593,595,597,598,599,600,601,602,603,604,605,606,607,608,609,610,611,612,613,614,615,616,617,618,619,620,621,622,623,624,625,626,628,630,632,634,636,638,640,642,644,646,648,650,651,652,653,654,655,656,657,658,659,660,661,662,663,664,665,666,667,668,669,670,671,672,673,674,675,676,682,688,694,700,706,712,718,720,724,728,730,734,736,738,740,742,744,746,748,750,752,754,756,758,759,760,761,762,763,764,765,766,767,768,769,770,771,772,773,774,775,776,777,778,779,780,781,782,783,784,1024,2048,3072,4096,5120,6144,7168,8192,9216,10240,11264,12288,13312,14336,15360,16384,17408,18432,19456,20480,21504,22528,23552,24576,25600,26624,27648,28672,29696,30720,31744};
 for(long i=0; i < 253; ++i) {
 //rightset.insert(ppos[i] );
 scheme.addRightRotKey(secretKey, ppos[i] );
 }
 timeutils.stop("Scheme generation");
 
 ////////////////////////////////////////////////// Data owner //////////////////////////////////////////////////

 ////////////////////////////////////////////////// Data owner //////////////////////////////////////////////////

 //////////////////////////////////////////////////////////////////////////////////
 Ciphertext encFirstData;                                                      //
 complex<double>* pzData = new complex<double>[slots]();                     //
 for (long j = 0; j < batch; ++j) {                                          //
 long rlen = (1 << fdimBits);                                              //
 for (long l = 0; l < factorDim; ++l) {                                    //
 pzData[(batch * 0 + j) * rlen + l].real(testdata[batch * 0 + j][l]);    //
 pzData[(batch * 0 + j) * rlen + l].imag(0);                             //
 }                                                                         //
 }                                                                           //
 scheme.encrypt(encFirstData, pzData, slots, wBits, logQ);                   //
 delete[] pzData;           
 SerializationUtils::writeCiphertext(encFirstData, "DataOwnerSend_encFirstData.txt");                                                 //
 //
 VolleyRevolverEncoding VRE(encFirstData, slots, 28, 28, NULL);                //
 ///////////////////////////////////////////////// Cloud server /////////////////////////////////////////////////
 long this_rownum = VRE.rownum;
 long this_colnum = VRE.colnum;
 long this_datacolnum = VRE.datacolnum;
 vector<int> indexsvec;
 for (long i=0; i < 51; i+=2) 
 indexsvec.push_back(i);
 indexsvec.push_back(1);
 indexsvec.push_back(56);

 std::map<long, Ciphertext> rotaterowincomplete_filterct;
 Ciphertext leftfilterct, rightfilterct;
 // Encoding Method: this->colnum * 10000 + pos*1000 + left(0)/right(1)
 for (auto it = indexsvec.begin() ; it != indexsvec.end(); ++it) {
 long idx = this_colnum * 10000 + *it * 1000 + 1;
 complex<double>* leftfilter = new complex<double> [slots]();
 for (int rowidx=0; rowidx < this_rownum; ++rowidx) {
 for (int idx=0; idx < this_colnum - *it; ++idx) {
 leftfilter[rowidx * this_datacolnum + idx] = 1;
 }
 }
 scheme.encrypt(leftfilterct, leftfilter, slots, 60, logQ);
 rotaterowincomplete_filterct[idx] = leftfilterct;

 idx = this_colnum * 10000 + *it * 1000 + 0;
 complex<double>* rightfilter = new complex<double> [slots]();
 for (int rowidx=0; rowidx < this_rownum; ++rowidx) {
 for (int idx=0; idx < *it; ++idx) {
 rightfilter[rowidx * this_datacolnum + this_colnum-*it +idx] = 1;
 }
 }
 scheme.encrypt(rightfilterct, rightfilter, slots, 60, logQ);
 rotaterowincomplete_filterct[idx] = rightfilterct;
 }

 Ciphertext zeroct;
 auto zeros = new complex<double>[slots]();
 scheme.encrypt(zeroct, zeros, slots, 60, logQ);

 long this__width = 28;
 long this__height = 28;
 long kernel_size = 3;
 long this__colnum = 784;
 long this__datacolnum = 1024;
 Ciphertext* cleanct = new Ciphertext[this__height - kernel_size + 1];
 for(long rowth=0; rowth < this__height - kernel_size + 1; ++rowth) { // ASSUME THAT: strides=(1, 1)
 // design a new matrix ....
 // rowfilter = [0] * ( rowth * (self.width - kernel_size + 1) ) +  [1] * (self.width - kernel_size + 1)
 // messagefilter = rowfilter + [0] * ( self.datacolnum - len(rowfilter) )
 auto messagefilter = new double[this__datacolnum]();
 for(long idx=rowth * (this__width - kernel_size + 1); idx < (rowth + 1) * (this__width - kernel_size + 1); ++idx ) {
 messagefilter[idx] = 1;
 }
 auto temp = VRE.fillfullof(messagefilter, this__datacolnum);
 
 scheme.encrypt(cleanct[rowth], temp, slots, 60, logQ);
 }
 Ciphertext* filterct = new Ciphertext[kernel_size*kernel_size];
 for(long idx=0; idx < kernel_size*kernel_size; ++idx) {
 long x = idx / kernel_size;
 long y = idx % kernel_size;
 long shift_0 = x;
 long shift_1 = y;

 // Build a new designed matrix to filter out the garbage values, with the help of shift(,)
 double* filtermatrix = new double[this__colnum]();
 for (int h=0; h < this__height; ++h) {
 for (int w=0; w < this__width; ++w) {
 if ( (w - shift_0) % kernel_size == 0 && w + kernel_size <= this__width )
 if ( (h - shift_1) % kernel_size == 0 && h + kernel_size <= this__height )
 filtermatrix[h * this__width + w] = 1;
 }
 }
 auto filtermessage = VRE.fillfullof(filtermatrix, this__colnum);
 
 scheme.encrypt(filterct[idx], filtermessage, slots, 60, logQ);
 }
 ///////////
 double* filtermatrix = new double[slots]();
 for (long rowidx=0; rowidx < this_rownum; ++rowidx)
 filtermatrix[rowidx * this__datacolnum] = 1;
 Ciphertext vre_mul_four_filterct;
 scheme.encrypt(vre_mul_four_filterct, filtermatrix, slots, 60, logQ);

 Ciphertext* vre_mul_four_filterct0 = new Ciphertext[this_rownum];
 long offset = 0;
 for (long r=0; r < this_rownum; ++r) {
 filtermatrix = new double[slots]();
 for (long rowidx=0; rowidx < this_rownum; ++rowidx)
 filtermatrix[offset + rowidx * this__datacolnum + (rowidx + r) % this_rownum] = 1;
 //Ciphertext filterct;
 scheme.encrypt(vre_mul_four_filterct0[r], filtermatrix, slots, 60, logQ);
 }

 Ciphertext* vre_mul_four_filterct1 = new Ciphertext[this_rownum];
 offset = 32;
 for (long r=0; r < this_rownum; ++r) {
 filtermatrix = new double[slots]();
 for (long rowidx=0; rowidx < this_rownum; ++rowidx)
 filtermatrix[offset + rowidx * this__datacolnum + (rowidx + r) % this_rownum] = 1;
 //Ciphertext filterct;
 scheme.encrypt(vre_mul_four_filterct1[r], filtermatrix, slots, 60, logQ);
 }
 ///////////
 filtermatrix = new double[slots]();
 for (long rowidx=0; rowidx < this_rownum; ++rowidx)
 filtermatrix[rowidx * this__datacolnum] = 1;
 Ciphertext MyDense_filterct;
 scheme.encrypt(MyDense_filterct, filtermatrix, slots, 60, logQ);
 
 Ciphertext* MyDense_filterct0 = new Ciphertext[this_rownum];
 offset = 0;
 for (long r=0; r < this_rownum; ++r) {
 filtermatrix = new double[slots]();
 for (long rowidx=0; rowidx < this_rownum; ++rowidx)
 filtermatrix[offset + rowidx * this__datacolnum + (rowidx + r) % this_rownum] = 1;
 //Ciphertext filterct;
 scheme.encrypt(MyDense_filterct0[r], filtermatrix, slots, 60, logQ);
 }
 ///////////////////////////////////////////////// Cloud server ///////////////////////////////////////////////// 


 //////////////////////////////////////////////// Model provider ////////////////////////////////////////////////
 cout << "Model provider : begin : CurrentRSS (MB): " << ( Tools::getCurrentRSS() /1024.0/1024.0 ) << endl;
 cout << "Model provider : begin : PeakRSS    (MB): " << ( Tools::getPeakRSS() /1024.0/1024.0 )    << endl;
 cout << "----------------------------------------" << endl;
 double* filter0 = new double[9](); 
 double bias0 = CNNdate[1][0];
 for(int i=0; i < 36; i+=4) { filter0[i/4] = CNNdate[0][i]; }

 double* filter1 = new double[9](); 
 double bias1 = CNNdate[1][1];
 for(int i=1; i < 36; i+=4) { filter1[i/4] = CNNdate[0][i]; }

 double* filter2 = new double[9](); 
 double bias2 = CNNdate[1][2];
 for(int i=2; i < 36; i+=4) { filter2[i/4] = CNNdate[0][i]; }

 double* filter3 = new double[9](); 
 double bias3 = CNNdate[1][3];
 for(int i=3; i < 36; i+=4) { filter3[i/4] = CNNdate[0][i]; }
 
 auto encKernelimage0 = VRE.spanKernelimage(scheme, filter0, 9, 3);
 cout << "encKernelimage0[0].logp = " << encKernelimage0[0].logp << endl;
 for (int i=0; i < 3*3; ++i)
 SerializationUtils::writeCiphertext(encKernelimage0[i], "ModelProviderSend_encKernelimage0["+ std::to_string(i) +"].txt");

 auto encKernelimage1 = VRE.spanKernelimage(scheme, filter1, 9, 3);
 cout << "encKernelimage1[0].logp = " << encKernelimage1[0].logp << endl;
 for (int i=0; i < 3*3; ++i)
 SerializationUtils::writeCiphertext(encKernelimage1[i], "ModelProviderSend_encKernelimage1["+ std::to_string(i) +"].txt");

 auto encKernelimage2 = VRE.spanKernelimage(scheme, filter2, 9, 3);
 cout << "encKernelimage2[0].logp = " << encKernelimage2[0].logp << endl;
 for (int i=0; i < 3*3; ++i)
 SerializationUtils::writeCiphertext(encKernelimage2[i], "ModelProviderSend_encKernelimage2["+ std::to_string(i) +"].txt");

 auto encKernelimage3 = VRE.spanKernelimage(scheme, filter3, 9, 3);
 cout << "encKernelimage3[0].logp = " << encKernelimage3[0].logp << endl;
 for (int i=0; i < 3*3; ++i)
 SerializationUtils::writeCiphertext(encKernelimage3[i], "ModelProviderSend_encKernelimage3["+ std::to_string(i) +"].txt");

 Ciphertext* convbias = new Ciphertext[4];
 for (long n=0; n < 4; ++n) {
 auto biases_zeros = new complex<double>[slots]();
 for(long rowidx=0; rowidx < VRE.rownum; ++rowidx) {
 long width=VRE.width;// - kernel_size + 1;
 long height=VRE.height;// - kernel_size + 1;
 //cout << "width = " << width << endl;
 //cout << "height= " << height << endl;
 for(long i=0; i < width*height; ++i) {
 biases_zeros[rowidx * VRE.datacolnum + i].real( CNNdate[1][n] );
 biases_zeros[rowidx * VRE.datacolnum + i].imag(0);
 }

 }
 scheme.encrypt(convbias[n], biases_zeros, VRE.slots_number, 60, logQ);
 SerializationUtils::writeCiphertext(convbias[n], "ModelProviderSend_convbias["+ std::to_string(n) +"].txt");
 }
 double** wmatrix = new double*[2704](); 
 for(int i=0; i < 2704; ++i) {
 double* temp = new double[64](); 
 for(int j=0; j < 64; ++j) {
 temp[j] = CNNdate[6][i * 64 + j];
 }
 wmatrix[i] = temp;
 }
 double** wmatrixT = new double*[64](); 
 for(int j=0; j < 64; ++j) {
 double* temp = new double[2704]();
 for(int i=0; i < 2704; ++i) {
 temp[i] = wmatrix[i][j]; 
 }
 wmatrixT[j] = temp;
 }
 long datacolnum = 1024;
 long rownum = 32;
 long colnum = 676;
 Ciphertext* wCT0 = new Ciphertext[4];
 for (long n=0; n < 4; ++n) {
 complex<double>* dense1_0 = new complex<double>[slots]();
 for (long r = 0; r < rownum; ++r) {

 for (long i=n; i < 2704; i+=4) {
 dense1_0[r * datacolnum + i/4].real( wmatrixT[r][i] );
 dense1_0[r * datacolnum + i/4].imag(0);

 //if (i < 32) cout << wmatrixT[r][i] << " " ;
 }
 //cout << endl;
 }
 scheme.encrypt(wCT0[n], dense1_0, slots, wBits, logQ);
 SerializationUtils::writeCiphertext(wCT0[n], "ModelProviderSend_wCT0["+ std::to_string(n) +"].txt");
 }

 Ciphertext* wCT1 = new Ciphertext[4];
 for (long n=0; n < 4; ++n) {
 complex<double>* dense1_0 = new complex<double>[slots]();
 for (long r = 0; r < rownum; ++r) {

 for (long i=n; i < 2704; i+=4) {
 dense1_0[r * datacolnum + i/4].real( wmatrixT[rownum + r][i] );
 dense1_0[r * datacolnum + i/4].imag(0);
 }
 }
 scheme.encrypt(wCT1[n], dense1_0, slots, wBits, logQ);
 SerializationUtils::writeCiphertext(wCT1[n], "ModelProviderSend_wCT1["+ std::to_string(n) +"].txt");
 }

 Ciphertext biasCT0;
 complex<double>* bias_0 = new complex<double>[slots]();
 for (long r = 0; r < rownum; ++r) {

 for (long i=0; i < 32; ++i) {
 bias_0[r * datacolnum + i] = CNNdate[7][i];
 }
 }
 scheme.encrypt(biasCT0, bias_0, slots, wBits, logQ);
 SerializationUtils::writeCiphertext(biasCT0, "ModelProviderSend_biasCT0.txt");

 Ciphertext biasCT1;
 complex<double>* bias_1 = new complex<double>[slots]();
 for (long r = 0; r < rownum; ++r) {

 for (long i=0; i < 32; ++i) {
 bias_1[r * datacolnum + 32 + i] = CNNdate[7][32 + i];
 }
 }
 scheme.encrypt(biasCT1, bias_1, slots, wBits, logQ);
 SerializationUtils::writeCiphertext(biasCT1, "ModelProviderSend_biasCT1.txt");


 double** wmatrix1 = new double*[64](); 
 for(int i=0; i < 64; ++i) {
 double* temp = new double[10](); 
 for(int j=0; j < 10; ++j) {
 temp[j] = CNNdate[12][i * 10 + j];
 }
 wmatrix1[i] = temp;
 }

 double** wmatrixT1 = new double*[10](); 
 for(int j=0; j < 10; ++j) {
 double* temp = new double[64]();
 for(int i=0; i < 64; ++i) {
 temp[i] = wmatrix1[i][j]; 
 }
 wmatrixT1[j] = temp;
 }
 Ciphertext weightct;
 complex<double>* dense2 = new complex<double>[slots]();
 for (long r = 0; r < 10; ++r) {
 for (long i=0; i < 64; ++i) {
 dense2[r * datacolnum + i].real( wmatrixT1[r][i] );
 dense2[r * datacolnum + i].imag(0);
 }
 }
 scheme.encrypt(weightct, dense2, slots, wBits, logQ);
 SerializationUtils::writeCiphertext(weightct, "ModelProviderSend_weightct.txt");

 Ciphertext biasct;
 complex<double>* biasvec = new complex<double>[slots]();
 for (long r = 0; r < 10; ++r) { 
 for (long i=0; i < 64; ++i) {
 biasvec[r * datacolnum + i].real( CNNdate[13][r] );
 biasvec[r * datacolnum + i].imag(0);
 }
 }
 scheme.encrypt(biasct, biasvec, slots, wBits, logQ);
 SerializationUtils::writeCiphertext(biasct, "ModelProviderSend_biasct.txt");

 cout << "----------------------------------------" << endl;
 cout << "Model provider : end   : CurrentRSS (MB): " << ( Tools::getCurrentRSS() /1024.0/1024.0 ) << endl;
 cout << "Model provider : end   : PeakRSS    (MB): " << ( Tools::getPeakRSS() /1024.0/1024.0 )    << endl;
 //////////////////////////////////////////////// Model provider ////////////////////////////////////////////////

 long counterr = 0;
 long countall = 0;
 double clientsendMem = 0;
 long long  FourConv2DTime = 0;
 long long  FourMyReLUTime = 0;
 long long  VreMulFourTime = 0;
 long long  MyReLU_Time = 0;
 long long  MyDenseTime = 0;
 long long  EachIterationTime = 0;
 for (long n = 0; n < rnum; ++n) {
 cout << "----------------------------------------" << endl;
 cout << "START : " << n+1 << "-th ROUND" << endl;
 cout << "----------------------------------------" << endl;
 cout << "CurrentRSS (MB): " << ( Tools::getCurrentRSS() /1024.0/1024.0 ) << endl;
 cout << "PeakRSS    (MB): " << ( Tools::getPeakRSS() /1024.0/1024.0 )    << endl;
 
 complex<double> *pzData = new complex<double> [slots]();

 if (n < rnum - 1) {

 for (long r = 0; r < batch; ++r) {
 long rlen = (1 << fdimBits);
 for (long i = 0; i < factorDim; ++i) {
 pzData[r * datacolnum + i].real( testdata[batch * n + r][i] );
 pzData[r * datacolnum + i].imag(0);
 }
 }

 } else {

 long rest = sampleDim - batch * (rnum - 1);
 for (long j = 0; j < rest; ++j) {
 long rlen = (1 << fdimBits);
 for (long l = 0; l < factorDim; ++l) {
 pzData[j * rlen + l].real( testdata[batch * n + j][l] );
 pzData[j * rlen + l].imag(0);
 }
 }
 }
 auto currmem0 = Tools::getCurrentRSS() /1024.0/1024.0;
 Ciphertext encTestData;
 scheme.encrypt(encTestData, pzData, slots, wBits, logQ);
 auto currmem1 = Tools::getCurrentRSS() /1024.0/1024.0;
 clientsendMem += double(currmem1 - currmem0);

 delete[] pzData;

 VolleyRevolverEncoding curVRE(encTestData, slots, 28, 28, NULL);
 encTestData.kill();

 
 Ciphertext** encK = new Ciphertext*[4];
 encK[0] = encKernelimage0;
 encK[1] = encKernelimage1;
 encK[2] = encKernelimage2;
 encK[3] = encKernelimage3;
 timeutils.start("Conv2D for four Kernels");
 auto convre0 = curVRE.Conv2D(scheme, secretKey, encK[0], 3, convbias[0], zeroct, cleanct, filterct, rotaterowincomplete_filterct);//.printImages(secretKey, scheme, 1).print(secretKey, scheme, 1);
 auto convre1 = curVRE.Conv2D(scheme, secretKey, encK[1], 3, convbias[1], zeroct, cleanct, filterct, rotaterowincomplete_filterct);//.printImages(secretKey, scheme, 1).print(secretKey, scheme, 1);
 auto convre2 = curVRE.Conv2D(scheme, secretKey, encK[2], 3, convbias[2], zeroct, cleanct, filterct, rotaterowincomplete_filterct);//.printImages(secretKey, scheme, 1).print(secretKey, scheme, 1);
 auto convre3 = curVRE.Conv2D(scheme, secretKey, encK[3], 3, convbias[3], zeroct, cleanct, filterct, rotaterowincomplete_filterct);//.printImages(secretKey, scheme, 1).print(secretKey, scheme, 1);
 timeutils.stop("Conv2D for four Kernels");
 FourConv2DTime += timeutils.timeElapsed;
 delete[] encK;

 
 timeutils.start("MyReLU for the outputs of four Kernels");
 auto myrelu0 = convre0->MyReLU(scheme, CNNdate[2][0], CNNdate[3][0], CNNdate[4][0], CNNdate[5][0]);
 auto myrelu1 = convre1->MyReLU(scheme, CNNdate[2][0], CNNdate[3][0], CNNdate[4][0], CNNdate[5][0]);
 auto myrelu2 = convre2->MyReLU(scheme, CNNdate[2][0], CNNdate[3][0], CNNdate[4][0], CNNdate[5][0]);
 auto myrelu3 = convre3->MyReLU(scheme, CNNdate[2][0], CNNdate[3][0], CNNdate[4][0], CNNdate[5][0]);
 timeutils.stop("MyReLU for the outputs of four Kernels");
 FourMyReLUTime += timeutils.timeElapsed;
 delete convre0;
 delete convre1;
 delete convre2;
 delete convre3;
 

 VolleyRevolverEncoding** v = new VolleyRevolverEncoding*[4];
 v[0] = myrelu0;
 v[1] = myrelu1;
 v[2] = myrelu2;
 v[3] = myrelu3;

 timeutils.start("vre_mul_four for dense1");
 auto vremydense1_0 = myrelu0->vre_mul_four(scheme, v, wCT0, 4, biasCT0, zeroct, vre_mul_four_filterct, vre_mul_four_filterct0, 0 );//.printImages(secretKey, scheme, 1).printparams().print(secretKey, scheme, 1);
 auto vremydense1_1 = myrelu0->vre_mul_four(scheme, v, wCT1, 4, biasCT1, zeroct, vre_mul_four_filterct, vre_mul_four_filterct1, 32 );//.printImages(secretKey, scheme, 1).printparams().print(secretKey, scheme, 1);
 auto vremydense10 = vremydense1_0->add(scheme, *vremydense1_1);
 timeutils.stop("vre_mul_four for dense1");
 VreMulFourTime += timeutils.timeElapsed;

 delete myrelu0;
 delete myrelu1;
 delete myrelu2;
 delete myrelu3;
 
 delete vremydense1_0;              delete vremydense1_1;
 auto vremydense1 = vremydense10->changeImagesize(64);
 delete vremydense10;


 timeutils.start("myrelu3");
 auto vremyrelu3 = vremydense1->MyReLU(scheme, CNNdate[8][0], CNNdate[9][0], CNNdate[10][0], CNNdate[11][0]);//.printImages(secretKey, scheme, 1).printparams();
 timeutils.stop("myrelu3");
 MyReLU_Time += timeutils.timeElapsed;

 delete vremydense1;


 timeutils.start("MyDense");
 auto vre_output = vremyrelu3->MyDense(scheme, weightct, 0, biasct, MyDense_filterct, MyDense_filterct0, 64);
 timeutils.stop("MyDense");
 MyDenseTime += timeutils.timeElapsed;

 delete vremyrelu3;


 //vreoutput->printargmax(secretKey, scheme);
 auto vreoutput = vre_output->changeImagesize(10);
 auto res = vreoutput->returnargmax(secretKey, scheme);
 delete vre_output;
 delete vreoutput;
 for (long idx=0; idx < batch; ++idx) {
 if ( n == rnum - 1 && idx == sampleDim - batch * (rnum - 1) )
 break;
 cout << res[idx] << " ";

 }
 
 cout << endl << "-----------real------------" << endl;

 for (long idx=batch*n; idx < batch*n+batch; ++idx) {
 if ( n == rnum - 1 && idx == sampleDim )
 break;
 cout << int(testlabel[idx]) << " ";
 }
 cout << endl;

 for (long idx=batch*n; idx < batch*n+batch; ++idx) {
 if ( n == rnum - 1 && idx == sampleDim )
 break;
 ++countall;
 if ( long(testlabel[idx]) != res[idx - batch*n] ) {
 ++counterr;
 cout << "x" << " ";
 }
 else 
 cout << "-" << " ";

 }
 delete[] res;
 cout << endl << endl << endl;

 cout << "----- countall = " << countall 
 << ", counterr = " << counterr 
 << ", (countall - counterr) / countall = " << double(countall - counterr) / countall << " -----" << endl; 
 cout << "-------------------- Four Conv2D Cost (ms):" <<  FourConv2DTime/ (n+1) << "--------------------" << endl;
 cout << "-------------------- Four MyReLU Cost (ms):" <<  FourMyReLUTime/ (n+1) << "--------------------" << endl;
 cout << "-------------------- Two VreMulFour Time = " <<  VreMulFourTime/ (n+1) << "--------------------" << endl;
 cout << "-------------------- MyReLU_Time         = " <<  MyReLU_Time   / (n+1) << "--------------------" << endl;
 cout << "-------------------- MyDenseTime         = " <<  MyDenseTime   / (n+1) << "--------------------" << endl;
 EachIterationTime = FourConv2DTime + FourMyReLUTime + VreMulFourTime + MyReLU_Time + MyDenseTime; 
 cout << "-------------------- EachIterationTime   = " <<EachIterationTime/(n+1) << "--------------------" << endl;
 cout << endl << "-------------------- END OF THE " << n+1 << "-th BATCH --------------------" << endl;

 }
 cout << endl << "--------------- END ---------------" << endl;

 return 0;
 }
 */

int CNNtraining(double **traindata, int *trainlabel, long trainfactorDim,
        long trainsampleDim, double **testdata, int *testlabel,
        long testfactorDim, long testsampleDim) {

    int labelnum = Tools::classnum(trainlabel, trainsampleDim);
    if (labelnum != Tools::classnum(testlabel, testsampleDim)) {
        cout
                << " Tools::classnum(trainlabel, trainsampleDim) != Tools::classnum(testlabel, testsampleDim) "
                << endl;
        exit(-1);
    }

    cout << "# of class label : " << endl;
    cout << labelnum << endl;

    int **trainlabel1hot, **testlabel1hot;
    Tools::oneHotEncoding(trainlabel, trainsampleDim, trainlabel1hot);
    Tools::oneHotEncoding(testlabel, testsampleDim, testlabel1hot);
    Tools::printData(trainlabel1hot, labelnum, trainsampleDim);

    cout << "X := " << endl;
    cout << "trainfactorDim = " << trainfactorDim << endl;
    Tools::printData(traindata, trainfactorDim, 2);
    double **zInvB = Tools::zInvBFromTraindata(traindata, trainfactorDim,
            trainsampleDim);
    cout << "after invert " << endl;
    Tools::printData(zInvB, trainfactorDim, 2);

    long fdimBits = (long) ceil(log2(trainfactorDim)); //ceil(x) : Rounds x upward, returning the smallest integral value that is not less than x.
    long cdimBits = (long) ceil(log2(labelnum)); //log2(x) : Returns the binary (base-2) logarithm of x.

    cout << "fdimBits = " << (fdimBits) << endl;
    long fdimNums = (1 << fdimBits);
    cout << "2^fdimBits = fdimNums = " << fdimNums << endl;

    cout << "cdimBits = " << (cdimBits) << endl;
    long cdimNums = (1 << cdimBits);
    cout << "2^cdimBits = cdimNums = " << cdimNums << endl;

    long logN = 16;
    long logQ = 990; // 991.300840336 > logQ  to ensure 128-bit security level. Security Parameter λ

    long slots = 1 << (logN - 1);
    long sBits = (long) ceil(log2(slots));

    long batch = long(slots / (1 << fdimBits)); // batch is the Number of samples encrypted in a single CipherText.
    long bBits = (long) ceil(log2(batch));

    long rnum = (long) ceil((double) trainsampleDim / batch); // To Divide the whole Data into Several CiphertTexts

    cout << endl << endl;
    cout << "trainfactorDim = " << trainfactorDim << endl << "trainsampleDim = "
            << trainsampleDim << endl;
    cout << "batch = " << batch << ", slots = " << slots << ", rnum = " << rnum
            << endl;
    cout << "logQ = " << logQ << ", logN = " << logN << ", cdimBits = "
            << cdimBits << ", fdimBits = " << fdimBits << endl;

    long wBits = 45;                             // Δ (delta)
    long pBits = 20;

    TimeUtils timeutils;
    timeutils.start("Scheme generating...");
    Ring ring;
    SecretKey secretKey(ring);
    Scheme scheme(secretKey, ring);
    //scheme.addLeftRotKeys(secretKey);
    //scheme.addRightRotKeys(secretKey);

    std::set<int> leftset, rightset;
    std::set<int>::iterator it;
    int pos[] = { +0, +1, +2, 3, +4, 5, 6, 7, +8, 9, 10, +16, +32, +64, +128, +256, +512, +1024,
            +1536, +2048, +2560, +3072, +3584, +4096, +4608, +5120, +5632,
            +6144, +6656, +7168, +7680, +8192, +8704, +9216, +9728, +10240,
            +10752, +11264, +11776, +12288, +12800, +13312, +13824, +14336,
            +14848, +15360, +15872, +16384, +16896, +17408, +17920, +18432,
            +18944, +19456, +19968, +20480, +20992, +21504, +22016, +22528,
            +23040, +23552, +24064, +24576, +25088, +25600, +26112, +26624,
            +27136, +27648, +28160, +28672, +29184, +29696, +30208, +30720,
            +31232, +31744, +32256  };
    for (long i = 0; i < 79; ++i) {
        //leftset.insert(pos[i] );
        scheme.addLeftRotKey(secretKey, pos[i]);
    }
    scheme.addRightRotKeys(secretKey);
    int ppos[] = {10, 20, 40, 80, 160, 320 };
    for(long i=0; i < 6; ++i) {
      //rightset.insert(ppos[i] );
      scheme.addRightRotKey(secretKey, ppos[i] );
    }
    timeutils.stop("Scheme generation");

//  long logT=3, logI=4;
//  timeutils.start("Bootstrap Key generating");
//  long bootlogq = 30  +10;
//  long lognslots = (long)ceil(log2(batch));  //batch = factorDim / cnum;
//  //scheme.addBootKey(secretKey, logn, logq+logI);
//  scheme.addBootKey(secretKey, lognslots, bootlogq +logI);
//  timeutils.stop("Bootstrap Key generated");

    ////////////////////////////////////////////////// Data owner ////////////////////////////////////////////////// 
    Ciphertext *encTrainData = new Ciphertext[rnum];
    Ciphertext *encTrainLab1 = new Ciphertext[rnum]; // Train Label One-Hot Encoding
    Ciphertext encW8ghtData; // = new Ciphertext[labelnum];   other method to achieve this goal might not work
    Ciphertext encVelocData;  // = new Ciphertext[labelnum];   //velocity
    Ciphertext encZInvBData;

    Ciphertext encZerosData;  // encZerosData = Enc(0)

    // encTrainData
    for (long i = 0; i < rnum - 1; ++i) {
        complex<double> *pzData = new complex<double> [slots](); //();

        for (long j = 0; j < batch; ++j) {
            for (long l = 0; l < trainfactorDim; ++l) {
                pzData[fdimNums * j + l].real(traindata[batch * i + j][l]);
            }
        }
        scheme.encrypt(encTrainData[i], pzData, slots, wBits, logQ);

        delete[] pzData;
    }
    complex<double> *pzData = new complex<double> [slots](); //();

    for (long j = 0; j < trainsampleDim - batch * (rnum - 1); ++j) {
        for (long l = 0; l < trainfactorDim; ++l) {
            pzData[fdimNums * j + l].real(traindata[batch * (rnum - 1) + j][l]);
        }
    }
    scheme.encrypt(encTrainData[rnum - 1], pzData, slots, wBits, logQ);

    delete[] pzData;


    // encTrainLab1   :   encrypt trainlabel1hot
    for (long i = 0; i < rnum - 1; ++i) {
        complex<double> *pzData = new complex<double> [slots](); //();

        for (long j = 0; j < batch; ++j) {
            for (long l = 0; l < labelnum; ++l) {
                pzData[fdimNums * j + l].real(trainlabel1hot[batch * i + j][l]);
            }
        }
        scheme.encrypt(encTrainLab1[i], pzData, slots, wBits, logQ);

        delete[] pzData;
    }
    complex<double> *pzDatae = new complex<double> [slots](); //();

    for (long j = 0; j < trainsampleDim - batch * (rnum - 1); ++j) {
        for (long l = 0; l < labelnum; ++l) {
            pzDatae[fdimNums * j + l].real(
                    trainlabel1hot[batch * (rnum - 1) + j][l]);
        }
    }
    scheme.encrypt(encTrainLab1[rnum - 1], pzDatae, slots, wBits, logQ);

    delete[] pzDatae;


    // encW8ghtData
    complex<double> *pzDatad = new complex<double> [slots](); //();

    for (long j = 0; j < labelnum; ++j) {   // Ensure: labelnum <= batch
        for (long l = 0; l < trainfactorDim; ++l) {
            pzDatad[fdimNums * j + l].real( 0.022 );
        }
    }
    scheme.encrypt(encW8ghtData, pzDatad, slots, wBits, logQ);

    delete[] pzDatad;

    // encVelocData
    encVelocData.copy(encW8ghtData);
    encVelocData.n = slots;


    // encZInvBData   :   trainfactorDim, trainsampleDim
    complex<double> *pzDataf = new complex<double> [slots](); //();

    for (long j = 0; j < labelnum; ++j) {
        for (long l = 0; l < trainfactorDim; ++l) {
            pzDataf[fdimNums * j + l].real(zInvB[j][l]);
        }
    }
    scheme.encrypt(encZInvBData, pzDataf, slots, wBits, logQ);

    delete[] pzDataf;

    complex<double> *dcw = scheme.decrypt(secretKey, encZInvBData);
    Tools::printData(dcw, fdimNums, 12);


    ////////////////////////////////////////////////// Data owner //////////////////////////////////////////////////

    //////////////////////////////////////////////// Model provider ////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////// Cloud server /////////////////////////////////////////////////

    // encZerosData
    complex<double> *pzDatag = new complex<double> [slots](); //();
    scheme.encrypt(encZerosData, pzDatag, slots, wBits, logQ);
    delete[] pzDatag;

    //ZZ* dummy;    // Generate dummy polynomial which is encoding of some special vector //
    //complex<double>* pvals = new complex<double>[slots]();
    //for (long i = 0; i < slots; i += fdimNums) {
    //  pvals[i].real(1.0);
    //}
    //dummy = new ZZ[1 << logN];
    //scheme.ring.encode(dummy, pvals, slots, pBits);
    //delete[] pvals;




    for (long iter = 0; iter < 1; ++iter) {

        Ciphertext *encProbData = new Ciphertext[rnum];


        Ciphertext *rotatedct = new Ciphertext[batch];
        NTL_EXEC_RANGE(batch, first, last);
        for (long r=first; r < last; ++r) {
            scheme.leftRotateFast(rotatedct[r], encVelocData, r * fdimNums);
        }
        NTL_EXEC_RANGE_END;


        int lg = int( ceil( log(trainfactorDim) / log(2) ) ); 
        Ciphertext* tempct = new Ciphertext[ batch ];


        cout << "Enter Into Outside Loop" << endl;
        //z = np.dot(X, velocity.T)
        for (int cti = 0; cti < rnum; ++cti) {

            // sub-matrix[X] * velocity.T
            // encTrainData[cti] * encVelocData
            encProbData[cti].copy(encZerosData);
            encProbData[cti].n = slots;
 

            cout << "Enter Into Inside Loop" << endl << endl;

            cout << "calculating the " << cti << "-th ciphertext encrypting sub-matrix  *   velcity.T" << endl << endl;

            // dealing with each rotatedct[rowidx]
            //NTL_EXEC_RANGE(batch, first, last);
            long first = 0, last = batch;
            for (long rowidx=first; rowidx < last; ++rowidx) {

                Ciphertext rotatedctidx;
                rotatedctidx.copy( rotatedct[rowidx] );
                rotatedctidx.n = slots;

                if(rotatedctidx.logq > encTrainData[cti].logq)
                    scheme.modDownToAndEqual(rotatedctidx, encTrainData[cti].logq);
                //if(rotatedctidx.logq < encTrainData[cti].logq)
                //    scheme.modDownToAndEqual(encTrainData[cti], rotatedctidx.logq);
                assert(encTrainData[cti].logq == rotatedctidx.logq);

                scheme.multAndEqual(rotatedctidx, encTrainData[cti]);


                for (long i=0; i < lg; ++i) {
                    //Ciphertext tempct;
                    scheme.leftRotateFast(tempct[rowidx], rotatedctidx, int( pow(2, i) ) );
                    if(rotatedctidx.logp > tempct[rowidx].logp) 
                        scheme.reScaleByAndEqual(rotatedctidx,  rotatedctidx.logp-tempct[rowidx].logp);
                    if(rotatedctidx.logp < tempct[rowidx].logp) 
                        scheme.reScaleByAndEqual(tempct[rowidx], tempct[rowidx].logp-rotatedctidx.logp);
                    if(rotatedctidx.logq > tempct[rowidx].logq) 
                        scheme.modDownToAndEqual(rotatedctidx,  tempct[rowidx].logq);
                    if(rotatedctidx.logq < tempct[rowidx].logq) 
                        scheme.modDownToAndEqual(tempct[rowidx], rotatedctidx.logq);

                    scheme.addAndEqual(rotatedctidx, tempct[rowidx]);
                }


                double* filtermatrix = new double[slots]();
                for (long idx=0; idx < batch; ++idx)
                    filtermatrix[idx * fdimNums] = 1;
                Ciphertext filterct;
                scheme.encrypt(filterct, filtermatrix, slots, wBits, rotatedctidx.logq);

                if(rotatedctidx.logq > filterct.logq) 
                    scheme.modDownToAndEqual(rotatedctidx,  filterct.logq);
                if(rotatedctidx.logq < filterct.logq) 
                    scheme.modDownToAndEqual(filterct, rotatedctidx.logq);
                scheme.multAndEqual(rotatedctidx, filterct);



                for (long i=0; i < lg; ++i) {
                    //Ciphertext tempct;
                    scheme.rightRotateFast(tempct[rowidx], rotatedctidx, int( pow(2, i) ) );
                    if(rotatedctidx.logp > tempct[rowidx].logp) 
                        scheme.reScaleByAndEqual(rotatedctidx,  rotatedctidx.logp-tempct[rowidx].logp);
                    if(rotatedctidx.logp < tempct[rowidx].logp) 
                        scheme.reScaleByAndEqual(tempct[rowidx], tempct[rowidx].logp-rotatedctidx.logp);
                    if(rotatedctidx.logq > tempct[rowidx].logq) 
                        scheme.modDownToAndEqual(rotatedctidx,  tempct[rowidx].logq);
                    if(rotatedctidx.logq < tempct[rowidx].logq) 
                        scheme.modDownToAndEqual(tempct[rowidx], rotatedctidx.logq);

                    scheme.addAndEqual(rotatedctidx, tempct[rowidx]);
                }



                filtermatrix = new double[slots]();
                for (long idx=0; idx < batch; ++idx)
                    if ((idx + rowidx) % batch < labelnum)
                        filtermatrix[idx * fdimNums + (idx + rowidx) % batch] = 1;

                for (long idx = 0; idx < batch; ++idx) {
                    for (long col = labelnum; col < fdimNums; ++col)
                        filtermatrix[idx * fdimNums + col] = 0;
                }

                scheme.encrypt(filterct, filtermatrix, slots, wBits, rotatedctidx.logq);

                if(rotatedctidx.logq > filterct.logq) 
                    scheme.modDownToAndEqual(rotatedctidx,  filterct.logq);
                if(rotatedctidx.logq < filterct.logq) 
                    scheme.modDownToAndEqual(filterct, rotatedctidx.logq);
                scheme.multAndEqual(rotatedctidx, filterct);
                filterct.kill();   


                if(encProbData[cti].logp > rotatedctidx.logp) 
                    scheme.reScaleByAndEqual(encProbData[cti],  encProbData[cti].logp-rotatedctidx.logp);
                if(encProbData[cti].logp < rotatedctidx.logp) 
                    scheme.reScaleByAndEqual(rotatedctidx, rotatedctidx.logp-encProbData[cti].logp);
                if(encProbData[cti].logq > rotatedctidx.logq) 
                    scheme.modDownToAndEqual(encProbData[cti],  rotatedctidx.logq);
                if(encProbData[cti].logq < rotatedctidx.logq) 
                    scheme.modDownToAndEqual(rotatedctidx, encProbData[cti].logq);  

            //}
            //NTL_EXEC_RANGE_END;
// Bugs Might Stem From This Block !
            //for (long rowidx=0; rowidx < batch; ++rowidx) {
                scheme.addAndEqual(encProbData[cti], rotatedctidx);
                rotatedctidx.kill();
            }
            //delete[] rotatedct;

cout << "Xsub * Velcity" << endl;
dcw = scheme.decrypt(secretKey, encProbData[cti]);
Tools::printData(dcw, fdimNums, 1);


            ////# Sigmoid(x) ~ poly3 = 0.5 + 0.10679534503216294.*x + -0.00038503259805075.*x.^3; (lambda = 128)
            ////# Sigmoid(x) ~ poly3 = 0.5 + 0.10679534503216294.*x * (1 − 0.003605331.*x.^2); (lambda = 128)
            ////P = 1.0 / ( 1.0 + np.exp(-z) ) 
            //P = hlambda(z)
            //z == encProbData[cti]   >>  -hlambda(encProbData[cti])
            Ciphertext ctx, ctx2, ctres;
            ctx.copy(encProbData[cti]); ctx.n = slots;

            scheme.square(ctx2, ctx);
            scheme.reScaleByAndEqual(ctx2,  ctx.logp);  

            scheme.multByConstAndEqual(ctx,   -0.106795345, ctx.logp);

            scheme.multByConstAndEqual(ctx2,  -0.003605331, ctx2.logp);
            scheme.addConstAndEqual(ctx2, 1.0, ctx2.logp);

            if(ctx.logq > ctx2.logq) 
                scheme.modDownToAndEqual(ctx,  ctx2.logq);
            if(ctx.logq < ctx2.logq) 
                scheme.modDownToAndEqual(ctx2, ctx.logq);
            assert(ctx.logq == ctx2.logq);            
            scheme.multAndEqual(ctx, ctx2);
            scheme.reScaleByAndEqual(ctx,  ctx2.logp);

            //Problems Stem from Here !!! 
            scheme.addConstAndEqual(ctx, -0.5, ctx.logp);

            encProbData[cti].copy(ctx);
            encProbData[cti].n = slots;

cout << "-P" << endl;
dcw = scheme.decrypt(secretKey, encProbData[cti]);
Tools::printData(dcw, fdimNums, 1);


            //gradient = (Y - P).T.dot(X)
            // (Y - P) = -P + Y   >>   encProbData[cti]
            if(encProbData[cti].logp > encTrainLab1[cti].logp)
                scheme.reScaleByAndEqual(encProbData[cti],  encProbData[cti].logp-encTrainLab1[cti].logp);
            if(encProbData[cti].logp < encTrainLab1[cti].logp) 
                scheme.reScaleByAndEqual(encTrainLab1[cti], encTrainLab1[cti].logp-encProbData[cti].logp);
            if(encProbData[cti].logq > encTrainLab1[cti].logq) 
                scheme.modDownToAndEqual(encProbData[cti],  encTrainLab1[cti].logq);
            if(encProbData[cti].logq < encTrainLab1[cti].logq) 
                scheme.modDownToAndEqual(encTrainLab1[cti], encProbData[cti].logq); 
            assert(encProbData[cti].logq == encTrainLab1[cti].logq);
            assert(encProbData[cti].logp == encTrainLab1[cti].logp);
            scheme.addAndEqual(encProbData[cti], encTrainLab1[cti]);

cout << "encTrainLab1[cti].logp " << encTrainLab1[cti].logp << endl;
cout << "encTrainLab1[cti].logq " << encTrainLab1[cti].logq << endl;


cout << "encProbData[cti].logp " << encProbData[cti].logp << endl;
cout << "encProbData[cti].logq " << encProbData[cti].logq << endl;

cout << "Y - P" << endl << endl;
dcw = scheme.decrypt(secretKey, encProbData[cti]);
Tools::printData(dcw, fdimNums, 1);



            //gradient = (Y - P).T.dot(X)
            double* filtermatrix = new double[slots]();
            for (long idx = 0; idx < batch; ++idx) {
                for (long col = 0; col < labelnum; ++col)
                    filtermatrix[idx * fdimNums + col] = 1;
            }
            Ciphertext filterct;
            scheme.encrypt(filterct, filtermatrix, slots, wBits, encProbData[cti].logq);

            if(encProbData[cti].logq > filterct.logq) 
                scheme.modDownToAndEqual(encProbData[cti],  filterct.logq);
            if(encProbData[cti].logq < filterct.logq) 
                scheme.modDownToAndEqual(filterct, encProbData[cti].logq);
            scheme.multAndEqual(encProbData[cti], filterct);

            if(encProbData[cti].logp > filterct.logp) 
                scheme.reScaleByAndEqual(encProbData[cti],  encProbData[cti].logp-filterct.logp);
            if(encProbData[cti].logp < filterct.logp) 
                scheme.reScaleByAndEqual(filterct, filterct.logp-encProbData[cti].logp);
            if(encProbData[cti].logq > filterct.logq) 
                scheme.modDownToAndEqual(encProbData[cti],  filterct.logq);
            if(encProbData[cti].logq < filterct.logq) 
                scheme.modDownToAndEqual(filterct, encProbData[cti].logq); 

cout << "Y - P" << endl << endl;
dcw = scheme.decrypt(secretKey, encProbData[cti]);
Tools::printData(dcw, fdimNums, 1);
cout << "encProbData[cti].logp " << encProbData[cti].logp << endl;
cout << "encProbData[cti].logq " << encProbData[cti].logq << endl;
            

            // // * labelnum +[1]) is to leave [1] column for an incomplete column shifting
            // if(fdimNums - labelnum > (ceil(trainfactorDim / labelnum) +1) * labelnum +1) {
            //     int lglabel = int( ceil( log(trainfactorDim) / log(2) ) ); 
            //     int lglabel = int( floor( log(512 / 10.0) / log(2) ) ); 
            //     cout << endl << "lglabel: " << lglabel <<  "\t" << 10 * int( pow(2, lglabel) ) << endl;

            //     int* rolist = new int[1+lglabel]();
            //     for (int i = 0; i < 1+lglabel; ++i) {
            //         cout << 10 * int( pow(2, i) ) << endl;
            //         rolist[i] = 10 * int( pow(2, i) );
            //     }

            //     int featureNum = 501;
            //     cout << "A" << (ceil(featureNum / 10.0) +1) * 10 << endl;
            //     auto newcolnum =   (ceil(featureNum / 10.0) +1) * 10;
            //     for(int i = lglabel; i>=0; --i) {
            //         while(newcolnum>=rolist[i]) {
            //             cout << rolist[i] << "\t";
            //             newcolnum -= rolist[i];
            //         }

            //     }
            // }

            //trainfactorDim = 441; labelnum = 10; 160 = 10 * 2^4; 320 = 10 * 2^5; 
            assert(160 + 320 > trainfactorDim + 2 * labelnum);

            Ciphertext tempct, tempct160;
            for (long i=0; i < 5; ++i) {

                if(i == 4) {
                    tempct160.copy(encProbData[cti]);
                    tempct160.n = slots;
                }
 
                scheme.rightRotateFast(tempct, encProbData[cti], labelnum * int( pow(2, i) ) );
                if(encProbData[cti].logp > tempct.logp) 
                    scheme.reScaleByAndEqual(encProbData[cti],  encProbData[cti].logp-tempct.logp);
                if(encProbData[cti].logp < tempct.logp) 
                    scheme.reScaleByAndEqual(tempct, tempct.logp-encProbData[cti].logp);
                if(encProbData[cti].logq > tempct.logq) 
                    scheme.modDownToAndEqual(encProbData[cti],  tempct.logq);
                if(encProbData[cti].logq < tempct.logq) 
                    scheme.modDownToAndEqual(tempct, encProbData[cti].logq);

                scheme.addAndEqual(encProbData[cti], tempct);
            }

            scheme.rightRotateFast(tempct, tempct160, 320 );
            if(encProbData[cti].logp > tempct.logp) 
                scheme.reScaleByAndEqual(encProbData[cti],  encProbData[cti].logp-tempct.logp);
            if(encProbData[cti].logp < tempct.logp) 
                scheme.reScaleByAndEqual(tempct, tempct.logp-encProbData[cti].logp);
            if(encProbData[cti].logq > tempct.logq) 
                scheme.modDownToAndEqual(encProbData[cti],  tempct.logq);
            if(encProbData[cti].logq < tempct.logq) 
                scheme.modDownToAndEqual(tempct, encProbData[cti].logq);
            scheme.addAndEqual(encProbData[cti], tempct);

cout << "Padding[Y - P]" << endl << endl;
dcw = scheme.decrypt(secretKey, encProbData[cti]);
Tools::printData(dcw, fdimNums, 1);
cout << "encProbData[cti].logp " << encProbData[cti].logp << endl;
cout << "encProbData[cti].logq " << encProbData[cti].logq << endl; 


        }

        cout << "gradient = (Y - P).T.dot(X)" << endl;
        Ciphertext encGradtData;
        encGradtData.copy(encZerosData);
        encGradtData.n = slots;

        Ciphertext* encTemptData = new Ciphertext[rnum];
        for (int rotidx = 0; rotidx < labelnum; ++rotidx) {
            Ciphertext encAccumData;
            encAccumData.copy(encZerosData);
            encAccumData.n = slots;
            for (int cti = 0; cti < rnum; ++cti) {
                // encProbData rotation
                // !!addLeftRotationKey
                scheme.leftRotateFast(encTemptData[cti], encProbData[cti], rotidx );


                if(encTemptData[cti].logq > encTrainData[cti].logq) 
                    scheme.modDownToAndEqual(encTemptData[cti],  encTrainData[cti].logq);
                if(encTemptData[cti].logq < encTrainData[cti].logq) 
                    scheme.modDownToAndEqual(encTrainData[cti], encTemptData[cti].logq);
                scheme.multAndEqual(encTemptData[cti], encTrainData[cti]);
                scheme.reScaleByAndEqual(encTemptData[cti],  encTrainData[cti].logp);


                if(encTemptData[cti].logp > encAccumData.logp) 
                    scheme.reScaleByAndEqual(encTemptData[cti],  encTemptData[cti].logp-encAccumData.logp);
                if(encTemptData[cti].logp < encAccumData.logp) 
                    scheme.reScaleByAndEqual(encAccumData, encAccumData.logp-encTemptData[cti].logp);
                if(encTemptData[cti].logq > encAccumData.logq) 
                    scheme.modDownToAndEqual(encTemptData[cti],  encAccumData.logq);
                if(encTemptData[cti].logq < encAccumData.logq) 
                    scheme.modDownToAndEqual(encAccumData, encTemptData[cti].logq);
                scheme.addAndEqual(encAccumData, encTemptData[cti]);
            }


            Ciphertext encResutData;
            encResutData.copy(encZerosData);
            encResutData.n = slots;
            for (int cti = 0; cti < batch; ++cti) {

                Ciphertext encTemppData;
                scheme.leftRotateFast(encTemppData, encAccumData, fdimNums * cti );

                if(encResutData.logp > encTemppData.logp) 
                    scheme.reScaleByAndEqual(encResutData,  encResutData.logp-encTemppData.logp);
                if(encResutData.logp < encTemppData.logp) 
                    scheme.reScaleByAndEqual(encTemppData, encTemppData.logp-encResutData.logp);
                if(encResutData.logq > encTemppData.logq) 
                    scheme.modDownToAndEqual(encResutData,  encTemppData.logq);
                if(encResutData.logq < encTemppData.logq) 
                    scheme.modDownToAndEqual(encTemppData, encResutData.logq);

                scheme.addAndEqual(encResutData, encTemppData);

                encTemppData.kill();

            }

    

            double* filtermatrix = new double[slots]();
            for (long idx = 0; idx < labelnum; ++idx) {
                for (long col = 0; col < fdimNums; ++col)
                    if( (col+labelnum-idx   +     rotidx)%labelnum == 0 )
                    filtermatrix[idx * fdimNums + col] = 1;
            }

            for (long idx = 0; idx < labelnum; ++idx) {
                for (long col = trainfactorDim; col < fdimNums; ++col)
                    filtermatrix[idx * fdimNums + col] = 0;
            }

            Ciphertext filterct;
            scheme.encrypt(filterct, filtermatrix, slots, wBits, encResutData.logq);

            if(encResutData.logq > filterct.logq) 
                scheme.modDownToAndEqual(encResutData,  filterct.logq);
            if(encResutData.logq < filterct.logq) 
                scheme.modDownToAndEqual(filterct, encResutData.logq);
            scheme.multAndEqual(encResutData, filterct);
            scheme.reScaleByAndEqual(encResutData,  filterct.logp);


            cout << endl << "part of gradient" << endl;
            dcw = scheme.decrypt(secretKey, encResutData);
            Tools::printData(dcw, fdimNums, 2);
            cout << "encResutData.logp " << encResutData.logp << endl;
            cout << "encResutData.logq " << encResutData.logq << endl; 




            if(encResutData.logp > encGradtData.logp) 
                scheme.reScaleByAndEqual(encResutData,  encResutData.logp-encGradtData.logp);
            if(encResutData.logp < encGradtData.logp) 
                scheme.reScaleByAndEqual(encGradtData, encGradtData.logp-encResutData.logp);
            if(encResutData.logq > encGradtData.logq) 
                scheme.modDownToAndEqual(encResutData,  encGradtData.logq);
            if(encResutData.logq < encGradtData.logq) 
                scheme.modDownToAndEqual(encGradtData,  encResutData.logq);

            scheme.addAndEqual(encGradtData, encResutData);


        }
        

            cout << endl << "Gradient:" << endl;
            dcw = scheme.decrypt(secretKey, encGradtData);
            Tools::printData(dcw, fdimNums, 12);















        ////P = 1.0 / ( 1.0 + np.exp(-z) ) 
        //P = hlambda(z)

        //gradient = (Y - P).T.dot(X)

        //MG = np.multiply(invBrow, gradient) 

        //    encProbData[r].kill();
        //delete[] encProbData;

    }

    return 0;
}
