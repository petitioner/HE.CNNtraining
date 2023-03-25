
#include "Tools.h"

#include <NTL/ZZ.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <iomanip>     //setw
#include <string.h>
#include <string>
#include <vector>
#include <algorithm>    // std::shuffle + std::copy
#include <array>        // std::array
#include <random>       // std::default_random_engine
#include <chrono>       // std::chrono::system_clock

using namespace NTL;

long Tools::suggestLogN(long lambda, long logQ) {
	long NBnd = ceil(logQ * (lambda + 110) / 3.6);
	double logNBnd = log2((double) NBnd);
	return (long) ceil(logNBnd);
}

bool Tools::is_number(string &s) {
	return (strspn(s.c_str(), "-.0123456789") == s.size());
}

/**
 * To Extract X and Y from data.
 *      Step 2. Extract X and Y from data
 *              get X and Y as follows:
 *              X = | X11 X12 X13 ... X1d|
 *                  | X21 X22 X23 ... X2d|
 *                  |  .   .   .  ...  . |
 *                  |  .   .   .  ...  . |
 *                  | Xn1 Xn2 Xn3 ... Xnd|
 *              Y = [ Y1  Y2  Y3  ... Yn ]
 *
 * @param  : sampleDim : the number of rows in the data
 * @return : zData = (Y{-1,+1}, Y@X)
 * @author : no one
 */
void Tools::dataFromFile(string &path, long &factorDim, long &sampleDim,
		double **&X, int *&Y) {
	vector<vector<double>> zline;
	factorDim = 0; 	// dimension of (Y,X)
	sampleDim = 0;	// number of samples
	ifstream openFile(path.data());
	if (openFile.is_open()) {
		string line, temp;

		while (getline(openFile, line)) {
			if (factorDim == 0)
				for (long i = 0; i < line.length(); ++i)
					if (line[i] == ',')
						factorDim++;

			size_t pos = line.find_first_of(',', 0);
			string first_word = line.substr(0, pos);
			if (!is_number(first_word))
				getline(openFile, line);

			size_t start = 0, end = 0;
			vector<double> vecline;
			do {
				end = line.find_first_of(',', start);
				temp = line.substr(start, end);
				vecline.push_back(atof(temp.c_str()));
				start = end + 1;
			} while (start);
			zline.push_back(vecline);
			sampleDim++;
		}
	} else {
		cout << "Error: cannot read file" << endl;
	}

	X = new double*[sampleDim];
	Y = new int[sampleDim];

	for (long j = 0; j < sampleDim; ++j) {
		double *xj = new double[factorDim];
		Y[j] = zline[j][0];

		for (long i = 0; i < factorDim; ++i) {
			xj[i] = zline[j][i];
		}
		X[j] = xj;
	}

	cout << endl << "dataFromFile is done !" << endl;
}

/**
 * To Extract X and Y from data.
 *      Step 2. Extract X and Y from data
 *              get X and Y as follows:
 *              X = | 1 X11 X12 X13 ... X1d|
 *                  | 1 X21 X22 X23 ... X2d|
 *                  | .  .   .   .  ...  . |
 *                  | .  .   .   .  ...  . |
 *                  | 1 Xn1 Xn2 Xn3 ... Xnd|
 *              Y = [ Y1  Y2  Y3  ... Yn ]
 *
 * @param  : sampleDim : the number of rows in the data
 * @author : no one
 */
void Tools::zDataFromFile(string &path, long &factorDim, long &sampleDim,
		double **&X, int *&Y) {
	vector<vector<double>> zline;
	factorDim = 1; 	// dimension of (Y,X)
	sampleDim = 0;	// number of samples
	ifstream openFile(path.data());
	if (openFile.is_open()) {
		string line, temp;

		while (getline(openFile, line)) {
			if (factorDim == 1)
				for (long i = 0; i < line.length(); ++i)
					if (line[i] == ',')
						factorDim++;

			size_t pos = line.find_first_of(',', 0);
			string first_word = line.substr(0, pos);
			if (!is_number(first_word))
				getline(openFile, line);

			size_t start = 0, end = 0;
			vector<double> vecline;
			do {
				end = line.find_first_of(',', start);
				temp = line.substr(start, end);
				vecline.push_back(atof(temp.c_str()));
				start = end + 1;
			} while (start);
			zline.push_back(vecline);
			sampleDim++;
		}
	} else {
		cout << "Error: cannot read file" << endl;
	}

	cout << "factorDim:= " << factorDim << endl;
	cout << "sampleDim:= " << sampleDim << endl;
	cout << "zline.size()" << zline.size() << endl;
	cout << "zline[0].size()" << zline[0].size() << endl;

	X = new double*[sampleDim];
	Y = new int[sampleDim];

	for (long j = 0; j < sampleDim; ++j) {
		double *xj = new double[factorDim];
		Y[j] = zline[j][0];

		xj[0] = 1.0;
		for (long i = 1; i < factorDim; ++i) {
			xj[i] = zline[j][i];
		}
		X[j] = xj;   // zline[j][0] = 1.0; X[j] = zline[j];
	}

	cout << endl << "zDataFromFile is done !" << endl;
}

int Tools::classnum(int *Y, long sampleDim) {
	std::set<int> labelset(Y, Y + sampleDim);

	cout << endl;
	for (auto iter = labelset.begin(); iter != labelset.end(); ++iter)
		cout << *iter << "\t";
	cout << endl;

	return labelset.size();

}

void Tools::oneHotEncoding(int *&labels, long sampleDim, int **&label1hot) {
	label1hot = new int*[sampleDim];
	int labelnum = Tools::classnum(labels, sampleDim);
	cout << "Make Sure: labels start from 0 rather than 1" << endl;

	cout << endl << endl << "labelnum := " << labelnum << endl;

	for (long j = 0; j < sampleDim; ++j) {
		int *xj = new int[labelnum];
		for (long i = 0; i < labelnum; ++i) {
			if (i == labels[j])
				xj[i] = 1;
			else
				xj[i] = 0;
		}
		label1hot[j] = xj;
	}

	cout << endl << "oneHotEncoding is done !" << endl;
}

void Tools::printData(int **zData, long factorDim, long sampleDim) {
	for (long i = 0; i < sampleDim; ++i) {
		cout << "--------------------------" << endl;
		cout << "The " << std::setw(3) << i << "-th Row:";
		for (long j = 0; j < 8; ++j) {
			cout << std::showpos << std::fixed << std::setw(12) << zData[i][j];
		}
		cout << endl;
	}
	cout << "-------------------------------------------------------------"
			<< endl << endl;
}

void Tools::printData(double **zData, long factorDim, long sampleDim) {
	for (long i = 0; i < sampleDim; ++i) {
		cout << "--------------------------" << endl;
		cout << "The " << std::setw(3) << i << "-th Row:";
		for (long j = 0; j < 8; ++j) {
			cout << std::showpos << std::fixed << setprecision(8)
					<< std::setw(18) << zData[i][j];
		}
		cout << endl;
	}
	cout << "-------------------------------------------------------------"
			<< endl << endl;
}

void Tools::printData(complex<double> *zData, long factorDim, long sampleDim) {
	for (long i = 0; i < sampleDim; ++i) {
		cout << "--------------------------" << endl;
		cout << "The " << std::setw(3) << i << "-th Row:";
		for (long j = 0; j < factorDim; ++j) {
			cout << std::showpos << std::fixed << setprecision(8)
					<< std::setw(18) << (zData[i * factorDim + j].real());
		}
		cout << endl;
	}
	cout << "-------------------------------------------------------------"
			<< endl << endl;
}

double** Tools::zInvBFromTraindata(double **traindata, long &factorDim,
		long &sampleDim, double epsilon) {
	vector<vector<double>> zline;

	double **zData = new double*[sampleDim];
	;
	double **zInvB = new double*[sampleDim];
	;

	for (long j = 0; j < sampleDim; ++j) {
		double *zj = new double[factorDim];
		double *Bj = new double[factorDim];
		for (long i = 0; i < factorDim; ++i) {
			zj[i] = std::abs(traindata[j][i]);
			Bj[i] = std::abs(traindata[j][i]);
		}
		zData[j] = zj;
		zInvB[j] = Bj;
	}

	for (long j = 0; j < sampleDim; ++j) {
		for (long i = 1; i < factorDim; ++i) {
			zInvB[j][0] += zInvB[j][i];
		}
	}

	for (long j = 0; j < sampleDim; ++j) {
		for (long i = 1; i < factorDim; ++i) {
			zInvB[j][i] = zInvB[j][0];
		}
	}

	for (long j = 0; j < sampleDim; ++j) {
		for (long i = 0; i < factorDim; ++i) {
			zInvB[j][i] *= zData[j][i];
		}
	}

	for (long i = 0; i < factorDim; ++i) {
		for (long j = 1; j < sampleDim; ++j) {
			zInvB[0][i] += zInvB[j][i];
		}
	}

	for (long i = 0; i < factorDim; ++i) {
		for (long j = 1; j < sampleDim; ++j) {
			zInvB[j][i] = zInvB[0][i];
		}
	}

	for (long i = 0; i < factorDim; ++i) {
		for (long j = 0; j < sampleDim; ++j) {
			zInvB[j][i] = 1.0 / (epsilon + .25 * zInvB[j][i]);
		}
	}

	//cout << "after invert" << endl;
	//printData(zInvB, factorDim, sampleDim);

	return zInvB;

}

/**
 * Shuffle the order of each row in zData.
 * @param  : zData 
 * @param  : factorDim : dimension of the raw data [X:y]
 * @param  : sampleDim : the number of rows in the data, except for the first row which is a interpretation.
 * @param  : isfirst :
 * @return :
 * @author : no one
 */
void Tools::shuffleZData(double **zData, long factorDim, long sampleDim) {
	srand(time(NULL));
	double *tmp = new double[factorDim];
	for (long i = 0; i < sampleDim; ++i) {
		long idx = i + rand() / (RAND_MAX / (sampleDim - i) + 1);
		copy(zData[i], zData[i] + factorDim, tmp);
		copy(zData[idx], zData[idx] + factorDim, zData[i]);
		copy(tmp, tmp + factorDim, zData[idx]);
	}
}
/**
 * Shuffle the order of each row in X and each element in Y in the same time.
 *
 * @param  : factorDim : dimension of the raw data X
 * @param  : sampleDim : the number of rows in the data.
 * @param  : isfirst :
 * @return :
 *
 * @author : no one
 */
void Tools::shuffleDataSync(double **X, long factorDim, long sampleDim,
		double *Y) {
	srand(time(NULL));
	double *tmp = new double[factorDim];
	double temp;
	for (long i = 0; i < sampleDim; ++i) {
		long idx = i + rand() / (RAND_MAX / (sampleDim - i) + 1);
		copy(X[i], X[i] + factorDim, tmp);
		temp = Y[i];
		copy(X[idx], X[idx] + factorDim, X[i]);
		Y[i] = Y[idx];
		copy(tmp, tmp + factorDim, X[idx]);
		Y[idx] = temp;
	}
}

double Tools::trueIP(double *a, double *b, long size) {
	double res = 0.0;
	for (long i = 0; i < size; ++i) {
		res += a[i] * b[i];
	}
	return res;
}
double Tools::calculateAUC(double **zData, double *wData, long factorDim,
		long sampleDim, double &correctness, double &auc) {
//	cout << "w:";
	for (long i = 0; i < factorDim; ++i) {
//		cout << wData[i] << ",";
	}
//	cout << endl;

	long TN = 0, FP = 0;

	vector<double> thetaTN;
	vector<double> thetaFP;

	for (int i = 0; i < sampleDim; ++i) {
		if (zData[i][0] > 0) {
			if (Tools::trueIP(zData[i], wData, factorDim) < 0)
				TN++;
			thetaTN.push_back(
					zData[i][0]
							* Tools::trueIP(zData[i] + 1, wData + 1,
									factorDim - 1));
		} else {
			if (Tools::trueIP(zData[i], wData, factorDim) < 0)
				FP++;
			thetaFP.push_back(
					zData[i][0]
							* Tools::trueIP(zData[i] + 1, wData + 1,
									factorDim - 1));
		}
	}

	correctness = 100.0 - (100.0 * (FP + TN) / sampleDim);

	if (thetaFP.size() == 0 || thetaTN.size() == 0) {
		cout << "n_test_yi = 0 : cannot compute AUC" << endl;
		auc = 0.0;
	} else {
		auc = 0.0;
		for (long i = 0; i < thetaTN.size(); ++i) {
			for (long j = 0; j < thetaFP.size(); ++j) {
				if (thetaFP[j] <= thetaTN[i])
					auc++;
			}
		}
		auc /= thetaTN.size() * thetaFP.size();
		cout << "AUC: " << auc << endl;
	}

	return auc;
}

double Tools::calculateACC(double **zData, double *wData, long factorDim,
		long sampleDim, double &correctness, double &auc) {
//	cout << "w:";
	for (long i = 0; i < factorDim; ++i) {
//		cout << wData[i] << ",";
	}
//	cout << endl;

	long TN = 0, FP = 0;

	vector<double> thetaTN;
	vector<double> thetaFP;

	for (int i = 0; i < sampleDim; ++i) {
		if (zData[i][0] > 0) {
			if (Tools::trueIP(zData[i], wData, factorDim) < 0)
				TN++;
			thetaTN.push_back(
					zData[i][0]
							* Tools::trueIP(zData[i] + 1, wData + 1,
									factorDim - 1));
		} else {
			if (Tools::trueIP(zData[i], wData, factorDim) < 0)
				FP++;
			thetaFP.push_back(
					zData[i][0]
							* Tools::trueIP(zData[i] + 1, wData + 1,
									factorDim - 1));
		}
	}

	correctness = 100.0 - (100.0 * (FP + TN) / sampleDim);

	if (thetaFP.size() == 0 || thetaTN.size() == 0) {
		cout << "n_test_yi = 0 : cannot compute AUC" << endl;
		auc = 0.0;
	} else {
		auc = 0.0;
		for (long i = 0; i < thetaTN.size(); ++i) {
			for (long j = 0; j < thetaFP.size(); ++j) {
				if (thetaFP[j] <= thetaTN[i])
					auc++;
			}
		}
		auc /= thetaTN.size() * thetaFP.size();
		//cout << "AUC: " << auc << endl;
	}

	return correctness;
}

double Tools::calculateMLE(double **zData, double *wData, long factorDim,
		long sampleDim, double &correctness, double &auc) {

	double *XW = new double[sampleDim]();
	for (int i = 0; i < sampleDim; ++i) {
		double res = 0;
		for (int j = 0; j < factorDim; ++j)
			res += zData[i][j] * wData[j];
		XW[i] = res;
	}
	double mle = 0.0;
	for (int i = 0; i < sampleDim; ++i) {
		double logsigm = -log(1 + exp(-XW[i]));
		mle += logsigm;
	}

	return mle;
}

/**
 * Returns the peak (maximum so far) resident set size (physical
 * memory use) measured in bytes, or zero if the value cannot be
 * determined on this OS.
 */
size_t Tools::getPeakRSS() {
#if defined(_WIN32)
    /* Windows -------------------------------------------------- */
    PROCESS_MEMORY_COUNTERS info;
    GetProcessMemoryInfo( GetCurrentProcess( ), &info, sizeof(info) );
    return (size_t)info.PeakWorkingSetSize;

#elif (defined(_AIX) || defined(__TOS__AIX__)) || (defined(__sun__) || defined(__sun) || defined(sun) && (defined(__SVR4) || defined(__svr4__)))
    /* AIX and Solaris ------------------------------------------ */
    struct psinfo psinfo;
    int fd = -1;
    if ( (fd = open( "/proc/self/psinfo", O_RDONLY )) == -1 )
        return (size_t)0L;      /* Can't open? */
    if ( read( fd, &psinfo, sizeof(psinfo) ) != sizeof(psinfo) )
    {
        close( fd );
        return (size_t)0L;      /* Can't read? */
    }
    close( fd );
    return (size_t)(psinfo.pr_rssize * 1024L);

#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
	/* BSD, Linux, and OSX -------------------------------------- */
	struct rusage rusage;
	getrusage(RUSAGE_SELF, &rusage);
#if defined(__APPLE__) && defined(__MACH__)
    return (size_t)rusage.ru_maxrss;
#else
	return (size_t)(rusage.ru_maxrss * 1024L);
#endif

#else
    /* Unknown OS ----------------------------------------------- */
    return (size_t)0L;          /* Unsupported. */
#endif
}

/**
 * Returns the current resident set size (physical memory use) measured
 * in bytes, or zero if the value cannot be determined on this OS.
 */
size_t Tools::getCurrentRSS() {
#if defined(_WIN32)
    /* Windows -------------------------------------------------- */
    PROCESS_MEMORY_COUNTERS info;
    GetProcessMemoryInfo( GetCurrentProcess( ), &info, sizeof(info) );
    return (size_t)info.WorkingSetSize;

#elif defined(__APPLE__) && defined(__MACH__)
    /* OSX ------------------------------------------------------ */
    struct mach_task_basic_info info;
    mach_msg_type_number_t infoCount = MACH_TASK_BASIC_INFO_COUNT;
    if ( task_info( mach_task_self( ), MACH_TASK_BASIC_INFO,
        (task_info_t)&info, &infoCount ) != KERN_SUCCESS )
        return (size_t)0L;      /* Can't access? */
    return (size_t)info.resident_size;

#elif defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
	/* Linux ---------------------------------------------------- */
	long rss = 0L;
	FILE *fp = NULL;
	if ((fp = fopen("/proc/self/statm", "r")) == NULL)
		return (size_t) 0L; /* Can't open? */
	if (fscanf(fp, "%*s%ld", &rss) != 1) {
		fclose(fp);
		return (size_t) 0L; /* Can't read? */
	}
	fclose(fp);
	return (size_t) rss * (size_t) sysconf(_SC_PAGESIZE);

#else
    /* AIX, BSD, Solaris, and Unknown OS ------------------------ */
    return (size_t)0L;          /* Unsupported. */
#endif
}

bool Tools::apowerof2(double floatnumber) {
	double l = log(floatnumber) / log(2);

	if (l > int(l))
		return false;
	if (l < int(l))
		return false;
	return true;
}
