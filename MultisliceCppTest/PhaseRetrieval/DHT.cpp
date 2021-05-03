// DHT.cpp : This file contains the 'main' function. Program execution begins and ends there.

#include <complex.h>
#include <chrono>
#include <omp.h>
#include <filesystem>

#include "IXAHWave.h"
#include "XAHWave.h"
#include "XArray2D.h"
#include "XArray3D.h"
#include "XA_data.h"
#include "XA_file.h"
#include "XA_fft2.h"
#include "fftwd3drc.h"
#include "XA_spln2.h"
#include "XA_spln3.h"
#include "XA_iwfr.h"
#include "XA_TIE.h"
#include "XA_tiff.h"
#include "XA_move2.h"

using namespace xar;

#define NATOMMAX 500000 // maximum number of atomic positions that can be saved in an XYZ file

struct Pair2 { 	double v; int n; }; // pair of values for storing peak values and their initial indexes
int Pair2comp(const void* pP1, const void* pP2) 
{ 
	double a = ((const Pair2*)pP1)->v;
	double b = ((const Pair2*)pP2)->v;
	return (a > b) ? 1 : ((a < b) ? -1 : 0); 
}

int main(int argc, char* argv[])
{
	// start the execution timer
	std::chrono::system_clock::time_point start_time = std::chrono::system_clock::now();

	try
	{
		printf("\nStarting Differential Holographic Tomography program ...");
		constexpr double dAtomWidth = 1.0; // "average width" of the potential distribution around a single atom in A, this may change in the future
		vector<Pair> v2angles;
		vector<vector <Pair> > vvdefocus;
		vector<Pair> v2shifts;
		vector<double> vastigm;

		//************************************ read input parameters from file
		// read input parameter file
		char cline[1024], ctitle[1024], cparam[1024], cparam1[1024], cparam2[1024], cparam3[1024], cparam4[1024];
		
		string sInputParamFile("PhaseRetrieval.txt");
		if (argc > 1) sInputParamFile = argv[1];
		FILE* ff0 = fopen(sInputParamFile.c_str(), "rt");
		if (!ff0) throw std::exception(string("Error: cannot open parameter file " + sInputParamFile + ".").c_str());
		else printf("\nReading input parameter file %s ...", sInputParamFile.c_str());

		// read and skip an arbitrary number of initial comment lines (i.e. the lines that start with // symbols)
		while (true)
		{
			fgets(cline, 1024, ff0);
			if (!(cline[0] == '/' && cline[1] == '/')) break;
		}

		strtok(cline, "\n"); // 1. Verbose output during execution? Yes = 1, No = 0
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading verbose output parameter from input parameter file.");
		bool bVerboseOutput(true); // if this is TRUE, additional information is printed during execution
		(atoi(cparam) == 0 || atoi(cparam) == 1) ? bVerboseOutput = (bool)atoi(cparam) : throw std::exception("Error: verbose output parameter must be 0 or 1 in input parameter file.");

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 2. Input file with rotation angles and defocus distances
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading file name with rotation angles and defocus distances from input parameter file.");
		printf("\nReading defocus parameter file %s ...", cparam);
		bool bRelion(false);
		if (GetFileExtension(string(cparam)) == string(".TXT"))
			ReadDefocusParamsFile(cparam, v2angles, vvdefocus, bVerboseOutput);
		else
			if (GetFileExtension(string(cparam)) == string(".RELION"))
			{
				ReadRelionDefocusParamsFile(cparam, v2angles, vvdefocus, vastigm, v2shifts, bVerboseOutput);
				bRelion = true;
			}
			else throw std::exception("Error: unrecognised filename extension in parameter 1 of the input parameter file.");
		index_t nangles = v2angles.size(); // number of rotation steps 
		vector<index_t> vndefocus(nangles); // vector of numbers of defocus planes at different illumination angles
		for (index_t i = 0; i < nangles; i++) vndefocus[i] = vvdefocus[i].size();

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 3. Input filename base of defocus series of the sample in TIFF, GRD, GRC or RAW format
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading defocus series file name base from input parameter file.");
		string filenamebaseIn = cparam;
		bool bRAWinput;
		bool bTIFFinput;
		if (GetFileExtension(filenamebaseIn) == string(".TIFF") || GetFileExtension(filenamebaseIn) == string(".TIF"))
		{
			bTIFFinput = true; bRAWinput = false;
		}
		else if (GetFileExtension(filenamebaseIn) == string(".GRD") || GetFileExtension(filenamebaseIn) == string(".GRC"))
		{
			bTIFFinput = false; ; bRAWinput = false;
		}
		else if (GetFileExtension(filenamebaseIn) == string(".RAW"))
		{
			bTIFFinput = false; bRAWinput = true;
		}
		else throw std::exception("Error: input filename extension must be TIF, GRD, GRC or RAW.");
		if (GetFileExtension(filenamebaseIn) == string(".GRC")) // check that there is only one input defocused complex amplitude file per each illumination angle
		{
			for (index_t i = 0; i < nangles; i++)
				if (vndefocus[i] != 1)
					throw std::exception("Error: only one input defocused complex amplitude file per each illumination angle is allowed.");
		}

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 4. Width(pix) Height(pix) HeaderLength(bytes) Endian(0 = little, 1 = big) ElementLength(bytes) (for RAW files only)
		bool bBigEndian;
		index_t nx, ny; // transverse array sizes in pixels
		index_t nHeaderLength, nElementLength;
		if (bRAWinput)
		{
			if (sscanf(cline, "%s %s %s %s %s %s", ctitle, cparam, cparam1, cparam2, cparam3, cparam4) != 6)
				throw std::exception("Error reading Width Height HeaderLength Endianness ElementLength line from input parameter file.");
			nx = atoi(cparam);
			ny = atoi(cparam1);
			nHeaderLength = atoi(cparam2);
			bBigEndian = (bool)atoi(cparam3);
			nElementLength = atoi(cparam4);
			if (nElementLength != sizeof(float) && nElementLength != sizeof(double))
				throw std::exception("Unrecognised ElementLength in input parameter file (only sizeof(float) or sizeof(double) are allowed).");
			printf("\nRAW input file parameters: Width = %zd, Height = %zd, HeaderLength = %zd, Endianness = %d, ElementLength = %zd", nx, ny, nHeaderLength, bBigEndian, nElementLength);
		}

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 5. Xmin Xmax Ymin Ymax in Angstroms (for TIFF and RAW input files only)
		double xlo, xhi, xst; // physical x-bounds in angstroms
		double ylo, yhi, yst; // physical y-bounds in angstroms
		if (bTIFFinput || bRAWinput)
		{
			if (sscanf(cline, "%s %s %s %s %s", ctitle, cparam, cparam1, cparam2, cparam3) != 5) 
				throw std::exception("Error reading Xmin, Xmax, Ymin or Ymax from input parameter file.");
			xlo = atof(cparam);
			xhi = atof(cparam1);
			ylo = atof(cparam2);
			yhi = atof(cparam3);
			printf("\nPhysical boundaries of input images: Xmin = %g, Xmax = %g, Ymin = %g, Ymax = %g (Angstroms)", xlo, xhi, ylo, yhi);
		}
				
		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 6.Input data normalization factors f1 f2 (input -> (input / f1) + f2)
		double dNormFactor1(1.0), dNormFactor2(0.0);
		if (sscanf(cline, "%s %s %s", ctitle, cparam, cparam1) != 3)
			throw std::exception("Error reading input data normalization factors from input parameter file.");
		dNormFactor1 = atof(cparam);
		dNormFactor2 = atof(cparam1);
		printf("\nInput data normalization factors: f1 = %g, f2 = %g", dNormFactor1, dNormFactor2);
		if (dNormFactor1 == 0) throw std::exception("The first normalization factor cannot be zero (leads to division by zero).");

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 7. Wavelength in Angstroms
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading wavelength from input parameter file.");
		double wl = atof(cparam); // wavelength in Angstroms
		printf("\nWavelength = %g (A)", wl);
		if (wl < 0 || wl > 1)
			throw std::exception("Error: wavelength value appears to be wrong.");
		double EE; // incident electron energy in volts (recalculated from the wavelength below)
		constexpr double hp = 6.62607004e-34; // Planck's constant (m2 kg / s)
		constexpr double cc = 299792458; // speed of light (m / s)
		constexpr double ee = 1.602176634e-19; // electron charge (coulomb)
		constexpr double m0 = 9.1093837015e-31; // electron rest mass (kg)
		constexpr long double mc2 = m0 * cc * cc; // mc^2
		long double chl2 = long double(cc * cc * hp * hp) / long double(wl * wl * 1.e-20);
		long double abra = sqrt(mc2 * mc2 + chl2);
		EE = double((-mc2 + abra) / long double(ee));
		printf("\nIncident electron energy E = %g (keV)", EE / 1000.0);

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 8. Objective aperture in mrad
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading objective aperture from input parameter file.");
		double aobj = atof(cparam); 
		if (aobj != 0) printf("\nObjective aperture = %g (mrad)", aobj);
		else  printf("\nObjective aperture is infinite.");
		if (aobj < 0 || aobj > 1000)
			throw std::exception("Error: objective aperture value appears to be wrong.");
		double k2maxo = pow(aobj * 0.001f / wl, 2.0); // Fourier space bandwidth
		
		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 9. Spherical aberrations Cs3 and Cs5 in mm
		if (sscanf(cline, "%s %s %s", ctitle, cparam, cparam1) != 3) throw std::exception("Error reading spherical aberrations from input parameter file.");
		double Cs3 = atof(cparam);
		double Cs5 = atof(cparam1);
		printf("\nSpherical aberrations: Cs3 = %g, Cs5 = %g (mm)", Cs3, Cs5);
		Cs3 *= 1.e+7; // mm --> Angstroms
		Cs5 *= 1.e+7; // mm --> Angstroms

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 10.Phase retrieval method: 1 = IWFR, 2 = CTFL2, 3 = MinLogAmp, 4 = PhaseB7
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading phase retrieval method from input parameter file.");
		int nPhaseRetrieval = atoi(cparam);
		
		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 11. Save or not phase phase-retrieved defocused complex amplitudes in files
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading save_or_not phase-retrieved defocused complex amplitudes switch from input parameter file.");
		int nSaveDefocCAmps = atoi(cparam);
		if (nSaveDefocCAmps == 1)
			printf("\nPhase-retrieved defocused complex amplitudes will be saved in GRC files");
		else if (nSaveDefocCAmps == 0) 
				printf("\nPhase-retrieved defocused complex amplitudes will not be saved in files");
			else
				throw std::exception("Error: save_or_not phase-retrieved defocused complex amplitudes switch must be 0 or 1.");

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 12. Maximal number of IWFR iterations
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading maximal number of iterations from input parameter file.");
		int itermax = atoi(cparam);
		printf("\nMaximal number of iterations = %d", itermax);
		if (itermax < 1)
			throw std::exception("Error: the maximal number of iterations should be >= 1.");

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 13. Minimal phase reconstruction error(IWFR) or Tikhonov regularization parameter alpha(CTFL2)
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading minimal phase reconstruction error/regulariztion parameter from input parameter file.");
		double epsilon = atof(cparam);
		printf("\nMinimal phase reconstruciton error (IWFR) or Tikhonov regularization parameter (CTFL2) = %g", epsilon);
		if (epsilon < 0)
			throw std::exception("Error: minimal phase reconstruction error must be non-negative.");

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 14. Output defocus distances min max and step in Angstroms
		if (sscanf(cline, "%s %s %s %s", ctitle, cparam, cparam1, cparam2) != 4) throw std::exception("Error reading output defocus distances from input parameter file.");
		double zlo = atof(cparam); // minimum output defocus in Angstroms - !!! will be corrected with dzextra below
		double zhi = atof(cparam1); // maximum output defocus in Angstroms - !!! will be corrected with dzextra below 
		double zst = abs(atof(cparam2)); // output defocus step in Angstroms   
		if (zlo > zhi) std::swap(zlo, zhi);

		fgets(cline, 1024, ff0); strtok(cline, "\n"); //15. Extra defocus for 3D reconstruction in_Angstroms
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading extra defocus for 3D reconstruction parameter from input parameter file.");
		double dzextra = atof(cparam);
		printf("\nExtra defocus for 3D reconstruction = %g (Angstroms)", dzextra);
		double zlodz(zlo), zhidz(zhi);
		zlodz += dzextra; zhidz += dzextra;

		printf("\nOutput defocus distances: min = %g, max = %g, step = %g (Angstroms)", zlo, zhi, zst);
		if (zst == 0)
			throw std::exception("Error: output defocus step cannot be zero.");
		int noutdefocus = int((zhidz - zlodz) / zst + 0.5); // number of defocus planes to propagate to
		if (noutdefocus <= 0)
			throw std::exception("Error: number of output defocus planes must be positive.");
		vector<double> voutdefocus(noutdefocus); // vector of output defocus distances
		printf("\nThere are %d output defocus plane positions", noutdefocus);
		for (index_t n = 0; n < noutdefocus; n++) voutdefocus[n] = zlodz + zst * n;

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 16. 3D Laplacian filter mode
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading 3D Laplacian filter mode from input parameter file.");
		int imodeInvLaplace = atoi(cparam);
		switch (imodeInvLaplace)
		{
		case 0:
			printf("\n3D Laplacian filter won't be applied.");
			break;
		case 1:
			printf("\n3D Laplacian filter will be applied.");
			break;
		case 2:
			printf("\nThe program will apply the 3D Laplacian filter to the 3D potential imported from the input files.");
			break;
		default:
			throw std::exception("Error: unknown value for 3D Laplacian filter mode in input parameter file.");
		}

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 17. Regularization parameter for inverse 3D Laplacian
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading regularization parameter for inverse 3D Laplacian from input parameter file.");
		double alpha = atof(cparam);
		printf("\nRegularization parameter for inverse 3D Laplacian = %g", alpha);
		if (alpha < 0 && (imodeInvLaplace == 1 || imodeInvLaplace == 2))
			throw std::exception("Error: regularization parameter for 3D Laplacian filter must be non-negative.");

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 18. Low-pass filter width in Angstroms, background subtraction value and lower threshold level in Volts
		if (sscanf(cline, "%s %s %s %s", ctitle, cparam, cparam1, cparam2) != 4) throw std::exception("Error reading low-pass filter width, background subtraction value and lower threshold level from input parameter file.");
		double dlpfiltersize = atof(cparam);
		double dBackground = atof(cparam1);
		double dThreshold = atof(cparam2);
		printf("\nLow-pass filter width for 3D potential = %g (Angstroms)", dlpfiltersize);
		printf("\nBackground subtraction value for 3D potential = %g (Volts)", dBackground);
		printf("\nLower threshold level for 3D potential = %g (Volts)", dThreshold);
		
		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 19. Multislice reprojection of the 3D electrostatic potential mode
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading reprojectin of 3D potential mode from input parameter file.");
		int imodeReproj = atoi(cparam);
		switch (imodeReproj)
		{
		case 0:
			printf("\nMultislice reprojection of the 3D electrostatic potential won't be applied.");
			break;
		case 1:
			printf("\nMultislice reprojection of the 3D electrostatic potential will be applied.");
			break;
		case 2:
			printf("\nThe program will apply multislice reprojection of the 3D electrostatic potential imported from the input files.");
			break;
		default:
			throw std::exception("Error: unknown value for the multislice reprojection mode in input parameter file.");
		}

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 20. Slice thickness for multislice reprojection in Angstroms
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading slice thickness for multislice reprojection from input parameter file.");
		double sliceTh = atof(cparam);
		printf("\nSlice thickness for multislice reprojection = %g (Angstroms)", sliceTh);
		if (sliceTh < (4.0 * zst) && (imodeReproj == 1 || imodeReproj == 2))
			throw std::exception("Error: slice thickness for multislice reprojection must be 4 x z_step or larger.");

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 21. Peak localization mode
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading peak localization mode from input parameter file.");
		int imodePeaks = atoi(cparam);
		switch (imodePeaks)
		{
		case 0:
			printf("\nPeak localization in the 3D electrostatic potential won't be applied.");
			break;
		case 1:
			printf("\nPeak localization in the 3D electrostatic potential will be applied.");
			break;
		case 2:
			printf("\nThe program will apply peak localization in the 3D electrostatic potential imported from the input files.");
			break;
		default:
			throw std::exception("Error: unknown value for the peak localization mode in input parameter file.");
		}

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 22. Cube side length for peak localization (in_Angstroms)
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading cube side length for peak localization from input parameter file.");
		double datomsize = atof(cparam);
		printf("\nCube side length for peak localization = %g (Angstroms)", datomsize);
		if (int(datomsize / zst + 0.5) < 2 && (imodePeaks == 1 || imodePeaks == 2))
			throw std::exception("Error: cubic box side length for peak localization must be 2 x z_step or larger.");

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 23. Output file name base in GRD or TIFF format
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading output file name base from input parameter file.");
		string filenamebaseOut = cparam;
		printf("\nFile name base for 3D_potential = %s", filenamebaseOut.c_str());
		if (!std::filesystem::exists(GetPathFromFilename(filenamebaseOut)))
			throw std::exception("Error: the specified file folder for 3D potential does not seem to exist.");
		bool bTIFFoutput;
		if (GetFileExtension(filenamebaseOut) == string(".TIFF") || GetFileExtension(filenamebaseOut) == string(".TIF")) bTIFFoutput = true;
		else if (GetFileExtension(filenamebaseOut) == string(".GRD")) bTIFFoutput = false;
		else throw std::exception("Error: output filename extension must be TIF ot GRD.");

		if (imodeInvLaplace == 2 || imodeReproj == 2 || imodePeaks == 2) // only read in pre-caculated 3D potential from "output" files and filter or reproject it
		{
			if (imodeInvLaplace == 1 || imodeReproj == 1 || imodePeaks == 1)
				throw std::exception("Error: inconsistency between the Laplacian filter, reporojection and peak localization mode parameters.");
			if (GetFileExtension(filenamebaseOut) == string(".TIFF") || GetFileExtension(filenamebaseOut) == string(".TIF"))
			{
				bTIFFinput = true; bRAWinput = false;
			}
			else if (GetFileExtension(filenamebaseOut) == string(".GRD") || GetFileExtension(filenamebaseOut) == string(".GRC"))
			{
				bTIFFinput = false; ; bRAWinput = false;
			}
			else if (GetFileExtension(filenamebaseOut) == string(".RAW"))
			{
				bTIFFinput = false; bRAWinput = true;
			}
			else throw std::exception("Error: input filename extension (in this mode - it is taken from the output filename template) must be TIF, GRD, GRC or RAW.");
		}

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // // 24. Folder name for auxiliary files
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading folder name for auxiliary files from input parameter file.");
		string folderAux = cparam;
		printf("\nFolder for auxiliary file output = %s", folderAux.c_str());
		if (folderAux.rfind('\\') != folderAux.size() - 1) // the last character is not '\'
			if (folderAux.rfind('\\') == string::npos)
				throw std::exception("Error: auxiliary file folder name must contain at least one '\' symbol.");
			else folderAux.append("\\");
		if(!std::filesystem::exists(folderAux))
			throw std::exception("Error: the specified auxiliary file folder does not seem to exist.");

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 25. Number of parallel threads
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading number of parallel threads from input parameter file.");
		int nThreads = atoi(cparam);
		printf("\nNumber of parallel threads = %d", nThreads);
		if (nThreads < 1)
			throw std::exception("Error: the number of parallel threads in input parameter file should be >= 1.");
		omp_set_num_threads(nThreads);

		fclose(ff0); // close input parameter file

		//************************************ end reading input parameters from file

		//************************************ create vectors of input and output file names
		
		vector<string> vinfilenamesTot, vinfilenames1Tot; // input filenames for the defocused intensities or complex amplitudes
		vector< vector<string> > vvinfilenames(nangles), vvinfilenames1(nangles); // same input filenames for the defocused series in the form of vector of vectors
		if (imodeInvLaplace == 2 || imodeReproj == 2 || imodePeaks == 2) // only read in pre-caculated 3D potential from "output" files and filter or reproject it
		{
			printf("\nInput file name base for pre-existing 3D potential = %s", filenamebaseOut.c_str());
			FileNames(1, noutdefocus, filenamebaseOut, vinfilenamesTot); // create 1D array of input filenames to read the input 2D slices of a previously reconstructed 3D object
			vvinfilenames.resize(noutdefocus);
			for (index_t na = 0; na < noutdefocus; na++)
			{
				vvinfilenames[na].resize(1); // number of defocus planes at the current illumination angle
				vvinfilenames[na][0] = vinfilenamesTot[na];
			}
		}
		else // read in the input files as directed
		{
			printf("\nInput defocus series file name base = %s", filenamebaseIn.c_str());
			FileNames2(vndefocus, filenamebaseIn, vinfilenamesTot); // create "total 2D array" of input filenames
			index_t ndefcurrent = 0;
			for (index_t na = 0; na < nangles; na++)
			{
				vvinfilenames[na].resize(vndefocus[na]); // number of defocus planes at the current illumination angle
				for (index_t n = 0; n < vndefocus[na]; n++) vvinfilenames[na][n] = vinfilenamesTot[ndefcurrent++];
			}
		}

		string filenamebaseOutNew("BAD_STRING"), filenamebaseOutCAmp("BAD_STRING"), filenamebaseOutDefocCAmp("BAD_STRING"); // don't use it, unless it is redefined later
		vector<string> voutfilenamesTot; // output filenames for the reconstructed 3D potential
		if (imodeInvLaplace == 2 || imodeReproj == 2 || imodePeaks == 2) // output the the renormalized 3D potential
		{
			std::filesystem::path apath(filenamebaseOut);
			filenamebaseOutNew = apath.filename().string();
			filenamebaseOutNew.insert(filenamebaseOutNew.find_last_of("."), "R");
			filenamebaseOutNew = folderAux + filenamebaseOutNew;
			printf("\nOutput file name base for the renormalized 3D potential = %s", filenamebaseOutNew.c_str());
			FileNames(1, noutdefocus, filenamebaseOutNew, voutfilenamesTot); // create 1D array of output filenames to save 2D slices of the filtered and/or rescaled 3D object
		}
		else
		{
			printf("\nOutput file name base for 3D potential = %s", filenamebaseOut.c_str());
			FileNames(1, noutdefocus, filenamebaseOut, voutfilenamesTot); // create 1D array of output filenames to save 2D slices of the reconstructed 3D object
		}

		vector<string> voutfilenamesPeaksTot; // output filenames for the peak-localized reconstructed 3D potential
		if (imodePeaks == 1 || imodePeaks == 2) // output filenames for the peak-localized 3D potential
		{
			std::filesystem::path apath(filenamebaseOut);
			filenamebaseOutNew = apath.filename().string();
			filenamebaseOutNew.insert(filenamebaseOutNew.find_last_of("."), "A");
			filenamebaseOutNew = folderAux + filenamebaseOutNew;
			printf("\nOutput file name base for the peak-localized 3D potential = %s", filenamebaseOutNew.c_str());
			FileNames(1, noutdefocus, filenamebaseOutNew, voutfilenamesPeaksTot); // create 1D array of output filenames to save 2D slices of the peak-localized 3D object
		}

		vector<string> voutfilenamesTotCAmp; // output filenames the reprojected defocused complex amplitudes
		vector< vector<string> > vvoutfilenamesCAmp(nangles);  // same output filenames for the reprojected defocused complex amplitudes in the form of vector of vectors
		if (imodeReproj == 1 || imodeReproj == 2) // create filenames for saving the reprojected intensities, based on the names of input defocused intensities
		{
			std::filesystem::path apath(filenamebaseIn);
			filenamebaseOutCAmp = apath.filename().string();
			filenamebaseOutCAmp.replace(filenamebaseOutCAmp.find_last_of("."), filenamebaseOutCAmp.length() - filenamebaseOutCAmp.find_last_of("."), "R.grc");
			filenamebaseOutCAmp = folderAux + filenamebaseOutCAmp;
			printf("\nOutput file name base for reprojected complex amplitudes = %s", filenamebaseOutCAmp.c_str());
			FileNames2(vndefocus, filenamebaseOutCAmp, voutfilenamesTotCAmp); // create "2D array" of output filenames for reprojected complex amplitudes
			
			// Also copy the filenames of the input defocused intensities for the purpose of comparing with the reprojected intensities
			// note that in the case imodeInvLaplace == 2 or imodeReproj == 2 or imodePeaks == 2, these files are different from input images (which come from the potential)
			printf("\nDefocus series file name base for comparing with re-projected images = %s", filenamebaseIn.c_str());
			FileNames2(vndefocus, filenamebaseIn, vinfilenames1Tot); // create "total 2D array" of input filenames

			// Create the vector of vectors of output filenames for reprojected complex amplitudes at different defocus distances at each illumination angle
			// and the vector of vectors of filenames of the input defocused images for the purpose of comparing with the reprojected intensities
			index_t ndefcurrent = 0;
			for (index_t na = 0; na < nangles; na++)
			{
				vvoutfilenamesCAmp[na].resize(vndefocus[na]); // number of defocus planes at the current illumination angle
				vvinfilenames1[na].resize(vndefocus[na]); // number of defocus planes at the current illumination angle
				for (index_t n = 0; n < vndefocus[na]; n++)
				{
					vvoutfilenamesCAmp[na][n] = voutfilenamesTotCAmp[ndefcurrent];
					vvinfilenames1[na][n] = vinfilenames1Tot[ndefcurrent];
					ndefcurrent++;
				}
			}
		}

		vector<string> voutfilenamesTotDefocCAmp; // output filenames the phase-retrieved defocused complex amplitudes
		if (nSaveDefocCAmps == 1) // create filenames for saving the phase-retrieved defocused complex amplitudes
		{
			std::filesystem::path apath(filenamebaseIn);
			filenamebaseOutDefocCAmp = apath.filename().string();
			filenamebaseOutDefocCAmp.replace(filenamebaseOutDefocCAmp.find_last_of("."), filenamebaseOutDefocCAmp.length() - filenamebaseOutDefocCAmp.find_last_of("."), "D.grc");
			filenamebaseOutDefocCAmp = folderAux + filenamebaseOutDefocCAmp;
			printf("\nOutput file name base for phase-retrieved defocused complex amplitudes = %s", filenamebaseOutDefocCAmp.c_str());
			FileNames(nangles, 1, filenamebaseOutDefocCAmp, voutfilenamesTotDefocCAmp);

		}

		//************************************ end creating vectors of input and output file names


		//*********************************** start main calculations
		bool bAbort(false);

		std::unique_ptr<IXAHead> pHead(nullptr);
		XArray3D<double> K3out; // big 3D reconstructed array (needs to fit into RAM alongside with with everything else)

		if (imodeInvLaplace != 2 && imodeReproj != 2 && imodePeaks != 2) // do phase retrieval and backpropagation prior to 3D filtering or reprojection and output
		{
			printf("\n\nPerforming phase retrieval and backpropagation into the 3D volume ...");

			if (bRAWinput) // in this case, the parameters in question cannot possibly change between different input files
			{
				pHead.reset(new Wavehead2D(wl, ylo, yhi, xlo, xhi));
				xst = (nx > 1) ? (xhi - xlo) / double(nx) : 1.0;
				yst = (ny > 1) ? (yhi - ylo) / double(ny) : 1.0;
			}

			// start of cycle over illumination angles
			//index_t ndefcurrent(0);
			for (index_t na = 0; na < nangles; na++)
			{
				double angleZ = v2angles[na].a * PI180;
				double angleY = v2angles[na].b * PI180;
				double cosangleY = cos(angleY);
				double sinangleY = sin(angleY);
				double cosangleZ = cos(angleZ);
				double sinangleZ = sin(angleZ);
				if (bVerboseOutput) printf("\n\n*** Illumination angle[%zd] = (%g, %g) (degrees)", na, angleZ / PI180, angleY / PI180);
				else printf("\n*** Illumination angle[%zd] = (%g, %g) (degrees)", na, angleZ / PI180, angleY / PI180);

				index_t ndefocus = vndefocus[na]; // number of defocus planes at the current illumination angle
				vector<Pair> vdefocus = vvdefocus[na]; // vector of input Z" rotation angles and defocus positions at the current illumination angle
				vector<string> vinfilenames = vvinfilenames[na]; // input filenames of defocused images at the current illumination angle
				vector<string> voutfilenames(noutdefocus);

				double zout(0.0); // z plane at which campOut will be defined
				XArray2D<dcomplex> campOut; // complex amplitude in the "middle" defocus plane

				// do phase retrieval (if required)
				if (GetFileExtension(vinfilenames[0]) == string(".GRC")) // phase retrieval is not required - work with the provided complex amplitude
				{
					if (bVerboseOutput) printf("\nz'' rotation angle = %g (degrees), defocus distance = %g (Angstroms)", vdefocus[0].a, vdefocus[0].b);
					zout = vdefocus[0].b;
					if (bVerboseOutput) printf("\nReading input file %s ...", vinfilenames[0].c_str());
					XArData::ReadFileGRC(campOut, vinfilenames[0].c_str()); //	read the only one (!) input GRC file
					if (dNormFactor1 != 1.0) campOut /= dcomplex(dNormFactor1, 0.0);
					if (dNormFactor2 != 0.0) campOut += dcomplex(dNormFactor2, 0.0);
					pHead.reset(campOut.GetHeadPtr()->Clone());
					IXAHWave2D* ph2 = GetIXAHWave2D(campOut);
					ny = campOut.GetDim1();
					nx = campOut.GetDim2();
					xlo = ph2->GetXlo();
					xhi = ph2->GetXhi();
					xst = ph2->GetXStep(nx);
					ylo = ph2->GetYlo();
					yhi = ph2->GetYhi();
					yst = ph2->GetYStep(ny);
					if (abs(xst - yst) > 0.1 * xst / nx || abs(xst - zst) > 0.1 * xst / nx)	throw std::exception("Error: the 3D reconstructed object is supposed to have cubic voxels");
					if (vdefocus[0].a != 0 || bRelion && (v2shifts[na].a != 0 || v2shifts[na].b != 0)) // rotate around Z'', then shift along XY, if needed
					{
						XArray2D<double> ampRe0, ampIm0, ampRe, ampIm;
						Re(campOut, ampRe0); Im(campOut, ampIm0);
						double averRe0 = ampRe0.NormAverEdge(5);
						double averIm0 = ampIm0.NormAverEdge(5);
						// shift along X and/or Y back to the unshifted position
						if (bRelion && (v2shifts[na].a != 0 || v2shifts[na].b != 0))
						{
							xar::XArray2DMove<double> xamoveRe(ampRe0);
							xamoveRe.Move((long)floor(-v2shifts[na].b / yst + 0.5), (long)floor(-v2shifts[na].a / xst + 0.5), averRe0);
							xar::XArray2DMove<double> xamoveIm(ampIm0);
							xamoveIm.Move((long)floor(-v2shifts[na].b / yst + 0.5), (long)floor(-v2shifts[na].a / xst + 0.5), averIm0);
						}
						// rotate input defocused complex amplitude around Z'' back to zero angle
						if (vdefocus[0].a != 0)
						{
							XArray2DSpln<double> xaSplnRe(ampRe0), xaSplnIm(ampIm0);
							xaSplnRe.Rotate(ampRe, -vdefocus[0].a, 0.5 * (yhi + ylo), 0.5 * (xhi + xlo), averRe0); // expecting uniform background
							xaSplnIm.Rotate(ampIm, -vdefocus[0].a, 0.5 * (yhi + ylo), 0.5 * (xhi + xlo), averIm0); // expecting uniform background
							MakeComplex(ampRe, ampIm, campOut, false);
						}
						else
							MakeComplex(ampRe0, ampIm0, campOut, false);

					}

					// optionally save rotated defocused and shifted complex amplitudes in GRC files
					if (nSaveDefocCAmps == 1)
						XArData::WriteFileGRC(campOut, voutfilenamesTotDefocCAmp[na].c_str(), eGRCBIN);
				}
				else // do phase retrieval and find the complex wave amplitude in the plane z = zout
				{
					if (ndefocus == 1)
					{
						if (nPhaseRetrieval != 3 && nPhaseRetrieval != 4)
							throw std::exception("Error: unsuitable phase retrieval method in the input parameter file for the case of single defocus distance.");
					}
					else
					{
						if (nPhaseRetrieval != 1 && nPhaseRetrieval != 2)
							throw std::exception("Error: unsuitable phase retrieval method in the input parameter file for the case of multiple defocus distances.");
					}

					// read defocused images from files
					vector<XArray2D<double>> vint0(ndefocus); // initial defocused intensities
					for (int n = 0; n < ndefocus; n++)
					{
						if (bVerboseOutput) printf("\nz'' rotation angle = %g (degrees), defocus distance = %g (Angstroms)", vdefocus[n].a, vdefocus[n].b);
						if (bVerboseOutput) printf("\nReading input file %s ...", vinfilenames[n].c_str());
						if (bRAWinput)
						{
							XArData::ReadFileRAW(vint0[n], vinfilenames[n].c_str(), ny, nx, nHeaderLength, nElementLength, bBigEndian);
							vint0[n].SetHeadPtr(pHead->Clone());
						}
						else if (bTIFFinput)
						{
							TIFFReadFile(vint0[n], vinfilenames[n].c_str()); // read input TIFF files
							if (n == 0)
							{
								pHead.reset(new Wavehead2D(wl, ylo, yhi, xlo, xhi));
								ny = vint0[0].GetDim1();
								nx = vint0[0].GetDim2();
								xst = (nx > 1) ? (xhi - xlo) / double(nx) : 1.0;
								yst = (ny > 1) ? (yhi - ylo) / double(ny) : 1.0;
							}
							else
							{
								ny = vint0[n].GetDim1();
								nx = vint0[n].GetDim2();
								if (ny != vint0[0].GetDim1() || nx != vint0[0].GetDim2())
									throw std::exception("Error: input TIFF files have different dimensions");
							}
							vint0[n].SetHeadPtr(pHead->Clone());
						}
						else
						{
							XArData::ReadFileGRD(vint0[n], vinfilenames[n].c_str(), wl); //	read input GRD files
							if (n == 0) pHead.reset(vint0[0].GetHeadPtr()->Clone());
							IXAHWave2D* ph2 = GetIXAHWave2D(vint0[n]);
							if (ph2 == 0) throw std::exception("Error: input GRD file does not have a Wave2D head.");
							ny = vint0[n].GetDim1();
							nx = vint0[n].GetDim2();
							xlo = ph2->GetXlo();
							xhi = ph2->GetXhi();
							xst = ph2->GetXStep(nx);
							ylo = ph2->GetYlo();
							yhi = ph2->GetYhi();
							yst = ph2->GetYStep(ny);
							if (abs(xst - yst) > 0.1 * xst / nx || abs(xst - zst) > 0.1 * xst / nx)	throw std::exception("Error: the 3D reconstructed object is supposed to have cubic voxels");
						}

						// renormalize input data
						if (dNormFactor1 != 1.0) vint0[n] /= dNormFactor1;
						if (dNormFactor2 != 0.0) vint0[n] += dNormFactor2;
						if (vint0[n].Norm(eNormMin) < 0)
						{
							printf("\nWARNING: negative values in the input intensity file.");
							vint0[n].ThresholdLow(0.0, 0.0);
						}

						// rotate input defocused image around Z'' back to zero angle
						if (vdefocus[n].a != 0 || bRelion && (v2shifts[na].a != 0 || v2shifts[na].b != 0))
						{
							double aver = vint0[n].NormAverEdge(5);
							// shift along X and/or Y back to the unshifted position
							if (bRelion && (v2shifts[na].a != 0 || v2shifts[na].b != 0))
							{
								XArray2DMove<double> xamove(vint0[n]);
								xamove.Move((long)floor(-v2shifts[na].b / yst + 0.5), (long)floor(-v2shifts[na].a / xst + 0.5), aver);
							}
							// rotate input defocused complex amplitude around Z'' back to zero angle
							if (vdefocus[n].a != 0)
							{
								XArray2D<double> vintTemp(vint0[n]);
								XArray2DSpln<double> xaSpln(vintTemp);
								xaSpln.Rotate(vint0[n], -vdefocus[n].a, 0.5 * (yhi + ylo), 0.5 * (xhi + xlo), aver); // expecting uniform background
							}
						}
					}

					// do phase retrieval
					switch (nPhaseRetrieval)
					{
					case 1:
					{
						vector<double> vdefocusdist(ndefocus);
						for (int n = 0; n < ndefocus; n++) vdefocusdist[n] = vdefocus[n].b;
						if (bVerboseOutput) printf("\nPerforming IWFR phase retrieval ...");
						zout = 0.0;
						for (size_t n = 0; n < ndefocus; n++) zout += vdefocus[n].b;
						zout /= double(ndefocus);
						XA_IWFR<double> xa_iwfr;
						xa_iwfr.Iwfr(vint0, campOut, vdefocusdist, zout, k2maxo, Cs3, Cs5, itermax, epsilon, bVerboseOutput);
					}
					break;
					case 2:
					{
						vector<double> vdefocusdist(ndefocus);
						for (int n = 0; n < ndefocus; n++) vdefocusdist[n] = vdefocus[n].b;
						if (bVerboseOutput) printf("\nPerforming CTFL2 phase retrieval ...");
						zout = vdefocus[0].b; //CTFL2 always retrieves the phase in the first defocused plane
						XA_IWFR<double> xa_iwfr;
						XArray2D<double> fiOut, ampOut;
						ampOut = vint0[0]; // use this before the call to xa_iwfr.CTFL2 spoils all vint0 images
						ampOut ^= 0.5; // preserve the amplitude at the zeros distance for reuse after phase retrieval
						double alphaCTF2 = epsilon; // alphaCTF2 will be changed by the call to CTFL2()
						xa_iwfr.CTFL2(vint0, fiOut, vdefocusdist, k2maxo, Cs3, Cs5, alphaCTF2);
						if (bVerboseOutput) printf("\nMinimal [sum(CTF_n^2)]^2 in denominator = %g", alphaCTF2);
						MakeComplex(ampOut, fiOut, campOut, true);
					}
					break;
					case 3:
					{
						zout = vdefocus[0].b;
						XArray2D<double> fiOut, ampOut;
						ampOut = vint0[0]; ampOut ^= 0.5;
						XA_IWFR<double> xa_iwfr;
						if (bVerboseOutput) printf("\nPerforming -0.5*Log(I) phase retrieval ...");
						xa_iwfr.MinLogAmp(vint0[0], fiOut);
						MakeComplex(ampOut, fiOut, campOut, true);
					}
					break;
					case 4:
					{
						zout = vdefocus[0].b;
						XArray2D<double> fiOut, ampOut;
						ampOut = vint0[0]; ampOut ^= 0.5;
						XA_IWFR<double> xa_iwfr;
						if (bVerboseOutput) printf("\nPerforming -0.5*sqrt(Kmax^2-K^2) phase retrieval ...");
						xa_iwfr.PhaseB7(vint0[0], fiOut);
						MakeComplex(ampOut, fiOut, campOut, true);
					}
					break;
					default:
						throw std::exception("Error: unknown phase retrieval method in the input parameter file.");
					}

					// optionally save phase retrieved defocused complex amplitudes in GRC files
					if (nSaveDefocCAmps == 1)
						XArData::WriteFileGRC(campOut, voutfilenamesTotDefocCAmp[na].c_str(), eGRCBIN);
				}

				// now start calculations of the output defocused images
				if (bVerboseOutput) printf("\nPropagating to output defocus planes ...");
				vector<XArray2D<dcomplex>> vcamp(noutdefocus); // defocused complex amplitudes
				#pragma omp parallel for
				for (int n = 0; n < noutdefocus; n++)
				{
					int iSign = voutdefocus[n] < zout ? iSign = -1 : iSign = 1;
					try
					{
						// propagate to the output defocused plane
						vcamp[n] = campOut;
						xar::XArray2DFFT<double> xafft(vcamp[n]);

						if (bRelion)
							xafft.FresnelA(voutdefocus[n] - vastigm[na] - zout, voutdefocus[n] + vastigm[na] - zout, false, k2maxo, iSign * Cs3, iSign * Cs5); // propagate to z_out[n]
						else
							xafft.Fresnel(voutdefocus[n] - zout, false, k2maxo, iSign * Cs3, iSign * Cs5); // propagate to z_out[n]
					}
					catch (std::exception& E)
					{
						printf("\n\n!!!Exception: %s\n", E.what());
						bAbort = true;
					}
				}
				if (bAbort) throw std::runtime_error("at least one thread has thrown an exception.");

				if (bVerboseOutput) printf("\nUpdating 3D reconstructed object ...");
				double xc = (xhi + xlo) / 2.0; // x-coordinate of the centre of rotation
				double yc = (yhi + ylo) / 2.0; // y-coordinate of the centre of rotation
				double zc = (zhidz + zlodz) / 2.0; // z-coordinate of the centre of rotation

				// calculate the coordinate illumination angle parameters
				double xxx;
				vector<double> x_sinangleY(nx), x_cosangleY(nx);
				for (index_t i = 0; i < nx; i++)
				{
					xxx = xlo + xst * i - xc;
					x_sinangleY[i] = xxx * sinangleY;
					x_cosangleY[i] = xxx * cosangleY;
				}

				double yyy;
				vector<double> y_sinangleZ(ny), y_cosangleZ(ny);
				for (index_t j = 0; j < ny; j++)
				{
					yyy = ylo + yst * j - yc;
					y_sinangleZ[j] = yyy * sinangleZ;
					y_cosangleZ[j] = yyy * cosangleZ;
				}
				double zzz;
				vector<double> z_sinangleY(noutdefocus), z_cosangleY(noutdefocus);
				for (index_t n = 0; n < noutdefocus; n++)
				{
					zzz = voutdefocus[n] - zc;
					z_sinangleY[n] = zzz * sinangleY;
					z_cosangleY[n] = zzz * cosangleY;
				}

				// allocate the large 3D output array
				if (na == 0)
				{
					K3out.Resize(noutdefocus, ny, nx, 0.0);
					K3out.SetHeadPtr(new xar::Wavehead3D(wl, zlo, zhi, ylo, yhi, xlo, xhi));
				}

				int ii, jj, nn;
				double xx;
				int nx2 = int(nx) - 2, ny2 = int(ny) - 2, noutdefocus2 = noutdefocus - 2;
				double xx_sinangleZ, xx_cosangleZ, dK, dx0, dx1, dy0, dy1, dz0, dz1, dz0K, dz1K;

				#pragma omp parallel for shared (K3out, vcamp) private(xx, xxx, yyy, zzz, xx_sinangleZ, xx_cosangleZ, dx1, dx0, dy1, dy0, dz1, dz0, ii, jj, nn, dK, dz0K, dz1K)
				for (int n = 0; n < noutdefocus; n++)
				{
					try
					{
						for (int i = 0; i < nx; i++)
						{
							// inverse rotation around Y' axis
							zzz = zc + z_cosangleY[n] + x_sinangleY[i]; // z coordinate after the inverse rotation around Y'
							dz1 = abs(zzz - zlodz) / zst;
							nn = (int)dz1; if (nn < 0 || nn > noutdefocus2) continue;
							dz1 -= nn; dz0 = 1.0 - dz1;

							xx = -z_sinangleY[n] + x_cosangleY[i]; // x - xc coordinate after the inverse rotation around Y'
							xx_sinangleZ = xx * sinangleZ;
							xx_cosangleZ = xx * cosangleZ;

							for (int j = 0; j < ny; j++)
							{

								// inverse rotation around Z axis
								xxx = xc + y_sinangleZ[j] + xx_cosangleZ; // x coordinate after the inverse rotation around Z
								dx1 = abs(xxx - xlo) / xst;
								ii = (int)dx1; if (ii < 0 || ii > nx2) continue;
								dx1 -= ii; dx0 = 1.0 - dx1;

								yyy = yc + y_cosangleZ[j] - xx_sinangleZ; // y coordinate after the inverse rotation around Z
								dy1 = abs(yyy - ylo) / yst;
								jj = (int)dy1; if (jj < 0 || jj > ny2) continue;
								dy1 -= jj; dy0 = 1.0 - dy1;

								dK = 1.0 - std::norm(vcamp[n][j][i]);
								dz0K = dz0 * dK;
								dz1K = dz1 * dK;
								K3out[nn][jj][ii] += dx0 * dy0 * dz0K;
								K3out[nn][jj][ii + 1] += dx1 * dy0 * dz0K;
								K3out[nn + 1][jj][ii] += dx0 * dy0 * dz1K;
								K3out[nn + 1][jj][ii + 1] += dx1 * dy0 * dz1K;
								K3out[nn][jj + 1][ii] += dx0 * dy1 * dz0K;
								K3out[nn][jj + 1][ii + 1] += dx1 * dy1 * dz0K;
								K3out[nn + 1][jj + 1][ii] += dx0 * dy1 * dz1K;
								K3out[nn + 1][jj + 1][ii + 1] += dx1 * dy1 * dz1K;
							}
						}
					}
					catch (std::exception& E)
					{
						printf("\n\n!!!Exception: %s\n", E.what());
						bAbort = true;
					}
				}
				if (bAbort) throw std::runtime_error("at least one thread has thrown an exception.");
			} // end of cycle over illumination angles

			// normalization of the reconstructed 3D distribution
			double dnorm = -2.0 * EE / (dzextra * dAtomWidth * nangles); // follows from eq.(11!) in our second Ultramicroscopy paper (except for the minus sign)
			if (imodeInvLaplace == 0) dnorm *= dAtomWidth * dAtomWidth / (4.0 * PI * PI); // we select the renormalization effect that the inverse Laplacian would have had on the frequency 1 / dAtomWidth.
			K3out *= dnorm;

		} // end of case if imodeInvLaplace != 2 && imodePeaks != 2 && imodeReproj != 2
		else // read in pre-existing 3D distribution of the electrostatic potential
		{
			printf("\n\nReading pre-existing 3D distribution of the electrostatic potential from files ...");

			// read input 2D slices
			XArray2D<double> ipIn;
			for (int n = 0; n < noutdefocus; n++)
			{
				if (bRAWinput)
				{
					XArData::ReadFileRAW(ipIn, vinfilenamesTot[n].c_str(), ny, nx, nHeaderLength, nElementLength, bBigEndian);
					if (dNormFactor1 != 1.0) ipIn /= dNormFactor1;
					if (dNormFactor2 != 0.0) ipIn += dNormFactor2;
					if (n == 0)
					{
						pHead.reset(new Wavehead2D(wl, ylo, yhi, xlo, xhi));
						xst = (nx > 1) ? (xhi - xlo) / double(nx) : 1.0;
						yst = (ny > 1) ? (yhi - ylo) / double(ny) : 1.0;
						if (abs(xst - yst) > 0.1 * xst / nx || abs(xst - zst) > 0.1 * xst / nx)	throw std::exception("Error: the input 3D object is supposed to have cubic voxels");
						K3out.Resize(noutdefocus, ny, nx, 0.0);
						K3out.SetHeadPtr(new xar::Wavehead3D(wl, zlo, zhi, ylo, yhi, xlo, xhi));
					}
				}
				else if (bTIFFinput)
				{
					TIFFReadFile(ipIn, vinfilenamesTot[n].c_str()); // read input TIFF files
					if (n == 0)
					{
						pHead.reset(new Wavehead2D(wl, ylo, yhi, xlo, xhi));
						ny = ipIn.GetDim1();
						nx = ipIn.GetDim2();
						xst = (nx > 1) ? (xhi - xlo) / double(nx) : 1.0;
						yst = (ny > 1) ? (yhi - ylo) / double(ny) : 1.0;
						if (abs(xst - yst) > 0.1 * xst / nx || abs(xst - zst) > 0.1 * xst / nx)	throw std::exception("Error: the input 3D object is supposed to have cubic voxels");
						K3out.Resize(noutdefocus, ny, nx, 0.0);
						K3out.SetHeadPtr(new xar::Wavehead3D(wl, zlo, zhi, ylo, yhi, xlo, xhi));
					}
					else
					{
						if (ny != ipIn.GetDim1() || nx != ipIn.GetDim2())
							throw std::exception("Error: input TIFF files have different dimensions");
					}
				}
				else
				{

					XArData::ReadFileGRD(ipIn, vinfilenamesTot[n].c_str(), wl); //	read input GRD files
					if (n == 0)
					{
						IXAHWave2D* ph2 = GetIXAHWave2D(ipIn);
						ny = ipIn.GetDim1();
						nx = ipIn.GetDim2();
						xlo = ph2->GetXlo();
						xhi = ph2->GetXhi();
						xst = ph2->GetXStep(nx);
						ylo = ph2->GetYlo();
						yhi = ph2->GetYhi();
						yst = ph2->GetYStep(ny);
						if (abs(xst - yst) > 0.1 * xst / nx || abs(xst - zst) > 0.1 * xst / nx) throw std::exception("Error: the input 3D object is supposed to have cubic voxels");
						pHead.reset(ipIn.GetHeadPtr()->Clone());
						K3out.Resize(noutdefocus, ny, nx, 0.0);
						K3out.SetHeadPtr(new xar::Wavehead3D(wl, zlo, zhi, ylo, yhi, xlo, xhi));
					}
				}
				for (index_t j = 0; j < ny; j++)
					for (index_t i = 0; i < nx; i++)
						K3out[n][j][i] = ipIn[j][i];
			}
		}

		// apply the regularized inverse 3D Laplacian
		if (imodeInvLaplace == 1 || imodeInvLaplace == 2)
		{
			printf("\nInverse Laplace filtering the reconstructed 3D object ...");
			if (imodeInvLaplace == 2) K3out *= (4.0 * PI * PI) / (dAtomWidth * dAtomWidth); // compensate for the renormalization that was applied when it was imodeInvLaplace == 0
			Fftwd3drc fft3(noutdefocus, int(ny), int(nx));
			fft3.InverseLaplacian(K3out, alpha);
		}

		// high-pass filtering, background subtraction and thresholding of the reconstructed 3D distribution
		printf("\nHigh-pass filtering and thresholding the reconstructed 3D object ...");
		if (dlpfiltersize > 0 && dlpfiltersize < (xhi - xlo))
		{
			XArray3D<double> K3outTemp = K3out;
			Fftwd3drc fft3(noutdefocus, int(ny), int(nx));
			fft3.GaussFilter(K3outTemp, dlpfiltersize / 2.355);
			K3out -= K3outTemp;
		}
		K3out.ThresholdLow(dBackground, dThreshold);
		 
		// output the 3D array
		printf("\n\nSaving the reconstructed 3D object into output files ...");
		if (bTIFFoutput) TIFFWriteFileStack(K3out, voutfilenamesTot, eTIFF32);
		else XArData::WriteFileStackGRD(K3out, voutfilenamesTot, eGRDBIN);

		// multislice reprojection
		if (imodeReproj == 1 || imodeReproj == 2)
		{
			printf("\nPerforming forward re-propagation through the 3D distribution of the electrostatic potential ...");
			double yc = (ylo + yhi) / 2.0, xc = (xlo + xhi) / 2.0;
			double dnorm3 = PI / (wl * EE) * zst;
			K3out *= dnorm3; // this converts the potential V into -2 * PI / lambda  * delta * zst, with zst required in the subsequent integration inside Multislice_eV
			xar::XArray3DSpln<double> spln3d(K3out);

			vector<vector<double> > vverr(nangles); // reconstruction errors in individual defocus planes at different illumination angles
			vector<double> verraver(nangles, 0.0); // reconstruction errors at different illumination angles averaged over the corresponding defocus distances
			vector<vector<double> > vvnorm0(nangles); // C0 norms of reprojected intensities-1 in individual defocus planes at different illumination angles
			vector<double> vnormaver0(nangles, 0.0); // C0 norms of reprojected intensities-1 different illumination angles averaged over the corresponding defocus distances
			vector<vector<double> > vvnorm1(nangles); // C0 norms of original intensities-1 in individual defocus planes at different illumination angles
			vector<double> vnormaver1(nangles, 0.0); // C0 norms of original intensities-1 different illumination angles averaged over the corresponding defocus distances

			#pragma omp parallel for shared (K3out, spln3d, vndefocus, vvdefocus, v2angles, vvoutfilenamesCAmp)
			for (int na = 0; na < nangles; na++)
			{
				try
				{
					vector<Pair> vdefocus;
					vector<string> voutfilenamesCAmp, vinfilenames1;
					double angleY(0);
					double angleZ(0);
					#pragma omp critical
					{
						vdefocus = vvdefocus[na];
						vverr[na].resize(vndefocus[na]);
						vvnorm0[na].resize(vndefocus[na]);
						vvnorm1[na].resize(vndefocus[na]);
						voutfilenamesCAmp = vvoutfilenamesCAmp[na];
						vinfilenames1 = vvinfilenames1[na];
						angleZ = v2angles[na].a * PI180;
						angleY = v2angles[na].b * PI180;
						if (!bVerboseOutput) printf("\n*** Illumination angle[%d] = (%g, %g) (degrees)", na, angleZ / PI180, angleY / PI180);
					}

					// multislice propagate through the reconstructed 3D distribution of the electrostatic potential
					XArray2D<dcomplex> campOut(ny, nx, dcomplex(1.0, 0.0)); // "incident" complex amplitude
					campOut.SetHeadPtr(pHead->Clone());
					spln3d.Multislice_eV(campOut, angleZ, angleY, sliceTh, k2maxo);

					// propagate to defocuse planes, rotate around the illumination axis and save the defocused intensities
					XArray2D<double > vint0n, K0n, vint1n, K1n, ampRe0, ampIm0, ampRe, ampIm;
					for (int n = 0; n < vndefocus[na]; n++)
					{
						XArray2D<dcomplex> campOut1(campOut);
						xar::XArray2DFFT<double> xafft(campOut1);
						int iSign = vdefocus[n].b < 0 ? iSign = -1 : iSign = 1;
						if (bRelion)
							xafft.FresnelA(vdefocus[n].b + vastigm[na], vdefocus[n].b - vastigm[na], false, k2maxo, iSign * Cs3, iSign * Cs5); // propagate to the current defocus distance
						else
							xafft.Fresnel(vdefocus[n].b, false, k2maxo, iSign * Cs3, iSign * Cs5); // propagate to the current defocus distance
						if (vdefocus[n].a != 0 || bRelion && (v2shifts[na].a != 0 || v2shifts[na].b != 0)) // rotate around Z'', then shift along XY, if needed
						{
							Re(campOut1, ampRe0); Im(campOut1, ampIm0);
							double averRe0 = ampRe0.NormAverEdge(5);
							double averIm0 = ampIm0.NormAverEdge(5);
							// rotate input defocused complex amplitude around Z''
							if (vdefocus[0].a != 0)
							{
								XArray2DSpln<double> xaSplnRe(ampRe0), xaSplnIm(ampIm0);
								xaSplnRe.Rotate(ampRe, vdefocus[n].a, yc, xc, averRe0); // expecting uniform background
								xaSplnIm.Rotate(ampIm, vdefocus[n].a, yc, xc, averIm0); // expecting uniform background
							}
							else
							{
								ampRe = ampRe0;
								ampIm = ampIm0;
							}
							// shift along X and/or Y
							if (bRelion && (v2shifts[na].a != 0 || v2shifts[na].b != 0))
							{
								xar::XArray2DMove<double> xamoveRe(ampRe);
								xamoveRe.Move((long)floor(v2shifts[na].b / yst + 0.5), (long)floor(v2shifts[na].a / xst + 0.5), averRe0);
								xar::XArray2DMove<double> xamoveIm(ampIm);
								xamoveIm.Move((long)floor(v2shifts[na].b / yst + 0.5), (long)floor(v2shifts[na].a / xst + 0.5), averIm0);
							}
							MakeComplex(ampRe, ampIm, campOut1, false);
						}
						Abs(campOut1, vint1n);
						Abs2(campOut1, K1n); K1n -= 1.0; 
						// read in the original intensity
						if (GetFileExtension(vinfilenames1[n]) == string(".GRC"))
						{
							XArray2D<dcomplex> campOut2;
							XArData::ReadFileGRC(campOut2, vinfilenames1[n].c_str());
							Abs(campOut2, vint0n);
							Abs2(campOut2, K0n); K0n -= 1.0;
						}
						else
						{
							if (GetFileExtension(vinfilenames1[n]) == string(".TIF") || GetFileExtension(vinfilenames1[n]) == string(".TIFF"))
							{
								TIFFReadFile(vint0n, vinfilenames1[n].c_str());
								vint0n.SetHeadPtr(pHead->Clone());
							}
							else
								XArData::ReadFileGRD(vint0n, vinfilenames1[n].c_str(), wl);
							K0n = vint0n; K0n -= 1.0;
							vint0n ^= 0.5;
						}
						// compare the reprojected and the original defocused intensities
						vvnorm0[na][n] = K0n.Norm(eNormC0);
						vnormaver0[na] += vvnorm0[na][n];
						vvnorm1[na][n] = K1n.Norm(eNormC0);
						vnormaver1[na] += vvnorm1[na][n];
						vint1n -= vint0n;
						vverr[na][n] = pow(vint1n.Norm(eNormL2), 2.0) / pow(vint0n.Norm(eNormL2), 2.0); // this error norm is modelled after the one in IWFR
						verraver[na] += vverr[na][n];

						//ReplaceModulus(campOut1, vint0n);
						//Conjg(campOut1); // @@@@@@@@@@@@@@@@ ?????

						#pragma omp critical
						{
							if (bVerboseOutput)
							{
								printf("\n\n*** Illumination angle[%d] = (%g, %g) (degrees)", na, angleZ / PI180, angleY / PI180);
								printf("\nz'' rotation angle = %g (degrees), defocus distance = %g (Angstroms)", vdefocus[n].a, vdefocus[n].b);
								printf("\nWriting output file %s ...", voutfilenamesCAmp[n].c_str());
							}
						}
						XArData::WriteFileGRC(campOut1, voutfilenamesCAmp[n].c_str(), eGRCBIN); //	write output GRC files
					}
					verraver[na] /= double(vndefocus[na]);
					vnormaver0[na] /= double(vndefocus[na]);
					vnormaver1[na] /= double(vndefocus[na]);
				}
				catch (std::exception& E)
				{
					printf("\n\n!!!Exception: %s\n", E.what());
					bAbort = true;
				}
			}
			if (bAbort) throw std::runtime_error("at least one thread has thrown an exception.");

			// calculate the current reconstruction error and
			// if the difference with the previous error is smaller than the defined minimum
			// or, if the error started to increase, interrupt the iterations
			double ssej = 0.0;
			for (index_t na = 0; na < nangles; na++) ssej += verraver[na];
			ssej /= double(nangles);

			printf("\n\nAverage relative error in reprojected defocused intensities = %g", ssej);
			for (index_t na = 0; na < nangles; na++) printf("\nIllum_angle_no. = %zd, C0(orig_contrast) = %g, C0(reproj_contrast) = %g, Rel.error = %g ", na, vnormaver0[na], vnormaver1[na], verraver[na]);

			// restore the 3D potential V in case we need to work with it further
			dnorm3 = 1.0 / dnorm3;
			K3out *= dnorm3;
		}

		// finding atomic positions
		if (imodePeaks == 1 || imodePeaks == 2)
		{
			int katom = int(datomsize / zst + 0.5), jatom = int(datomsize / yst + 0.5), iatom = int(datomsize / xst + 0.5);
			double datomsize2 = datomsize * datomsize;
			int natom(0);
			vector<int> vimax, vjmax, vkmax, vimax1, vjmax1, vkmax1, Znum1;
			vector<float> xa, ya, za, xa1, ya1, za1, occ1, wobble1;

			printf("\nSearching for peak positions in the 3D distribution of the electrostatic potential ...");
			#pragma omp parallel for shared(K3out, natom, vimax, vjmax, vkmax)
			for (int k = 0; k < noutdefocus - katom; k += katom)
			{
				for (int j = 0; j < ny - jatom; j += jatom)
					for (int i = 0; i < nx - iatom; i += iatom)
					{
						double K3max(0.0);
						int kmax(0), jmax(0), imax(0);
						for (int kk = k; kk < k + katom; kk++)
							for (int jj = j; jj < j + jatom; jj++)
								for (int ii = i; ii < i + iatom; ii++)
									if (K3out[kk][jj][ii] > K3max)
									{
										K3max = K3out[kk][jj][ii];
										K3out[kmax][jmax][imax] = 0.0;
										kmax = kk; jmax = jj; imax = ii;
									}
									else K3out[kk][jj][ii] = 0.0;
						#pragma omp critical
						if (K3max > 0)
						{
							natom++;
							vimax.push_back(imax);
							vjmax.push_back(jmax);
							vkmax.push_back(kmax);
						}
					}
			}

			if (natom == 0) printf("WARNING: no peaks have been found!");
			else
			{

				// create vectors of peak location coordinates
				xa.resize(natom); ya.resize(natom); za.resize(natom);
				vector<Pair2> K3maxPair(natom);
				for (int n = 0; n < natom; n++)
				{
					xa[n] = float(xlo + xst * vimax[n]);
					ya[n] = float(ylo + yst * vjmax[n]);
					za[n] = float(zst * vkmax[n]);
					K3maxPair[n].v = K3out[vkmax[n]][vjmax[n]][vimax[n]];
					K3maxPair[n].n = n;
				}

				// sort the located peaks according to their height, in order to improve the subsequent procedure of elimination of closely located peaks
				printf("\nSorting the located peaks in descending order according to their heights ...");
				std::qsort(&K3maxPair[0], (size_t)natom, sizeof(Pair2), Pair2comp);

				// exclude the smaller one from each pair of peak positions that are located closer than datomsize to each other (e.g. in adjacent corners of neigbouring cubes)
				printf("\nEliminating adjacent peaks in the 3D distribution of the electrostatic potential ...");
				int n, m;
				for (int nn = natom - 1; nn >= 0; nn--)
				{
					if (K3maxPair[nn].v != 0.0)
					{
						// we can already count this maximum in, as it is guaranteed to be larger than all subsequent ones
						n = K3maxPair[nn].n;
						Znum1.push_back(6);
						vimax1.push_back(vimax[n]);
						vjmax1.push_back(vjmax[n]);
						vkmax1.push_back(vkmax[n]);
						xa1.push_back(xa[n]);
						ya1.push_back(ya[n]);
						za1.push_back(za[n]);
						occ1.push_back((float)K3maxPair[nn].v);
						wobble1.push_back(0);
						#pragma omp parallel for private(m)
						for (int mm = nn - 1; mm >= 0; mm--)
						{
							if (K3maxPair[mm].v != 0.0)
							{
								m = K3maxPair[mm].n;
								if ((xa[n] - xa[m]) * (xa[n] - xa[m]) + (ya[n] - ya[m]) * (ya[n] - ya[m]) + (za[n] - za[m]) * (za[n] - za[m]) < datomsize2)
									K3maxPair[mm].v = 0.0;
							}
						}
					}
				}
				natom = (int)Znum1.size();
				printf("\n%d peak positions have been found in total.", natom);

				// copy the located isolated peaks back into the 3D potential distribution
				K3out.Fill(0.0);
				for (int n = 0; n < natom; n++)
					K3out[vkmax1[n]][vjmax1[n]][vimax1[n]] = occ1[n];

				// exclude the peaks located outside the minimal ball
				//double xac(50.0), yac(50.0), zac(50.0), xarad(41.0);
				//double xarad2 = xarad * xarad;
				//for (int n = 0; n < natom; n++)
				//		if ((xa[n] - xac) * (xa[n] - xac) + (ya[n] - yac) * (ya[n] - yac) + (za[n] - zac) * (za[n] - zac) > xarad2) 
				//			K3out[vkmax[n]][vjmax[n]][vimax[n]] = 0.0f;

				printf("\nSaving the reconstructed peak-localized 3D object into output files ...");
				if (bTIFFoutput) TIFFWriteFileStack(K3out, voutfilenamesPeaksTot, eTIFF32);
				else XArData::WriteFileStackGRD(K3out, voutfilenamesPeaksTot, eGRDBIN);

				if (natom > NATOMMAX)
					printf("\nNot saving the reconstructed peak-localized 3D object into an XYZ file, since too many peaks (%d) have been found.", natom);
				else
				{
					printf("\nSaving the reconstructed peak-localized 3D object into an XYZ file ...");
					size_t dotpos = filenamebaseOutNew.find_last_of(".");
					string fileoutXYZ = filenamebaseOutNew.replace(dotpos + 1, filenamebaseOutNew.size() - dotpos - 1, "xyz");
					SaveXYZfile(fileoutXYZ, float(xhi - xlo), float(zhi - zlo), (int)natom, &Znum1[0], &xa1[0], &ya1[0], &za1[0], &occ1[0], &wobble1[0]);
				}
			}
		} // end of peak localization module

	}
	catch (std::exception& E)
	{
		printf("\n\n!!!Exception: %s\n", E.what());
	}

    std::chrono::system_clock::time_point end_time = std::chrono::system_clock::now();
	printf("\n\nMain program finished. Execution time = %I64d s.", std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count());

	printf("\nPress any key to exit..."); getchar();
	return 0;

}