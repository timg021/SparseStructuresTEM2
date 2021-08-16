// Implementation of the fftwd3drc class

#include "fftwd3drc.h"
#include "XA_head3.h"

void Fftwd3drc::Cleanup()
{
	if (aplan != 0) { fftw_destroy_plan(aplan); aplan = 0; }
	if (bplan != 0) { fftw_destroy_plan(bplan); bplan = 0; }
	if (pin != 0) { fftw_free(pin); pin = 0; }
	if (pout != 0) { fftw_free(pout); pout = 0; }
	nz = ny = nx = 0;
}


void Fftwd3drc::GetRealXArray3D(xar::XArray3D<double>& aaa)
{
	if ((int)aaa.GetDim1() != nz || (int)aaa.GetDim2() != ny || (int)aaa.GetDim3() != nx)
		throw std::runtime_error("dimensions of the output real array are different from the internal one in Fftwd3drc");
	int m = 0;
	double fnorm = 1.0 / double(GetNr());
	for (index_t k = 0; k < nz; k++)
		for (index_t j = 0; j < ny; j++)
			for (index_t i = 0; i < nx; i++)
				aaa[k][j][i] = pin[m++] * fnorm;
}


void Fftwd3drc::SetRealXArray3D(xar::XArray3D<double> aaa)
{
	if ((int)aaa.GetDim1() != nz || (int)aaa.GetDim2() != ny || (int)aaa.GetDim3() != nx)
		throw std::runtime_error("dimensions of the input real array are different from the internal one in Fftwd3drc");
	int m = 0;
	for (index_t k = 0; k < nz; k++)
		for (index_t j = 0; j < ny; j++)
			for (index_t i = 0; i < nx; i++)
				pin[m++] = aaa[k][j][i];
}


void Fftwd3drc::GetComplexXArray3D(xar::XArray3D<xar::dcomplex>& aaa)
{
	int nx2 = GetNx2();
	if ((int)aaa.GetDim1() != nz || (int)aaa.GetDim2() != ny || (int)aaa.GetDim3() != nx2)
		throw std::runtime_error("dimensions of the output complex array are different from the internal one in Fftwd3drc");
	int m = 0;
	for (index_t k = 0; k < nz; k++)
		for (index_t j = 0; j < ny; j++)
			for (index_t i = 0; i < nx2; i++)
			{ 
				aaa[k][j][i] = xar::dcomplex(pout[m][0], pout[m][1]); 
				m++; 
			}
}


void Fftwd3drc::SetComplexXArray3D(xar::XArray3D<xar::dcomplex> aaa)
{
	int nx2 = GetNx2();
	if ((int)aaa.GetDim1() != nz || (int)aaa.GetDim2() != ny || (int)aaa.GetDim3() != nx2)
		throw std::runtime_error("dimensions of the input complex array are different from the internal one in Fftwd3drc");
	int m = 0;
	for (index_t k = 0; k < nz; k++)
		for (index_t j = 0; j < ny; j++)
			for (index_t i = 0; i < nx2; i++)
			{
				pout[m][0] = aaa[k][j][i].real();
				pout[m][1] = aaa[k][j][i].imag();
				m++;
			}
}


void Fftwd3drc::PrintRealArray(const char* message)
{
	double* pin = GetReal();
	printf(message);
	int m = 0;
	for (index_t k = 0; k < nz; k++)
		for (index_t j = 0; j < ny; j++)
			for (index_t i = 0; i < nx; i++)
				printf("\npin[%zd,%zd,%zd] = %g", k, j, i, pin[m++]);
}


void Fftwd3drc::PrintComplexArray(const char* message)
{
	printf(message);
	int m = 0;
	int nx2 = GetNx2();
	for (index_t k = 0; k < nz; k++)
		for (index_t j = 0; j < ny; j++)
			for (index_t i = 0; i < nx2; i++)
			{
				printf("\npout[%zd,%zd,%zd] = (%g, %g)", k, j, i, pout[m][0], pout[m][1]); 
				m++; 
			}
}


void Fftwd3drc::InverseMLaplacian(xar::XArray3D<double>& xa3, double alpha)
// inverse 3D (-Laplacian) via multiplication of the FFT by 4 * PI^2 * (ksi^2 + eta^2 + dzeta^2)
// alpha is the usual Tikhonov regularization parameter
{
	if (nz != xa3.GetDim1() || ny != xa3.GetDim2() || nx != xa3.GetDim3())
		throw std::exception("Error: input XArray3D in InverseLaplacian() has wrong dimensions.");

	IXAHWave3D* ph3 = GetIXAHWave3D(xa3);
	if (ph3 == 0) throw std::exception("Error: input XArray3D in InverseLaplacian() does not have a Wave3D head.");
	double xlo = ph3->GetXlo();
	double xhi = ph3->GetXhi();
	double ylo = ph3->GetYlo();
	double yhi = ph3->GetYhi();
	double zlo = ph3->GetZlo();
	double zhi = ph3->GetZhi();

	double dnorm1 = xa3.Norm(xar::eNormL1);
	if (dnorm1 == 0.0) return;

	double alpha1 = alpha / (1.0 + alpha); // this is the weitht applied to the identity operator for regularization
	double alpha2 = 1.0 / (1.0 + alpha); // this is the weight applied to the inverse Laplacian for regulalization

	// FFT of input array
	SetRealXArray3D(xa3);
	ForwardFFT();

	/// multiply FFT of input array by the FFT version of regularized inverse 3D Laplacian
	double fact = alpha2 * 4.0 * xar::PI * xar::PI;
	double dksi2 = fact / ((xhi - xlo) * (xhi - xlo));
	double deta2 = fact / ((yhi - ylo) * (yhi - ylo));
	double dzeta2 = fact / ((zhi - zlo) * (zhi - zlo));
	double dtemp, dtemp1, dk2, djk2;

	fftw_complex* pout = GetComplex();
	int m = 0, nc2 = GetNx2();
	index_t k1, j1, nyd2 = ny / 2, nzd2 = nz / 2;
	for (index_t k = 0; k < nz; k++)
	{
		k <= nzd2 ? k1 = k : k1 = nz - k;
		dk2 = k1 * k1 * dzeta2 + alpha1;
		for (index_t j = 0; j < ny; j++)
		{
			j <= nyd2 ? j1 = j : j1 = ny - j;
			djk2 = j1 * j1 * deta2 + dk2;
			for (index_t i = 0; i < nc2; i++)
			{
				dtemp = i * i * dksi2 + djk2;
				dtemp != 0 ? dtemp1 = 1.0 / dtemp : dtemp1 = 0.0; // protection against division by zero
				pout[m][0] *= dtemp1;
				pout[m][1] *= dtemp1;
				m++;
			}
		}
	}

	// inverse FFT of the product
	InverseFFT();
	GetRealXArray3D(xa3);

	// restore L1 norm
	//xa3 *= dnorm1 / xa3.Norm(xar::eNormL1);
}


void Fftwd3drc::GaussFilter(xar::XArray3D<double>& xa3, double sigma)
// Convolution with a 3D Gaussian with the standard deviation equal to sigma
{
	if (nz != xa3.GetDim1() || ny != xa3.GetDim2() || nx != xa3.GetDim3())
		throw std::exception("Error: input XArray3D in GaussFilter() has wrong dimensions.");

	IXAHWave3D* ph3 = GetIXAHWave3D(xa3);
	if (ph3 == 0) throw std::exception("Error: input XArray3D in GaussFilter() does not have a Wave3D head.");
	double xlo = ph3->GetXlo();
	double xhi = ph3->GetXhi();
	double ylo = ph3->GetYlo();
	double yhi = ph3->GetYhi();
	double zlo = ph3->GetZlo();
	double zhi = ph3->GetZhi();

	double dnorm1 = xa3.Norm(xar::eNormL1);
	if (dnorm1 == 0.0) return;

	// FFT of input array
	SetRealXArray3D(xa3);
	ForwardFFT();

	/// multiply FFT of input array by the FFT transform of the 3D Gaussian
	double fact = 2.0 * xar::PI * xar::PI * sigma * sigma;
	double dksi2 = fact / ((xhi - xlo) * (xhi - xlo));
	double deta2 = fact / ((yhi - ylo) * (yhi - ylo));
	double dzeta2 = fact / ((zhi - zlo) * (zhi - zlo));
	double dtemp, dk2, djk2;

	fftw_complex* pout = GetComplex();
	int m = 0, nc2 = GetNx2();
	index_t k1, j1, nyd2 = ny / 2, nzd2 = nz / 2;
	for (index_t k = 0; k < nz; k++)
	{
		k <= nzd2 ? k1 = k : k1 = nz - k;
		dk2 = k1 * k1 * dzeta2;
		for (index_t j = 0; j < ny; j++)
		{
			j <= nyd2 ? j1 = j : j1 = ny - j;
			djk2 = j1 * j1 * deta2 + dk2;
			for (index_t i = 0; i < nc2; i++)
			{
				dtemp = exp(-(i * i * dksi2 + djk2));
				pout[m][0] *= dtemp;
				pout[m][1] *= dtemp;
				m++;
			}
		}
	}

	// inverse FFT of the product
	InverseFFT();
	GetRealXArray3D(xa3);

	// restore L1 norm
	xa3 *= dnorm1 / xa3.Norm(xar::eNormL1);
}