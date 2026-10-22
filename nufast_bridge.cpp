/*
NuFast-LBL bridge used by the Python benchmark.

One source file builds TWO shared libraries:
    nufast_common.dll
        Hamiltonian parameter set and conventions
    nufast_native.dll
        Original NuFast parameter set and conventions
*/

#include <array>
#include <vector>
#include "Benchmark.h"
#include "NuFast_LBL.h"
using std::array;
using std::vector;

double const eVsqkm_to_GeV_over4 = (1.0e3 / 197.3269804e-9) / (4.0e9);
#if defined(NUFAST_VCC_COMMON) && defined(NUFAST_VCC_NATIVE)
#error "Select only one NuFast matter-conversion mode."
#elif defined(NUFAST_VCC_COMMON)
// VCC_PER_RHOYE = np.sqrt(2.0) * GF_GEV * NA * DENSITY_CONVERSION 
// YERHOE2A = 2.0e9 * VCC_PER_RHOYE
double const YerhoE2a = 1.526493243e-4;

#elif defined(NUFAST_VCC_NATIVE)
double const YerhoE2a = 1.52588e-4;
#else
#error "Define NUFAST_VCC_COMMON or NUFAST_VCC_NATIVE."
#endif
double const oneThird = 1.0 / 3.0;


#if defined(_WIN32)
#define EXPORT extern "C" __declspec(dllexport)
#else
#define EXPORT extern "C"
#endif


EXPORT void nufast_matter(
    double s12sq,
    double s13sq,
    double s23sq,
    double delta,
    double dm21,
    double dm31,
    double L,
    const double* E,
    int nE,
    double rhoYe,
    int N_Newton,
    double* out
)
{
    vector<double> energies(E, E + nE);
    vector<array<array<double, 3>, 3>> probabilities(nE);

    NuFast::Probability_Matter_LBL(
        s12sq,
        s13sq,
        s23sq,
        delta,
        dm21,
        dm31,
        L,
        energies,
        rhoYe,
        N_Newton,
        probabilities
    );

    for (int i = 0; i < nE; ++i)
    {
        for (int alpha = 0; alpha < 3; ++alpha)
        {
            for (int beta = 0; beta < 3; ++beta)
            {
                out[(i * 3 + alpha) * 3 + beta] =
                    probabilities[i][alpha][beta];
            }
        }
    }
}


EXPORT void nufast_vacuum(
    double s12sq,
    double s13sq,
    double s23sq,
    double delta,
    double dm21,
    double dm31,
    double L,
    const double* E,
    int nE,
    double* out
)
{
    vector<double> energies(E, E + nE);
    vector<array<array<double, 3>, 3>> probabilities(nE);

    /*
    Current NuFast vacuum interface accepts two arguments that it does not use, 
    so its input structure stays similar to the matter function.
    */
    double unused_density = 0.0;
    int unused_newton = 0;

    NuFast::Probability_Vacuum_LBL(
        s12sq,
        s13sq,
        s23sq,
        delta,
        dm21,
        dm31,
        L,
        energies,
        unused_density,
        unused_newton,
        probabilities
    );

    for (int i = 0; i < nE; ++i)
    {
        for (int alpha = 0; alpha < 3; ++alpha)
        {
            for (int beta = 0; beta < 3; ++beta)
            {
                out[(i * 3 + alpha) * 3 + beta] =
                    probabilities[i][alpha][beta];
            }
        }
    }
}
