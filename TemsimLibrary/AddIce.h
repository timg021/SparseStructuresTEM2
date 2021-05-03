#pragma once
#include <string>
#include <vector>
#include "pdb.h"

int AddIce(float iceThick, float ctblength, int natom, int** pZnum, float** px, float** py, float** pz, float** pocc, float** pwobble, float wobbleaver, unsigned long* piseed);
int AddCarbon(double cThick, double cWidth, pdbdata& pd, unsigned long* piseed, double& xmin, double& xmax, double& ymin, double& ymax, double& zmin, double& zmax);
