#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <xstring>
#include "pdb.h"


int main(int argc, char* argv[])
{
	char pdbfile[1024], pdbfile1[1024];
	char outfile[1024];
	char cfileinfo[1024];
	char cline[1024], ctitle[1024], cparam[1024]; // auxiliary storage

	// Read input parameter file pdb.txt
	printf("\nStarting pdb-compare program ...\n");

	std::string sInputParamFile("pdb-compare.txt");
	if (argc > 1) sInputParamFile = argv[1];
	FILE* ffpar = fopen(sInputParamFile.c_str(), "rt");
	if (!ffpar)
	{
		printf("\nError: cannot open parameter file %s!\n", sInputParamFile.c_str());
		return -1;
	}
	else printf("\nReading input parameter file %s ...\n", sInputParamFile.c_str());

	// read and skip an arbitrary number of initial comment lines (i.e. the lines that start with // symbols)
	while (true)
	{
		fgets(cline, 1024, ffpar);
		if (!(cline[0] == '/' && cline[1] == '/')) break;
	}

	// line 1
	strtok(cline, "\n"); // 1nd parameter: input test file name
	if (sscanf(cline, "%s %s", ctitle, pdbfile) != 2)
	{
		printf("\n!!!Error reading input test file name from input parameter file.");
		return -1;
	}

	// line 2
	fgets(cline, 1024, ffpar); strtok(cline, "\n"); // 2nd parameter: input test file type
	if (sscanf(cline, "%s %s", ctitle, cparam) != 2)
	{
		printf("\n!!!Error reading input test file type from input parameter file.");
		return -1;
	}
	int nfiletype = 0; // input file type 0 - for PDB input file, 1 - for Vesta XYZ input file, 2 - for Kirkland XYZ file.
	nfiletype = atoi(cparam);
	if (nfiletype != 0 && nfiletype != 1 && nfiletype != 2)
	{
		printf("\n!!!Unknown input test file type in input parameter file.");
		return -1;
	}

	// line 3
	fgets(cline, 1024, ffpar); strtok(cline, "\n"); // 3nd parameter: input reference file name
	if (sscanf(cline, "%s %s", ctitle, pdbfile1) != 2)
	{
		printf("\n!!!Error reading input reference file name from input parameter file.");
		return -1;
	}

	// line 4
	fgets(cline, 1024, ffpar); strtok(cline, "\n"); // 4nd parameter: input reference file type
	if (sscanf(cline, "%s %s", ctitle, cparam) != 2)
	{
		printf("\n!!!Error reading input reference file type from input parameter file.");
		return -1;
	}
	int nfiletype1 = 0; // input file type 0 - for PDB input file, 1 - for Vesta XYZ input file, 2 - for Kirkland XYZ file.
	nfiletype1 = atoi(cparam);
	if (nfiletype1 != 0 && nfiletype1 != 1 && nfiletype1 != 2)
	{
		printf("\n!!!Unknown input reference file type in input parameter file.");
		return -1;
	}

	// line 5
	fgets(cline, 1024, ffpar); strtok(cline, "\n"); // 5th parameter: output file name
	if (sscanf(cline, "%s %s", ctitle, outfile) != 2)
	{
		printf("\n!!!Error reading output file name from input parameter file.");
		return -1;
	}

	// line 8
	fgets(cline, 1024, ffpar); strtok(cline, "\n"); // 8th parameter: free form first line for the output file
	if (sscanf(cline, "%s %s", ctitle, cfileinfo) != 2)
	{
		printf("\n!!!Error reading free-form line for the output file from input parameter file.");
		return -1;
	}

	fclose(ffpar);


	// read input file1
	pdbdata pd;
	pdbdata_init(&pd);
	if (nfiletype == 0)
	{
		printf("Reading input file1 %s in PDB format ...\n", pdbfile);
		if (read_pdb(pdbfile, &pd) == -1) return -1; // read PDB file
	}
	else if (nfiletype == 1)
	{
		printf("Reading input file1 %s in Vesta XYZ export format ...\n", pdbfile);
		if (data_from_VestaXYZfile(pdbfile, &pd) == -1) return -1; // read Vesta export XYZ file
	}
	else if (nfiletype == 2)
	{
		printf("Reading input file1 %s in Kirkland XYZ format ...\n", pdbfile);
		if (data_from_KirklandXYZfile(pdbfile, &pd) == -1) return -1; // read Kirkland XYZ file
	}

	// read input file2
	pdbdata pd1;
	pdbdata_init(&pd1);
	if (nfiletype1 == 1)
	{
		printf("Reading input file2 %s in Vesta XYZ export format ...\n", pdbfile1);
		if (data_from_VestaXYZfile(pdbfile1, &pd1) == -1) return -1; // read Vesta export XYZ file
	}
	else if (nfiletype1 == 0)
	{
		printf("Reading input file2 %s in PDB format ...\n", pdbfile1);
		if (read_pdb(pdbfile1, &pd1) == -1) return -1; // read PDB file
	}
	else if (nfiletype1 == 2)
	{
		printf("Reading input file2 %s in Kirkland XYZ format ...\n", pdbfile1);
		if (data_from_KirklandXYZfile(pdbfile1, &pd1) == -1) return -1; // read Kirkland XYZ file
	}
	

	//translate element symbols into atomic weights
	int* ia = (int*)malloc(pd.natoms * sizeof(int));
	if (nfiletype != 2)
	{
		if (pdb_atomnumbers(&pd, ia))
		{
			printf("\n!!!Error encountered while finding atomic weight for a given element name!!!\n");
			return -1;
		}
	}
	else
		for (int i = 0; i < pd.natoms; i++) ia[i] = pd.adata[i].serial;

	int* ia1 = (int*)malloc(pd1.natoms * sizeof(int));
	if (nfiletype1 != 2)
	{
		if (pdb_atomnumbers(&pd1, ia1))
		{
			printf("\n!!!Error encountered while finding atomic weight for a given element name!!!\n");
			return -1;
		}
	}
	else
		for (int i = 0; i < pd1.natoms; i++) ia1[i] = pd1.adata[i].serial;

	// check that the reference data does not contain hydrogen atoms (matching of hydrogen atoms is considered counter-productive since the test positions will be "wasted" on them)
	for (int i = 0; i < pd1.natoms; i++)
		if (ia1[i] == 1)
		{
			printf("\n!!!Error: hydrogen atom found in the reference data - matching of hydrogen atoms in not supported!!!\n");
			return -1;
		}

	// sort entries by atom weight in descending order
	pdb_bubbleSort1(&pd, ia);
	pdb_bubbleSort1(&pd1, ia1);


	printf("\nComparing atoms positions from %s test file with those in %s reference file ...", pdbfile, pdbfile1);
	printf("\n");
	//	if (pd.natoms > pd1.natoms)
	//	{
	//		printf("\n!!!Error: the test file contains more atom entries than the reference file.!!!\n");
	//		return -1;
	//	}


		// open the output file early in order to start writing some results there as they are obtained
	FILE* ff = fopen(outfile, "wt");
	if (ff == NULL)
	{
		printf("\nERROR: cannot open %s file for writing!!!", outfile);
		return -2;
	}
	fprintf(ff, "For each atom from %s file, locating the closest atom in %s file.\n", pdbfile, pdbfile1);
	fprintf(ff, "Each entry contains the following data:\n");
	fprintf(ff, "1st column contains atomic weights (reference data)\n");
	fprintf(ff, "2nd, 3rd, and 4th columns contain x, y and z atomic coordinates, respectively (reference data)\n");
	fprintf(ff, "5th column contains the L2 distance between atoms in the test and the matched reference data.\n");
	fprintf(ff, "6th column contains atomic weights (matched test data)\n");
	fprintf(ff, "7th, 8th, and 9th columns contain x, y and z atomic coordinates, respectively (matched test data)\n");
	fprintf(ff, "10th column contains original data for the test file\n");
	fprintf(ff, "%s", cfileinfo); // free-form file info line
	fprintf(ff, "\n");

	// clean the entries for the future matched indexes
	for (int j = 0; j < pd.natoms; j++) pd.adata[j].tempFactor = -1; // this parameter will contain the index of the found best match
	for (int i = 0; i < pd1.natoms; i++) pd1.adata[i].tempFactor = -1; // this parameter will contain the index of the found best match

	// for each atom in the reference structure pd1 find the closest atom in the test structure pd, and save the results in the new "double" structure pd2 
	int i0, jmin;
	double r2, r2min;
	for (int i = 0; i < pd1.natoms; i++)
	{
		jmin = 0;
		r2min = (pd1.adata[i].x - pd.adata[0].x) * (pd1.adata[i].x - pd.adata[0].x) + (pd1.adata[i].y - pd.adata[0].y) * (pd1.adata[i].y - pd.adata[0].y) + (pd1.adata[i].z - pd.adata[0].z) * (pd1.adata[i].z - pd.adata[0].z);
		for (int j = 1; j < pd.natoms; j++)
		{
			r2 = (pd1.adata[i].x - pd.adata[j].x) * (pd1.adata[i].x - pd.adata[j].x) + (pd1.adata[i].y - pd.adata[j].y) * (pd1.adata[i].y - pd.adata[j].y) + (pd1.adata[i].z - pd.adata[j].z) * (pd1.adata[i].z - pd.adata[j].z);
			if (r2 < r2min) { r2min = r2;  jmin = j; }
		}

		if (pd.adata[jmin].tempFactor == -1) // this test atom has not been matched yet
		{
			pd.adata[jmin].tempFactor = i;
			pd1.adata[i].tempFactor = jmin;
			pd1.adata[i].occupancy = sqrt(r2min);
		}
		else // this test atom has been already matched with another reference atom before
		{
			// let us find out if the new match is better than the old one, and if so - reassign the previous assignment to the new one
			i0 = (int)pd.adata[jmin].tempFactor;
			r2 = (pd1.adata[i0].x - pd.adata[jmin].x) * (pd1.adata[i0].x - pd.adata[jmin].x) + (pd1.adata[i0].y - pd.adata[jmin].y) * (pd1.adata[i0].y - pd.adata[jmin].y) + (pd1.adata[i0].z - pd.adata[jmin].z) * (pd1.adata[i0].z - pd.adata[jmin].z);
			if (r2min < r2) // new match is better
			{
				pd1.adata[i0].tempFactor = -1; // invalidate the old match
				pd.adata[jmin].tempFactor = i;
				pd1.adata[i].tempFactor = jmin;
				pd1.adata[i].occupancy = sqrt(r2min);
			}
			// else do nothing, i.e. ignore this newly found match, since it is worse than the old one
		}
	}

	// calculate average distance and std
	int nodup = 0;
	double adist = 0.0, maxdist = 0.0;
	for (int i = 0; i < pd1.natoms; i++)
	{
		if (pd1.adata[i].tempFactor == -1) continue;
		nodup++;
		adist += pd1.adata[i].occupancy;
		if (pd1.adata[i].occupancy > maxdist) maxdist = pd1.adata[i].occupancy;
	}
	if (nodup == 0) adist = 0.0; else adist /= (double)nodup;
	double stddist = 0.0;
	for (int i = 0; i < pd1.natoms; i++)
	{
		if (pd1.adata[i].tempFactor == -1) continue;
		stddist += (pd1.adata[i].occupancy - adist) * (pd1.adata[i].occupancy - adist);
	}
	if (nodup <= 1) stddist = 0.0; else stddist = sqrt(stddist / ((double)(nodup - 1)));
	printf("\n%d out of total %d atoms in the reference structure (i.e. %g percent) have been uniquely matched with atoms of the test structure.", nodup, pd1.natoms, (double)nodup / (double)pd1.natoms * 100.0);
	printf("\nThe result contains %d false negatives (%g percent) and %d false positives (%g percent).", pd1.natoms - nodup, 100.0 * (double)(pd1.natoms - nodup) / (double)pd1.natoms, pd.natoms - nodup, 100.0 * (double)(pd.natoms - nodup) / (double)pd.natoms);
	printf("\nAverage distance between the matched test and template atoms = %g.", adist);
	printf("\nMaximum distance between the matched test and template atoms = %g.", maxdist);
	printf("\nStandard deviation of the distance between the matched test and template atoms = %g.", stddist);
	printf("\n");
	fprintf(ff, "\n%d atoms (out of total %d atoms) from the reference structure (i.e. %g percent) have been uniquely matched with atoms of the test structure.", nodup, pd1.natoms, (double)nodup / (double)pd1.natoms * 100.0);
	fprintf(ff, "\nThe result contains %d false negatives (%g percent) and %d false positives (%g percent).", pd1.natoms - nodup, 100.0 * (double)(pd1.natoms - nodup) / (double)pd1.natoms, pd.natoms - nodup, 100.0 * (double)(pd.natoms - nodup) / (double)pd.natoms);
	fprintf(ff, "\nAverage distance between the matched test and template atoms = %g.", adist);
	fprintf(ff, "\nMaximum distance between the matched test and template atoms = %g.", maxdist);
	fprintf(ff, "\nStandard deviation of the distance between the matched test and template atoms = %g.", stddist);
	fprintf(ff, "\n");

	// list all atoms from the reference file which have not been matched by any of the atoms from the test file
	for (int i = 0; i < pd1.natoms; i++)
		if (pd1.adata[i].tempFactor == -1)
		{
			//printf("\n@@@ Reference atom no. %d, atom type = %d, atom positions = (%g, %g, %g) has not been matched.", i, ia1[i], pd1.adata[i].x, pd1.adata[i].y, pd1.adata[i].z);
			fprintf(ff, "\n@@@ Reference atom no. %d, atom type = %d, atom positions = (%g, %g, %g) has not been matched.", i, ia1[i], pd1.adata[i].x, pd1.adata[i].y, pd1.adata[i].z);
		}

	// output to the target file
	printf("\nWriting output comparison file %s ...\n", outfile);
	fprintf(ff, "\n");
	for (int i = 0; i < pd1.natoms; i++)
	{
		if (pd1.adata[i].tempFactor == -1) continue;
		jmin = (int)pd1.adata[i].tempFactor;
		fprintf(ff, "\n%d %f %f %f %f", ia1[i], pd1.adata[i].x, pd1.adata[i].y, pd1.adata[i].z, pd1.adata[i].occupancy);
		fprintf(ff, " %d %f %f %f %f ", ia[jmin], pd.adata[jmin].x, pd.adata[jmin].y, pd.adata[jmin].z, pd.adata[jmin].occupancy);
	}
	fclose(ff);
	free(ia1);
	free(ia);

	printf("Finished!");

	return 0;
}