Brief description of the PDB-compare module "Pdb-compare.exe"
of the Differential Holographic Tomography (DHT) software package

This program is based on A.V. Martin's PDB code.

Pdb-compare.exe module compares PDB and XYZ atomic structure files, producing the output data 
in ASCII text format (see details below). In particular, this program can be used to assess
the accuracy of 3D reconstruction of atomic structures (e.g. performed by PhaseRetrieval program)
by comparing the reconstructed XYZ file with a known XYZ file of the same strcture. 

This program has been written in C++ and can in principle be 
compiled under any OS supporting standard C++ compilers and execution environment. 
However, at present the executable module is only available for 64-bit Windows OS. Any 
Windows x64 PC can in principle be used for running this code.

The control of execution of Pdb-compare.exe program is managed via an editable text file 
Pdb-compare.txt which contains all modifiable input parameters for Pdb-compare.exe. In 
order for the program to run, the input parameter file Pdb-compare.txt must be present in the 
same folder where Pdb-compare.exe is started from. Alternatively, the name of a suitable input
parameter file can be presented as the first command-line argument (following the pdb-compare.exe
name in the command line) when starting pdb-compare.exe program. In this case, the parameter
file can be given an arbitrary name and can be located in any place accessible from the started
instance of pdb-compare.exe. The format of these input parameter files is fixed and must
comply exactly with the structure described below.

The format of Pdb-compare.txt file allows any number of "comment" lines to be present at the beginning and/or
at the end of the file, each such line starting with the double forward slash symbol, "//". The comment lines 
are ignored by Pdb-compare.exe program. Other (non-comment) lines must all have the fixed structure 
(described below) and the fixed sequential order with respect to each other. Each non-comment line consists of:
(a) An arbitrary "expression" with no white spaces, which typically represents the "name" (brief description) of a 
parameter. The contents of this parameter name expression are ignored by Pdb-compare.exe. 
An example of the parameter name is "2.Type_of_the_input_test_file_(see_below):".
(b) One or more alphanumerical entries separated by a single white space from each other 
and from the parameter name, each such entry representing one parameter value, for example 
" 2", which corresponds to one parameter value equal to 2. Each line should be terminated by the usual "end-
of-line / new line / caret-return", etc., symbols, as common for a given OS. The number and 
types of parameters in each parameter line of the file Pdb-compare.txt is pre-determined and 
cannot be changed, only the values of these parameters can be edited. The parameter names 
and comment lines are supposed to give explicit hints about the type and number of parameters expected
in each line. 

Here is a typical example of a valid Pdb-compare.txt parameter file.

=======================================================================================
*** PDB compare program
1.Input_test_file_name: outAmg2000.xyz
2.Type_of_the_input_test_file_(see_below): 2
3.Rotate_test_structure_around_z,_y'_and_z"_axes_by_(degrees): 0 0 0
4.Input_reference_file_name: 7kod_DW_noH.xyz
5.Type_of_the_input_reference_file_(see_below): 2
6.Output_file_name: zzzmg2000_compare.txt
7.Free-form_line_for_output_file: ***

//Input parameter file for pdb.exe program
//1st line contains input test file name
//2nd line: type of the 1st input file: 0 - PDB input file, 1 - Vesta XYZ input file, 2 - Kirkland XYZ file.
//3rd line contains Euler rotation angles in degrees that are applied to the test structure before comparison
//4th line contains input reference file name - must be equal to or longer than the test file
//5th line: type of the 2nd input file: 0 - PDB input file, 1 - Vesta XYZ input file, 2 - Kirkland XYZ file.
//6th line contains output file name
//7th line contains the free-form info line that will be recorded into the first line of the output file
========================================================================================

Note that the "numbering" of parameters in the parameter names, such as "2" in 
"2.Type_of_the_input_test_file_(see_below): 2", is inconsequential and is ignored by 
Pdb-compare.exe. On the other hand, the order of the parameter lines in this file is very 
important and must be the same as in the above example. In other words, Pdb-compare.exe 
interprets the input parameters according to the order of the non-comment lines in the 
parameter file, while ignoring any numbering that may or may not appear in the parameter 
names.

In the above example, Parameter 1 contains the name of a "test" atomic structure file in 
PDB (https://en.wikipedia.org/wiki/Protein_Data_Bank_(file_format), Vesta XYZ (https://jp-minerals.org/vesta/en/) or
Kirkland XYZ (35. E. J. Kirkland, Advanced Computing in Electron Microscopy, second edition, Springer, New York, 2010) format. 
This input file must be present in the same folder where Pdb-compare.exe is started from, or, 
alternatively, the filename can include a fully specified pathname (OS specific). 

Parameter 2 can be equal to 0, 1 or 2. It defines the format of the first input file.
"0" corresponds to PDB input file format. 
"1" corresponds to Vesta XYZ input file format.
"2" corresponds to Kirkland XYZ input file format.

Parameter 3 contains Euler rotation angles (around z, y and z" axes, in this order) in degrees that are applied
to the test structure before comparison with the reference structure. This can be helpful in checking known symmetries
of the structures.

Parameter 4 contains the name of a "reference" atomic structure file in PDB , Vesta XYZ or
Kirkland XYZ format. This input file must be present in the same folder where Pdb-compare.exe is started from, or, 
alternatively, the filename can include a fully specified pathname (OS specific). 

Parameter 5 can be equal to 0, 1 or 2. It defines the format of the second input file.
"0" corresponds to PDB input file format. 
"1" corresponds to Vesta XYZ input file format.
"2" corresponds to Kirkland XYZ input file format.

Parameter 6 contains the name of the output file. This file will be saved in ASCII text format.

Parameter 7 contains a free-form information line (with no white spaces) that is written 
as the first line in the output file.

================================================================================================== 
Example of a valid XYZ file in Kirkland's XYZ format (see detailed description in the 
references given above):

Aspartate_molecule
10.000000 10.000000 10.000000
7 4.386000 6.358500 5.031500 1.000000 0.1
6 4.233000 4.956500 4.621500 1.000000 0.1
6 2.835000 4.490500 4.936500 1.000000 0.1
8 2.169000 5.085500 5.751500 1.000000 0.1
6 5.242000 4.090500 5.378500 1.000000 0.1
6 6.641000 4.475500 4.969500 1.000000 0.1
8 6.812000 5.351500 4.155500 1.000000 0.1
8 7.695000 3.844500 5.508500 1.000000 0.1
8 2.329000 3.414500 4.313500 1.000000 0.1
1 3.775000 6.959500 4.498500 1.000000 0.1
1 4.225000 6.465500 6.021500 1.000000 0.1
1 4.411000 4.869500 3.549500 1.000000 0.1
1 5.122000 4.245500 6.450500 1.000000 0.1
1 5.070000 3.040500 5.141500 1.000000 0.1
1 8.572000 4.125500 5.215500 1.000000 0.1
1 1.428000 3.153500 4.549500 1.000000 0.1
-1

=========================================================================================================
Example of a valid XYZ file in Vesta's XYZ format (see detailed description in the 
references given above):

16
Aspartate_molecule
O 2.169000 6.446000 5.869000 0.000000
O 6.812000 6.712000 4.273000 0.000000
O 7.695000 5.205000 5.626000 0.000000
O 2.329000 4.775000 4.431000 0.000000
N 4.386000 7.719000 5.149000 0.000000
C 4.233000 6.317000 4.739000 0.000000
C 2.835000 5.851000 5.054000 0.000000
C 5.242000 5.451000 5.496000 0.000000
C 6.641000 5.836000 5.087000 0.000000
H 3.775000 8.320000 4.616000 0.000000
H 4.225000 7.826000 6.139000 0.000000
H 4.411000 6.230000 3.667000 0.000000
H 5.122000 5.606000 6.568000 0.000000
H 5.070000 4.401000 5.259000 0.000000
H 8.572000 5.486000 5.333000 0.000000
H 1.428000 4.514000 4.667000 0.000000

===========================================================================================================
Example of a valid file in muSTEM's XTL format (see detailed description in the 
references given above):

Aspartate_molecule
10.000000 10.000000 10.000000 90.0 90.0 90.0
200.000000 
4 
O
4 8.000000 1.0 0.0
0.216900 0.644600 0.586900
0.681200 0.671200 0.427300
0.769500 0.520500 0.562600
0.232900 0.477500 0.443100
N
1 7.000000 1.0 0.0
0.438600 0.771900 0.514900
C
4 6.000000 1.0 0.0
0.423300 0.631700 0.473900
0.283500 0.585100 0.505400
0.524200 0.545100 0.549600
0.664100 0.583600 0.508700
H
7 1.000000 1.0 0.0
0.377500 0.832000 0.461600
0.422500 0.782600 0.613900
0.441100 0.623000 0.366700
0.512200 0.560600 0.656800
0.507000 0.440100 0.525900
0.857200 0.548600 0.533300
0.142800 0.451400 0.466700

============================================================================================================