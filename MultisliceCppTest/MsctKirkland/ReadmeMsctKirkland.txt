Brief description of the "forward" simulation module "MsctKirkland.exe"
of the Differential Holographic Tomography (DHT) software package

This program is based on E.J. Kirkland's temsim code which is described in the following 
references [https://github.com/jhgorse/kirkland/tree/master/temsim. 
See also https://sourceforge.net/projects/computem/. 
Description of the code and the underlying theory can be found in E. J. Kirkland, Advanced 
Computing in Electron Microscopy, second edition, Springer, New York, 2010.] 

MsctKirkland.exe module performs "forward" simulation tasks, producing the data (images) 
that can be subsequently used as input for the DHT reconstruction or for comparison with the 
DHT reconstruction results. This program has been written in C++ and can in principle be 
compiled under any OS supporting standard C++ compilers and execution environment. 
However, at present the executable module is only available for 64-bit Windows OS. Any 
Windows x64 PC can in principle be used for running this code, but multi-core CPU systems 
are recommended for faster execution.

The two main "modes" of execution of MsctKirkland are:  
(1) generation of a set of 2D defocused images (distributions of the intensity, phase or 
complex amplitude) at a set of rotational orientations of the specified 
atomic structure and a set of defocus distances at each rotational position; or  
(2) generation of a 3D distribution of the electrostatic potential, due to the atoms in the 
specified atomic structure, in an arbitrary volume of space.

The control of execution of MsctKirkland.exe program is managed via an editable text file 
MsctKirkland.txt which contains all modifiable input parameters for MsctKirkland.exe. In 
order for the program to run, the input parameter file MsctKirkland.txt must be present in the 
same folder where MsctKirkland.exe is started from. Alternatively, the name of a suitable input
parameter file can be presented as the first command-line argument (following the MsctKirkland.exe
name in the command line) when starting MsctKirkland.exe program. In this case, the parameter
file can be given an arbitrary name and can be located in any place accessible from the started
instance of MsctKirkland.exe. The format of these input parameter files is fixed and must
comply exactly with the structure described below.

The format of MsctKirkland.txt file allows any number of "comment" lines to be present at the beginning and/or
at the end of the file, each such line starting with the double forward slash symbol, "//". The comment lines 
are ignored by MsctKirkland.exe program. Other (non-comment) lines must all have the fixed structure 
(described below) and the fixed sequential order with respect to each other. Each non-comment line consists of:
(a) An arbitrary "expression" with no white spaces, which typically represents the "name" (brief description) of a 
parameter. The contents of this parameter name expression are ignored by MsctKirkland.exe. 
An example of the parameter name is "6.Wavefunction_size_in_pixels,_Nx,Ny:".
(b) One or more alphanumerical entries separated by a single white space from each other 
and from the parameter name, each such entry representing one parameter value, for example 
" 256 256", which corresponds to two parameter values equal to 256 each.
Note that the number of white spaces separating the parameters must be equal to one, 
otherwise a "zero" parameter may be wrongly read in. There should be also no further white 
spaces after the last parameter in the line. Each line should be terminated by the usual "end-
of-line / new line / caret-return", etc., symbols, as common for a given OS. The number and 
types of parameters in each parameter line of the file MsctKirkland.txt is pre-determined and 
cannot be changed, only the values of these parameters can be edited. The parameter names 
and comment lines are supposed to give explicit hints about the type and number of parameters expected
in each line. 

The Euler angle convention corresponds to that in Heymann et al, Journal of Structural Biology 151 (2005) 196–207.

Here is a typical example of a valid MsctKirkland.txt parameter file.

//*****Input parameter file for MsctKirkland.exe program*****
1.Input_file_with_atomic_numbers_and_coordinates_in_XYZ_format: ASP_NEW_Kirck_DW.xyz
2.Output_files_shall_contain_intensity(0),_phase(1),_complex_amplitude(2)_or_3D_potential(3): 0
3.Output_TIFF/GRD/GRC_filename_template: C:\Users\tgureyev\Downloads\Temp\asp.tif
4.Use_multislice(0),_projection(1),_or_1st_Born(2)_approximation: 0
5.Incident__electron_beam_energy_in_keV: 200.0
6.Wavefunction_size_in_pixels,_Nx,Ny: 256 256
7.Slice_thickness_in_Angstroms: 1.25
8.Objective_aperture_in_mrad: 45.0
9.Spherical_aberration_Cs3_and_Cs5_in_mm: 0 0
10.Include_thermal_vibrations(1)_or_not(0): 1
11.____Temperature_in_degrees_K: 300.0
12.____Number_of_configurations_to_average_over: 64
13.Ice_layer_thickness_in_Angstroms: 100
14.Text_file_with_output_rotation_angles_in_degrees_and_defocus_distances_in_Angstroms: DefocusRand36_1_NEW.txt
15.Number_of_worker_threads_to_launch: 20
//*****
//*****!!!Don't forget to save this file after editing the parameters

Note that the "numbering" of parameters in the parameter names, such as "6" in 
"6.Wavefunction_size_in_pixels,_Nx,Ny: ", is inconsequential and is ignored by 
MsctKirkland.exe. On the other hand, the order of the parameter lines in this file is very 
important and must be the same as in the above example. In other words, MsctKirkland.exe 
interprets the input parameters according to the order of the non-comment lines in the 
parameter file, while ignoring any numbering that may or may not appear in the parameter 
names.

In the above example, Parameter 1 contains the name of an "atomic structure" file in 
Kirkland's XYZ format (see the references above for the description of this format). The 
XYZ file must be present in the same folder where MsctKirkland.exe is started from, or, 
alternatively, the filename can include a fully specified pathname (OS specific). 

Parameter 2 can be equal to 0, 1, 2 or 3. It defines what type of output data is generated and 
saved by MsctKirkland.exe. "0" corresponds to the output of defocused intensities,
"1" corresponds to the output of 2D phase distributions in the defocus planes, 
"2" corresponds to the output of complex defocused amplitudes. 
"3" corresponds to the direct generation and output of the 3D distribution of the electrostatic potential
(in volts) within the planes defined in the text file given in Parameter 14 (usually, these are equispaced 
planes inside the object). This last option uses Kirkland's function "vatom()" to calculate 
the contribution of the potential from each atom given in the file in Parameter 1 to each
spatial position on the 3D grid specified by Parameters 1, 6 and 14.
Note that when Parameter 2 is equal to "3" (the case of 3D potential output), the format of the
text file with rotation angles and defocused distances in Parameter 14 must be somewhat special (see
an example below), defining the positions of transverse "defocus" planes along a fixed "illumination direction",
so that the set of such planes would fill a 3D volume inside which the 3D potential needs to be evaluated.
The value of Parameter 2 should be compatible with the value of Parameter 3 below. 

Parameter 3 contains a fully specified pathname template for the output files containing either 
(a) defocused images (in GRD format or uncompressed 32-bit floating-point TIFF format), 
(b) defocused phases (in GRD format or uncompressed 32-bit floating-point TIFF format), 
(c) defocused complex amplitudes (in GRC format), or 
(d) 2D cross-sections (in GRD format or uncompressed 32-bit floating-point TIFF format) of the 
3D distribution of the electrostatic potential corresponding to the atomic structure specified
in Parameter 1. 
The GRD and GRC file formats are described below.
Note that here the pathname is actually supposed to contain a "template" for the
eventual output filenames which MsctKirkland.exe creates by inserting one or a pair of numbers
that correspond to the index of a defocus distance (first inserted number) and the index of a 
rotational orientation (second inserted number). As a result, when a filename template 
"ag.grd" is given in this parameter, the output filenames may look like "ag0_00.grd, 
ag0_01.grd, ..., ag0_35.grd, ag1_00.grd, ..., ag1_35.grd, ag2_00.grd, ..., ag2_35.grd" (36 x 3 
= 108 files in total), if 36 different rotational positions with 3 defocus distances at each 
rotational position are specified in a file whose name is given in Parameter 14. 
See an example DefocusRand36_NEW.txt file below.
Note that in the case of GRC file output (i.e. in the case of output complex amplitudes), only 
a single defocused distance is allowed per each illumination direction. This is dictated by the
fact that when such files are used as input for the 3D object reconstruction module, "PhaseRetrieval.exe" 
(see its description in the ReadmePhaseRetrieval.txt file), the complex amplitude is assumed to
contain the correct complex wavefunction at a given defocused distance, which in principle can be propagated
to any other defocused distance by calculating Fresnel diffraction integrals. In particular,
no phase retrieval is required in this case. Obviously, such GRC (complex amplitude) output is intended
for testing and debugging purposes only, since only defocused intensities are expected to be available
under experimentally-relevant conditions. In the latter case, multiple intensity files can be
generated by this program at each illumination direction and later used by PhaseRetrieval.exe module
to perform 2D phase retrieval at each illumination direction.

Parameter 4 can be equal to 0, 1, 2 or 3. This parameter defines what method is used to calculate the 
propagation of the incident plane electron wave with unit intensity through the object defined 
by Parameter 1. "0" corresponds to the multislice approximation (both the Fresnel diffraction 
and multiple scattering are taken into account). "1" corresponds to the projection (line 
integrals) approximation (multiple scattering is taken into account, but the Fresnel diffraction 
inside the object is ignored). "2" corresponds to the 1st Born approximation (multiple 
scattering between different atoms is ignored, but the Fresnel diffraction inside the object is 
taken into account). 

Parameter 5 contains the incident electron beam energy in keV, e.g. 200.0.

Parameter 6 defines the number of numerical grid points along the x and y directions in the 
transverse (xy) planes (which are orthogonal to the optical axis z), e.g. " 256 256". These 
numbers must be integer powers of 2.

Parameter 7 defines the thickness (in angstroms) of the slices in multislice calculations, e.g. 
1.25.

Parameter 8 defines the objective aperture in milliradians. Zero or negative value here is 
interpreted as an infinite aperture.

Parameter 9 contains the values of the spherical aberrations Cs3 and Cs5 (in millimetres) of 
the imaging system. Zero values correspond to the absence of aberrations.

Parameter 10 can be equal to 0 or 1. "0" corresponds to ignoring the effect of the thermal 
vibrations of atoms in the object structure. This case is physically unrealistic and is intended
to be used only for testing purposes. "1" corresponds to taking the thermal vibrations 
into account while calculating the defocused images. This is achieved by repeatedly calculating
the defocused images (if Parameter 2 is equal to 0) or the 3D electrostatic potential 
(if Parameter 2 is equal to 3) for multiple spatial configurations with 
pseudo-random shifts of positions of the atoms present in the XYZ file in Parameter 1, 
with the magnitude of the shifts depending on other parameters (see below), and averaging over
these configurations. Note that the computational time increases proportionally to the number
of these configurations. The number of different spatial configurations to average over is defined 
by Parameter 12 below. If Parameter 12 is equal to 1, when Parameter 10 is also equal to 1, then,
instead of averaging over different spatial configurations of the atoms, the effect of thermal
vibrations is calculated according to a simple "negative exponent" model of the Debye-Waller factor,
according to the theory found e.g. in C.J. Humphreys, Rep.Prog.Phys., vol.42, 1979, pp.1827-1887
(see eqs.(4.4)-(4.6) there). This approach is implemented both for the defocused images (when
Parameter 2 is equal to 0) and for generation of the 3D electrostatic potential (Parameter 2 is 
equal to 3). This method is much faster, compared to averaging over multiple spatial configurations,
but it can be less accurate.

Parameter 11 contains the temperature (in degrees of Kelvin), e.g. 290, which is used for the 
thermal vibration calculations in Kirkland's code. This parameter is used to scale the standard 
deviation of thermal vibrations of different atoms. These standard deviation values (corresponding
to the temperature of 300 K) must be present in column 6 of the XYZ file specified in Parameter 1.
The value of 300 in Parameter 11 corresponds to the scaling equal to 1.

Parameter 12 contains the number of configurations to average over, e.g. 64, while 
calculating the effect of thermal vibrations of atoms on the defocused images or on the 3D 
electrostatic potential. When this parameter is equal to 1, and the thermal vibrations are
enabled in Parameter 10, an alternative method based on average Debye-Waller factor, is used
for calculation of the effect of thermal vibrations on the defocused images or the generated
3D electrostatic potential (see the description of Parameter 10 above for more details).

Parameter 13 contains the thickness (in angstroms) of the optional layer of amorphous ice to
be added to the atoms present in the input XYZ file specified in Parameter 1. If the ice
thickness specified here is larger than zero, H20 molecules will be added at random locations
around the atoms of the structure given in the XYZ file, according to the following constraints.
(a) All H and O atoms will be added within the parallelepiped with the X and Y dimensions equal to
the X and Y dimensions of the structure in the XYZ file (i.e. the first two parameters in line 2
of the XYZ file) and the Z dimension equal to the specified ice layer thickness. Note that the
specified non-zero ice layer thickness must be equal to or larger than the Z dimension of the
structure in the XYZ file (i.e. larger than the third parameter in line 2 of the XYZ file).
(b) All added H2O molecules will be separated by at least 1.4 angstroms from each other and from
all atoms in the XYZ file.
(c) The average density of the ice will be equal to 0.9167 g/cm^3.
This program will also write out the modified XYZ files with added ice for each rotational 
position of the structure from the input XYZ file specified in Parameter 1. The locations 
and names of these files will be similar to those of the output image files (see description
of Parameter 3 above), but the output XYZ file will have the ".xyz" extension, and only one
XYZ file will be written per each illumination direction. Note that in these output XYZ files,
the input XYZ structure will be present in its original rotational orientation with respect to
the illumination axis, but rotated according to the first two rotation angles specified in the
input file in Parameter 14 (see the description below). Note also that the ice is added only
once per illumination direction, so all images at different defocus distances at this direction
see the same ice.

Parameter 14 contains the name of a text file with the object rotation angles in degrees and 
defocus distances in angstroms for the calculation of the defocused intensities or complex 
amplitudes. An example of such a file is given below. This file may have one of two distinct
formats, which is determined by the file extension. 
(a) If the file extenstion is ".txt", then this file may contain an arbitrary number of leading
comment lines, each starting with a pair of forward slashes "//". Non-comment lines all have the
same structure, but may contain different number of entries, according to the following scheme. 
Columns one and two always contain the rotation angles (in degrees) around the Z and Y' axes, 
respectively (where Y' axis is the Y axis after the initial rotation of the 3D space around the Z axis).
These two angles are followed by an arbitrary number of pairs of values, with each odd column 
(starting from the third one) containing the angle (in degrees) of rotation around the Z" axis
(Z" axis is the Z axis  after the initial rotations of the 3D space around the Z and Y' axes,
which corresponds to the optic axis = the illumination direction). 
Each even column (starting from the fourth one) contains a defocus distance in angstroms along the Z
axis. The number of lines in this file is not limited. Note that this file
has a somewhat "special" form in the case of output of the 3D electrostatic potential. A suitable example 
is also given below.
(b) If the file extension is ".RELION", then this file contains an arbitrary number of lines,
each line containing exactly nine entries, separated by white spaces, with the following contents:
1st entry contains the sequential number of the line - this parameter is not used in MsctKirkland.exe
2nd entry contains the name of the file where the data was sourced from - this parameter is not used in MsctKirkland.exe
3rd entry is the defocus distance dx corresponding to the x coordinate
4th entry is the defocus distance dy corresponding to the y coordinate
5th entry is the image shift along x in angstroms 
6th entry is the image shift along y in angstroms 
7th entry is the rotation angle "rot" around the Z coordinate in radians
8th entry is the rotation angle "tilt" around the Y' coordinate in radians
9th entry is the rotation angle "psi" around the Z" coordinate in radians
A suitable example is given below.
Note that in the case of .RELION format file, the defocus distances given in this file are considered
to be measured from the centre of the molecule, while in the case of .TXT format files the defocus
distances are measured from the "exit" plane of the cube containing the molecule. In order to reconcile
this, the defocus distances given in .RELION file are adjusted by subtracting one-half of the 
z-extent of the molecule as found in the input XYZ file given in Parameter 1.
The text file given in Parameter 14 must be present in the same folder where MsctKirkland.exe is started from,
or, alternatively, the filename can include a fully specified pathname (OS specific). 

Parameter 15 contains the desired number of worker threads to launch during the execution of 
MsctKirkland.exe. The recommend number is equal to the hyper-threading capacity of one's 
computer (often, this number is equal to the number of available CPU cores times 2), minus 1 
(one thread may be reserved for the "main user thread" that leaves the computer responsive to 
user interactions during prolonged calculations). For example, if you run MsctKirkland.exe 
on a PC with an Intel Core i5 CPU with 6 cores, the recommended value of this parameter is 11.

================================================================================================== 
Example of a valid XYZ file in Kirkland's XYZ format (see detailed description in the 
references given above):

Aspartate molecule from PDB site via Vesta in xyz format
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

======================================================================================= 
** Example of a text file with ".txt" extension whose name may appear in Parameter 14 when generating 
defocused images, phases or complex amplitudes:

//Z_rotation	Y'_rotation	Z''_rotation1	Defocus1	Z''_rotation2	Defocus2
//phi=360.0*rand()	theta=180/PI()*ACOS(2*RAND()-1)	psi=360.0*rand()	dz=5 + 10 * RAND()	psi=360.0*rand()	dz=5 + 10 * RAND()
134.545593	20.8879562	250.2507883	8.222943558	100.118147	5.766113138
328.9524395	68.27488858	343.0435375	14.69256862	155.8774949	13.27846134
341.2124109	58.14192331	259.6780077	9.394382628	159.9422219	7.886505485
196.5382573	35.93093368	155.7524815	10.66208447	132.2968222	10.23461056
259.1734396	92.41853487	191.961092	13.05924461	21.64331534	14.06469824
100.7594369	63.88481208	291.9033034	8.891866815	289.5271391	12.0210506
323.1600184	167.2336091	276.145737	5.124502934	154.6881345	12.70999187
259.7074261	45.77533335	183.20903	12.21241235	302.7163076	6.172852555
207.5392464	68.30142143	199.7200929	9.077328511	238.9836551	11.8385373
171.1598766	89.68242302	92.48191277	6.655790719	347.0470774	11.80187741
60.55292977	25.93491786	211.3969304	7.568231724	226.5960793	5.605451246
48.85497938	31.68942623	256.624174	11.403842	87.85325344	6.75988388
357.5837845	132.619438	349.2236099	11.97751674	289.0335287	13.37860964
162.3045554	79.02550822	257.9751296	9.164385773	326.2744611	7.944083188
266.5438071	119.8127464	208.8295398	7.316264206	198.8033863	5.846068438
347.4672745	43.82225966	170.2100489	8.462136081	257.2294484	6.048567208

======================================================================================= 
** Example of a text file with ".RELION" extenstion whose name may appear in Parameter 14 when generating 
defocused images, phases or complex amplitudes:
0 b'J19/extract/16300063732981320253_FoilHole_24015405_Data_24016401_24016403_20200225_0014_Fractions.mrc_rigid_aligned_particles.mrc' 13067.9072265625 12987.794921875 -3.5859375 4.3828125 -2.7313623428344727 -1.9179142713546753 -0.15076839923858643 
1 b'J19/extract/16300063732981320253_FoilHole_24015405_Data_24016401_24016403_20200225_0014_Fractions.mrc_rigid_aligned_particles.mrc' 13067.9072265625 12987.794921875 -0.557812511920929 -8.685937881469727 0.5294424891471863 -0.038568660616874695 -3.0679616928100586 
2 b'J19/extract/16300063732981320253_FoilHole_24015405_Data_24016401_24016403_20200225_0014_Fractions.mrc_rigid_aligned_particles.mrc' 13067.9072265625 12987.794921875 -0.557812511920929 -5.498437404632568 -3.285348653793335 -0.859029233455658 0.08765604346990585 
3 b'J19/extract/16300063732981320253_FoilHole_24015405_Data_24016401_24016403_20200225_0014_Fractions.mrc_rigid_aligned_particles.mrc' 13067.9072265625 12987.794921875 -3.9046874046325684 2.629687547683716 1.1114786863327026 1.8127269744873047 -1.9529767036437988 
4 b'J19/extract/16300063732981320253_FoilHole_24015405_Data_24016401_24016403_20200225_0014_Fractions.mrc_rigid_aligned_particles.mrc' 13067.9072265625 12987.794921875 0.3984375 -5.339062690734863 -2.198413610458374 -1.8407769203186035 0.19284330308437347 
5 b'J19/extract/16300063732981320253_FoilHole_24015405_Data_24016401_24016403_20200225_0014_Fractions.mrc_rigid_aligned_particles.mrc' 13067.9072265625 12987.794921875 2.629687547683716 -6.454687595367432 0.14375591278076172 3.215223789215088 0.34711793065071106 
6 b'J19/extract/16300063732981320253_FoilHole_24015405_Data_24016401_24016403_20200225_0014_Fractions.mrc_rigid_aligned_particles.mrc' 13067.9072265625 12987.794921875 7.251562595367432 -3.426562547683716 -1.076416254043579 -2.815512180328369 0.5785298943519592 
7 b'J19/extract/16300063732981320253_FoilHole_24015405_Data_24016401_24016403_20200225_0014_Fractions.mrc_rigid_aligned_particles.mrc' 13067.9072265625 12987.794921875 1.514062523841858 -6.614062309265137 1.8618143796920776 -0.9221416115760803 -1.4971652030944824 
8 b'J19/extract/16300063732981320253_FoilHole_24015405_Data_24016401_24016403_20200225_0014_Fractions.mrc_rigid_aligned_particles.mrc' 13067.9072265625 12987.794921875 -6.295312404632568 0.23906250298023224 1.055378794670105 3.2993736267089844 0.01753120869398117 
9 b'J19/extract/16300063732981320253_FoilHole_24015405_Data_24016401_24016403_20200225_0014_Fractions.mrc_rigid_aligned_particles.mrc' 13067.9072265625 12987.794921875 2.151562452316284 1.9921875 1.6374149322509766 1.9670016765594482 2.1843886375427246 
10 b'J19/extract/16300063732981320253_FoilHole_24015405_Data_24016401_24016403_20200225_0014_Fractions.mrc_rigid_aligned_particles.mrc' 13067.9072265625 12987.794921875 7.410937309265137 -5.498437404632568 -2.4228131771087646 -2.394763231277466 -0.6626796722412109 
11 b'J19/extract/16300063732981320253_FoilHole_24015405_Data_24016401_24016403_20200225_0014_Fractions.mrc_rigid_aligned_particles.mrc' 13067.9072265625 12987.794921875 -0.7171875238418579 -8.048437118530273 -1.8337644338607788 -2.7874622344970703 -0.5855423808097839 
12 b'J19/extract/16300063732981320253_FoilHole_24015405_Data_24016401_24016403_20200225_0014_Fractions.mrc_rigid_aligned_particles.mrc' 13067.9072265625 12987.794921875 -4.3828125 6.135937690734863 0.859029233455658 -3.2082111835479736 -0.5434674620628357 
13 b'J19/extract/16300063732981320253_FoilHole_24015405_Data_24016401_24016403_20200225_0014_Fractions.mrc_rigid_aligned_particles.mrc' 13067.9072265625 12987.794921875 3.1078124046325684 11.395312309265137 -2.577087640762329 0.2980305552482605 0.17881833016872406 
14 b'J19/extract/16300063732981320253_FoilHole_24015405_Data_24016401_24016403_20200225_0014_Fractions.mrc_rigid_aligned_particles.mrc' 13067.9072265625 12987.794921875 -2.4703125953674316 3.5859375 2.773437261581421 0.10869349539279938 1.6654648780822754 
15 b'J19/extract/16300063732981320253_FoilHole_24015405_Data_24016401_24016403_20200225_0014_Fractions.mrc_rigid_aligned_particles.mrc' 13067.9072265625 12987.794921875 1.673437476158142 6.932812690734863 2.5280003547668457 0.8379917740821838 -1.8337644338607788 
16 b'J19/extract/16300063732981320253_FoilHole_24015405_Data_24016401_24016403_20200225_0014_Fractions.mrc_rigid_aligned_particles.mrc' 13067.9072265625 12987.794921875 -5.020312309265137 0.3984375 -0.6346297860145569 -3.1731488704681396 -0.49438008666038513 
17 b'J19/extract/16300063732981320253_FoilHole_24015405_Data_24016401_24016403_20200225_0014_Fractions.mrc_rigid_aligned_particles.mrc' 13067.9072265625 12987.794921875 8.048437118530273 7.410937309265137 3.0048491954803467 0.5785298943519592 -0.8520167469978333 
18 b'J19/extract/16300063732981320253_FoilHole_24015405_Data_24016401_24016403_20200225_0014_Fractions.mrc_rigid_aligned_particles.mrc' 13067.9072265625 12987.794921875 -1.0359375476837158 -2.151562452316284 -1.5252151489257812 2.8225245475769043 -0.14375591278076172 
19 b'J19/extract/16300063732981320253_FoilHole_24015405_Data_24016401_24016403_20200225_0014_Fractions.mrc_rigid_aligned_particles.mrc' 13067.9072265625 12987.794921875 -4.064062595367432 -0.23906250298023224 2.226463556289673 1.8688268661499023 0.49438008666038513 
20 b'J19/extract/16300063732981320253_FoilHole_24015405_Data_24016401_24016403_20200225_0014_Fractions.mrc_rigid_aligned_particles.mrc' 13067.9072265625 12987.794921875 -2.4703125953674316 -0.23906250298023224 -2.885637044906616 0.08765604346990585 2.4158005714416504 

========================================================================================== 
** Example of a text file whose name may appear in Parameter 14 when generating a 3D 
electrostatic potential on a set of 256 equispaced transverse (xy) planes uniformly distributed 
between z=-10 A and z=0 A:

//Z rotation	Y' rotation	Z'' rotation1	Defocus1
//defocus = -10 + 0.0390625 * (ROW() - 3)			
0	0	0	-10
0	0	0	-9.9609375
0	0	0	-9.921875
0	0	0	-9.8828125
0	0	0	-9.84375
0	0	0	-9.8046875
0	0	0	-9.765625
0	0	0	-9.7265625
0	0	0	-9.6875
0	0	0	-9.6484375
0	0	0	-9.609375
0	0	0	-9.5703125
0	0	0	-9.53125
0	0	0	-9.4921875
0	0	0	-9.453125
0	0	0	-9.4140625
0	0	0	-9.375
0	0	0	-9.3359375
0	0	0	-9.296875
0	0	0	-9.2578125
0	0	0	-9.21875
0	0	0	-9.1796875
0	0	0	-9.140625
0	0	0	-9.1015625
0	0	0	-9.0625
0	0	0	-9.0234375
0	0	0	-8.984375
0	0	0	-8.9453125
0	0	0	-8.90625
0	0	0	-8.8671875
0	0	0	-8.828125
0	0	0	-8.7890625
0	0	0	-8.75
0	0	0	-8.7109375
0	0	0	-8.671875
0	0	0	-8.6328125
0	0	0	-8.59375
0	0	0	-8.5546875
0	0	0	-8.515625
0	0	0	-8.4765625
0	0	0	-8.4375
0	0	0	-8.3984375
0	0	0	-8.359375
0	0	0	-8.3203125
0	0	0	-8.28125
0	0	0	-8.2421875
0	0	0	-8.203125
0	0	0	-8.1640625
0	0	0	-8.125
0	0	0	-8.0859375
0	0	0	-8.046875
0	0	0	-8.0078125
0	0	0	-7.96875
0	0	0	-7.9296875
0	0	0	-7.890625
0	0	0	-7.8515625
0	0	0	-7.8125
0	0	0	-7.7734375
0	0	0	-7.734375
0	0	0	-7.6953125
0	0	0	-7.65625
0	0	0	-7.6171875
0	0	0	-7.578125
0	0	0	-7.5390625
0	0	0	-7.5
0	0	0	-7.4609375
0	0	0	-7.421875
0	0	0	-7.3828125
0	0	0	-7.34375
0	0	0	-7.3046875
0	0	0	-7.265625
0	0	0	-7.2265625
0	0	0	-7.1875
0	0	0	-7.1484375
0	0	0	-7.109375
0	0	0	-7.0703125
0	0	0	-7.03125
0	0	0	-6.9921875
0	0	0	-6.953125
0	0	0	-6.9140625
0	0	0	-6.875
0	0	0	-6.8359375
0	0	0	-6.796875
0	0	0	-6.7578125
0	0	0	-6.71875
0	0	0	-6.6796875
0	0	0	-6.640625
0	0	0	-6.6015625
0	0	0	-6.5625
0	0	0	-6.5234375
0	0	0	-6.484375
0	0	0	-6.4453125
0	0	0	-6.40625
0	0	0	-6.3671875
0	0	0	-6.328125
0	0	0	-6.2890625
0	0	0	-6.25
0	0	0	-6.2109375
0	0	0	-6.171875
0	0	0	-6.1328125
0	0	0	-6.09375
0	0	0	-6.0546875
0	0	0	-6.015625
0	0	0	-5.9765625
0	0	0	-5.9375
0	0	0	-5.8984375
0	0	0	-5.859375
0	0	0	-5.8203125
0	0	0	-5.78125
0	0	0	-5.7421875
0	0	0	-5.703125
0	0	0	-5.6640625
0	0	0	-5.625
0	0	0	-5.5859375
0	0	0	-5.546875
0	0	0	-5.5078125
0	0	0	-5.46875
0	0	0	-5.4296875
0	0	0	-5.390625
0	0	0	-5.3515625
0	0	0	-5.3125
0	0	0	-5.2734375
0	0	0	-5.234375
0	0	0	-5.1953125
0	0	0	-5.15625
0	0	0	-5.1171875
0	0	0	-5.078125
0	0	0	-5.0390625
0	0	0	-5
0	0	0	-4.9609375
0	0	0	-4.921875
0	0	0	-4.8828125
0	0	0	-4.84375
0	0	0	-4.8046875
0	0	0	-4.765625
0	0	0	-4.7265625
0	0	0	-4.6875
0	0	0	-4.6484375
0	0	0	-4.609375
0	0	0	-4.5703125
0	0	0	-4.53125
0	0	0	-4.4921875
0	0	0	-4.453125
0	0	0	-4.4140625
0	0	0	-4.375
0	0	0	-4.3359375
0	0	0	-4.296875
0	0	0	-4.2578125
0	0	0	-4.21875
0	0	0	-4.1796875
0	0	0	-4.140625
0	0	0	-4.1015625
0	0	0	-4.0625
0	0	0	-4.0234375
0	0	0	-3.984375
0	0	0	-3.9453125
0	0	0	-3.90625
0	0	0	-3.8671875
0	0	0	-3.828125
0	0	0	-3.7890625
0	0	0	-3.75
0	0	0	-3.7109375
0	0	0	-3.671875
0	0	0	-3.6328125
0	0	0	-3.59375
0	0	0	-3.5546875
0	0	0	-3.515625
0	0	0	-3.4765625
0	0	0	-3.4375
0	0	0	-3.3984375
0	0	0	-3.359375
0	0	0	-3.3203125
0	0	0	-3.28125
0	0	0	-3.2421875
0	0	0	-3.203125
0	0	0	-3.1640625
0	0	0	-3.125
0	0	0	-3.0859375
0	0	0	-3.046875
0	0	0	-3.0078125
0	0	0	-2.96875
0	0	0	-2.9296875
0	0	0	-2.890625
0	0	0	-2.8515625
0	0	0	-2.8125
0	0	0	-2.7734375
0	0	0	-2.734375
0	0	0	-2.6953125
0	0	0	-2.65625
0	0	0	-2.6171875
0	0	0	-2.578125
0	0	0	-2.5390625
0	0	0	-2.5
0	0	0	-2.4609375
0	0	0	-2.421875
0	0	0	-2.3828125
0	0	0	-2.34375
0	0	0	-2.3046875
0	0	0	-2.265625
0	0	0	-2.2265625
0	0	0	-2.1875
0	0	0	-2.1484375
0	0	0	-2.109375
0	0	0	-2.0703125
0	0	0	-2.03125
0	0	0	-1.9921875
0	0	0	-1.953125
0	0	0	-1.9140625
0	0	0	-1.875
0	0	0	-1.8359375
0	0	0	-1.796875
0	0	0	-1.7578125
0	0	0	-1.71875
0	0	0	-1.6796875
0	0	0	-1.640625
0	0	0	-1.6015625
0	0	0	-1.5625
0	0	0	-1.5234375
0	0	0	-1.484375
0	0	0	-1.4453125
0	0	0	-1.40625
0	0	0	-1.3671875
0	0	0	-1.328125
0	0	0	-1.2890625
0	0	0	-1.25
0	0	0	-1.2109375
0	0	0	-1.171875
0	0	0	-1.1328125
0	0	0	-1.09375
0	0	0	-1.0546875
0	0	0	-1.015625
0	0	0	-0.9765625
0	0	0	-0.9375
0	0	0	-0.8984375
0	0	0	-0.859375
0	0	0	-0.8203125
0	0	0	-0.78125
0	0	0	-0.7421875
0	0	0	-0.703125
0	0	0	-0.6640625
0	0	0	-0.625
0	0	0	-0.5859375
0	0	0	-0.546875
0	0	0	-0.5078125
0	0	0	-0.46875
0	0	0	-0.4296875
0	0	0	-0.390625
0	0	0	-0.3515625
0	0	0	-0.3125
0	0	0	-0.2734375
0	0	0	-0.234375
0	0	0	-0.1953125
0	0	0	-0.15625
0	0	0	-0.1171875
0	0	0	-0.078125
0	0	0	-0.0390625

===========================================================================
Specifications of the GRD file format:

GRD files contain:
(1) string "DSAA" or "DSBB" (in ASC and BIN GRD files, respectively)
(2) nx ny - array dimensions as "short"s (in ASC format, separated by a white space)
(3) xlo xhi - lower and higher X boundaries as "double"s (in ASC format, separated by a 
white space)
(4) ylo yhi - lower and higher Y boundaries as "double"s (in ASC format, separated by a 
white space)
(5) zlo zhi - minimum and maximum array value as "double"s (in ASC format, separated by a 
white space)
(6) u[i][j] - array values as "float"s (in ASC format, separated by a white space), j index 
changes most rapidly
 
===========================================================================
Specifications of the GRC file format

GRC files contain:
(1) -5 - as a "short" value (this is the GRC file format identifier)
(2) wl - wavelength as a "double"
(3) nx ny - array dimensions as "short"s (in ASC format, separated by a white space)
(4) xlo xhi - lower and higher X boundaries as "double"s (in ASC format, separated by a 
white space)
(5) ylo yhi - lower and higher Y boundaries as "double"s (in ASC format, separated by a 
white space)
(6) real(u[i][j]) imag(u[i][j]) - real and imaginary parts of array values as "float"s (in ASC 
format, separated by a white space), j index changes most rapidly.