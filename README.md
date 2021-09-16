# IncBicrystal
Incompatibility stresses in an infinite bi-crystal with a planar interface

This MATLAB code computes the stress fields in an infinite bi-crystal with a planar interface. 	

The bi-crystal is submitted to a macroscopic homogeneous and remotely applied stress as well as to piecewise uniform plastic strains. 

The bi-crystal can have different grain volume fractions.

The two grains are assumed perfectly bonded  along the planar interface whose normal is along e2. 

The code provides also the effective stiffness tensor, the macroscopic strain tensor and the effective plastic strain tensor.

References for the model can be found in:

T. Richeton, S. Berbenni, Eur. J. Mech. A. Solids 37 (2013) 231–247

T. Richeton, S. Berbenni, Int. J. Sol. Struct. 51 (2014) 794-807

T. Richeton, I. Tiba, S. Berbenni, O. Bouaziz, Phil. Mag. 95 (2015) 12-31

I. Tiba, T. Richeton, C. Motz, H. Vehoff, S. Berbenni, Acta Mater. 83 (2015) 227-238

T. Richeton, Crystals 7 (2017) 203 (1-14)

T. Richeton, Scripta Mater. 169 (2019) 14-18

The contracted Voigt notation (11:1, 22: 2, 33: 3, 23: 4, 31:5, 12: 6) is used for stresses and strains in vector notation and stiffnesses and compliances in 6x6 matrix notation.
