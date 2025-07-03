'#Major' N (number of PSLs)
Write_PSL_info(1st_PSL)
Write_PSL_info(2nd_PSL)
:
:
:
Write_PSL_info(Nth__PSL)
'#Medium' N (number of PSLs)
Write_PSL_info(1st_PSL)
Write_PSL_info(2nd_PSL)
:
:
:
Write_PSL_info(Nth__PSL)
'#Minor' N (number of PSLs)
Write_PSL_info(1st_PSL)
Write_PSL_info(2nd_PSL)
:
:
:
Write_PSL_info(Nth__PSL)
#Outline "Cartesian" or "Unstructured"
#Vertices numVertices
xx_1 yy_1 zz_1
xx_2 yy_2 zz_2
:
:
:
xx_numVertices yy_numVertices zz_numVertices
#Faces numFaces //vertex indices range from 0 ~ numVertices-1
f_1_1 f_1_2 f_1_3
f_2_1 f_2_2 f_2_3
:
:
:
f_numFaces_1 f_numFaces_2 f_numFaces_3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Comment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Comment 1
function Write_PSL_info(ith_PSL)
	 L	IM1 (geo-based) IM2 (PS-based) IM3 (vM-based) IM4 (Length-based) AppearanceOrder Seed(x, y, z)
	'spatial coordinates:' X_1 Y_1 Z_1	X_2 Y_2 Z_2 ..... X_L Y_L Z_L
	'band vertex coordinates (unSmoothed) :' x_1 y_1 z_1	x_2 y_2 z_2 ..... x_2L y_2L z_2L
	'band vertex coordinates (smoothed) :' x_1 y_1 z_1	x_2 y_2 z_2 ..... x_2L y_2L z_2L
	'Scalar Field:' s_1_1 s_1_2 ..... s_1_L //principal stress (Sigma in my test code)
					s_2_1 s_2_2 ..... s_2_L	//von Mises stress (Sigma_vM)
					s_3_1 s_3_2 ..... s_3_L	//normal stress on x-x dir (Sigma_xx)
					s_4_1 s_4_2 ..... s_4_L //normal stress on y-y dir (Sigma_yy)
					s_5_1 s_5_2 ..... s_5_L //normal stress on z-z dir (Sigma_zz)
					s_6_1 s_6_2 ..... s_6_L //shear stress on y-z dir (Sigma_yz)
					s_7_1 s_7_2 ..... s_7_L //shear stress on z-x dir (Sigma_zx)
					s_8_1 s_8_2 ..... s_8_L //shear stress on x-y dir (Sigma_xy)
end
%%

%% Comment 2
L: 'length of ith major PSL'
IM: 'Importance Metirc'
%%

%% Comment 3
	band geometry	
	p_1 ----- p_3 ----- p_5 ----- p_7 ..... p_(2L-1) 
	|
	P_1 ----- P_2 ----- P_3 ----- P_4 ..... P_L
	|
	p_2 ----- p_4 ----- p_6 ----- p_8 ..... p_2L
	
	P_i = (X_i, Y_i, Z_i)	
	p_i = (x_i, y_i, z_i)	
	
	In my practice, I draw the band via patches like (p_1, p_2, p_4, p_3), (p_3, p_4, p_6, p_5) .....
	P_(2i-1), P_(2i) and p_i share same scalar value s_i 
	
	In case you want to adjust the width of band:
	p_1_new = P_1 + scalingFactor * (p_1-P_1);
	p_2_new = P_1 + scalingFactor * (p_2-P_1);
%%

%% Comment 3
	"Cartesian" == needs to be smoothed
	"Unstructured" == no need to be smoothed
%%



