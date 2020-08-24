#include "colors.inc"
#include "finish.inc"

global_settings {assumed_gamma 1 max_trace_level 6}
background {color White transmit 1.0}
camera {orthographic
  right -26.76*x up 27.65*y
  direction 1.00*z
  location <0,0,50.00> look_at <0,0,0>}
light_source {<  2.00,   3.00,  40.00> color White
  area_light <0.70, 0, 0>, <0, 0.70, 0>, 3, 3
  adaptive 1 jitter}

#declare simple = finish {phong 0.7}
#declare pale = finish {ambient 0.5 diffuse 0.85 roughness 0.001 specular 0.200 }
#declare intermediate = finish {ambient 0.3 diffuse 0.6 specular 0.1 roughness 0.04}
#declare vmd = finish {ambient 0.0 diffuse 0.65 phong 0.1 phong_size 40.0 specular 0.5 }
#declare jmol = finish {ambient 0.2 diffuse 0.6 specular 1 roughness 0.001 metallic}
#declare ase2 = finish {ambient 0.05 brilliance 3 diffuse 0.6 metallic specular 0.7 roughness 0.04 reflection 0.15}
#declare ase3 = finish {ambient 0.15 brilliance 2 diffuse 0.6 metallic specular 1.0 roughness 0.001 reflection 0.0}
#declare glass = finish {ambient 0.05 diffuse 0.3 specular 1.0 roughness 0.001}
#declare glass2 = finish {ambient 0.01 diffuse 0.3 specular 1.0 reflection 0.25 roughness 0.001}
#declare Rcell = 0.070;
#declare Rbond = 0.100;

#macro atom(LOC, R, COL, TRANS, FIN)
  sphere{LOC, R texture{pigment{color COL transmit TRANS} finish{FIN}}}
#end
#macro constrain(LOC, R, COL, TRANS FIN)
union{torus{R, Rcell rotate 45*z texture{pigment{color COL transmit TRANS} finish{FIN}}}
      torus{R, Rcell rotate -45*z texture{pigment{color COL transmit TRANS} finish{FIN}}}
      translate LOC}
#end

atom(< -0.34,  -1.14,  -6.77>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #0 
atom(< -1.52,  -1.66,  -6.16>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #1 
atom(<  0.90,  -1.65,  -6.43>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #2 
atom(< -1.43,  -2.68,  -5.23>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #3 
atom(< -2.24,  -2.64,  -4.04>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #4 
atom(< -2.43,  -0.57,  -5.94>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #5 
atom(< -3.21,  -0.53,  -4.80>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #6 
atom(< -3.11,  -1.59,  -3.83>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #7 
atom(<  1.86,   0.60,  -6.39>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #8 
atom(<  0.57,   1.13,  -6.74>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #9 
atom(<  2.02,  -0.77,  -6.23>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #10 
atom(< -0.51,   0.28,  -6.93>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #11 
atom(< -1.80,   0.63,  -6.42>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #12 
atom(<  0.39,   2.36,  -6.03>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #13 
atom(< -0.85,   2.70,  -5.53>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #14 
atom(< -1.97,   1.82,  -5.73>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #15 
atom(<  3.40,  -0.41,  -4.25>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #16 
atom(<  3.23,   1.01,  -4.41>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #17 
atom(<  2.81,  -1.28,  -5.15>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #18 
atom(<  2.47,   1.50,  -5.46>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #19 
atom(<  1.57,   2.60,  -5.24>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #20 
atom(<  3.11,   1.59,  -3.11>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #21 
atom(<  2.24,   2.64,  -2.89>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #22 
atom(<  1.45,   3.16,  -3.98>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #23 
atom(<  2.16,  -2.77,  -3.32>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #24 
atom(<  2.78,  -1.86,  -2.39>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #25 
atom(<  2.18,  -2.48,  -4.68>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #26 
atom(<  3.39,  -0.70,  -2.85>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #27 
atom(<  3.21,   0.53,  -2.14>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #28 
atom(<  1.97,  -1.82,  -1.20>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #29 
atom(<  1.80,  -0.63,  -0.52>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #30 
atom(<  2.43,   0.57,  -1.00>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #31 
atom(< -0.14,  -3.21,  -4.88>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #32 
atom(< -0.16,  -3.51,  -3.47>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #33 
atom(<  1.00,  -2.71,  -5.47>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #34 
atom(<  0.97,  -3.29,  -2.71>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #35 
atom(<  0.85,  -2.70,  -1.40>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #36 
atom(< -1.45,  -3.16,  -2.95>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #37 
atom(< -1.57,  -2.60,  -1.69>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #38 
atom(< -0.39,  -2.36,  -0.90>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #39 
atom(<  0.34,   1.14,  -0.16>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #40 
atom(<  0.51,  -0.28,   0.00>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #41 
atom(< -0.57,  -1.13,  -0.19>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #42 
atom(< -0.90,   1.65,  -0.50>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #43 
atom(<  0.14,   3.21,  -2.06>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #44 
atom(<  1.43,   2.68,  -1.71>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #45 
atom(<  1.52,   1.66,  -0.77>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #46 
atom(< -1.00,   2.71,  -1.47>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #47 
atom(< -2.16,   2.77,  -3.62>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #48 
atom(< -0.97,   3.29,  -4.23>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #49 
atom(<  0.16,   3.51,  -3.47>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #50 
atom(< -2.18,   2.48,  -2.26>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #51 
atom(< -3.40,   0.41,  -2.68>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #52 
atom(< -3.39,   0.70,  -4.08>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #53 
atom(< -2.78,   1.86,  -4.54>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #54 
atom(< -2.81,   1.28,  -1.79>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #55 
atom(< -1.86,  -0.60,  -0.54>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #56 
atom(< -2.47,  -1.50,  -1.47>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #57 
atom(< -3.23,  -1.01,  -2.52>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #58 
atom(< -2.02,   0.77,  -0.70>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #59 
