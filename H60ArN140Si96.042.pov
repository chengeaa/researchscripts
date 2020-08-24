#include "colors.inc"
#include "finish.inc"

global_settings {assumed_gamma 1 max_trace_level 6}
background {color White transmit 1.0}
camera {orthographic
  right -63.09*x up 44.84*y
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

cylinder {< -1.55, -16.23,  -5.22>, < 13.88, -15.20,  -7.55>, Rcell pigment {Black}}
cylinder {<-11.41, -14.43, -17.22>, <  4.02, -13.40, -19.56>, Rcell pigment {Black}}
cylinder {<-12.52,  14.32, -11.99>, <  2.90,  15.34, -14.33>, Rcell pigment {Black}}
cylinder {< -2.66,  12.52,   0.01>, < 12.77,  13.54,  -2.33>, Rcell pigment {Black}}
cylinder {< -1.55, -16.23,  -5.22>, <-11.41, -14.43, -17.22>, Rcell pigment {Black}}
cylinder {< 13.88, -15.20,  -7.55>, <  4.02, -13.40, -19.56>, Rcell pigment {Black}}
cylinder {< 12.77,  13.54,  -2.33>, <  2.90,  15.34, -14.33>, Rcell pigment {Black}}
cylinder {< -2.66,  12.52,   0.01>, <-12.52,  14.32, -11.99>, Rcell pigment {Black}}
cylinder {< -1.55, -16.23,  -5.22>, < -2.66,  12.52,   0.01>, Rcell pigment {Black}}
cylinder {< 13.88, -15.20,  -7.55>, < 12.77,  13.54,  -2.33>, Rcell pigment {Black}}
cylinder {<  4.02, -13.40, -19.56>, <  2.90,  15.34, -14.33>, Rcell pigment {Black}}
cylinder {<-11.41, -14.43, -17.22>, <-12.52,  14.32, -11.99>, Rcell pigment {Black}}
atom(< -4.09, -13.73,  -9.75>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #0 
atom(< -1.12, -13.36, -10.52>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #1 
atom(< -1.89, -13.93,  -7.59>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #2 
atom(<  1.34, -12.10,  -9.15>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #3 
atom(<  0.56, -12.68,  -6.19>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #4 
atom(<  3.55, -12.33,  -6.94>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #5 
atom(< -4.20, -10.75,  -9.16>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #6 
atom(< -1.94, -11.03,  -7.08>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #7 
atom(< -1.24, -10.55,  -9.99>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #8 
atom(<  1.02,  -9.08,  -8.68>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #9 
atom(<  0.38,  -9.84,  -5.66>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #10 
atom(<  3.31,  -9.46,  -6.44>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #11 
atom(< -4.52,  -7.85,  -8.71>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #12 
atom(< -2.06,  -8.07,  -6.48>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #13 
atom(< -1.44,  -7.68,  -9.47>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #14 
atom(<  0.96,  -6.33,  -8.07>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #15 
atom(<  0.45,  -6.92,  -5.16>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #16 
atom(<  3.32,  -6.54,  -5.83>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #17 
atom(< -4.20,  -5.12,  -8.10>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #18 
atom(< -2.16,  -5.25,  -5.89>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #19 
atom(< -1.50,  -4.69,  -8.85>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #20 
atom(<  0.82,  -3.48,  -7.47>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #21 
atom(<  0.08,  -3.89,  -4.40>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #22 
atom(<  3.11,  -3.53,  -5.51>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #23 
atom(<  0.74, -14.26,  -5.50>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #24 
atom(<  0.47, -13.50,  -9.73>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #25 
atom(< -2.36, -13.65,  -9.27>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #26 
atom(<  4.34, -13.70,  -7.81>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #27 
atom(<  1.84, -12.16,  -7.41>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #28 
atom(< -4.96, -12.36,  -8.86>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #29 
atom(< -1.00, -11.76, -11.29>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #30 
atom(< -1.00, -12.51,  -7.02>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #31 
atom(<  0.89, -11.39,  -4.97>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #32 
atom(<  0.36, -10.62,  -9.27>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #33 
atom(< -2.47, -10.81,  -8.77>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #34 
atom(<  4.36, -10.69,  -7.13>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #35 
atom(<  1.61,  -9.61,  -6.97>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #36 
atom(< -5.09,  -9.43,  -8.32>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #37 
atom(< -1.45,  -8.93, -10.65>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #38 
atom(< -1.14,  -9.57,  -6.51>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #39 
atom(<  0.49,  -8.54,  -4.34>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #40 
atom(< -0.08,  -7.74,  -8.39>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #41 
atom(< -2.79,  -8.07,  -8.23>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #42 
atom(<  3.88,  -7.85,  -6.83>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #43 
atom(<  1.60,  -6.38,  -6.39>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #44 
atom(< -5.02,  -6.59,  -7.63>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #45 
atom(< -1.61,  -5.99, -10.11>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #46 
atom(< -1.15,  -6.70,  -5.88>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #47 
atom(< -2.61,  -5.01,  -7.53>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #48 
atom(<  0.37,  -5.62,  -3.88>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #49 
atom(< -0.02,  -4.90,  -7.98>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #50 
atom(<  3.97,  -4.97,  -6.18>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #51 
atom(<  1.34,  -3.58,  -5.76>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #52 
atom(< -5.19,  -3.88,  -7.46>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #53 
atom(< -1.84,  -3.06,  -9.54>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #54 
atom(< -1.52,  -3.67,  -5.18>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #55 
atom(< -0.03,  -1.94,  -7.75>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #56 
atom(<  0.64,  -2.74,  -3.03>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #57 
atom(<  3.98,  -2.13,  -6.16>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #58 
atom(< -4.80, -14.96,  -9.33>, 0.28, rgb <1.00, 1.00, 1.00>, 0.0, ase2) // #59 
atom(< -1.24, -14.46, -11.59>, 0.28, rgb <1.00, 1.00, 1.00>, 0.0, ase2) // #60 
atom(< -1.02, -15.17,  -7.44>, 0.28, rgb <1.00, 1.00, 1.00>, 0.0, ase2) // #61 
atom(<  1.65, -14.65,  -5.88>, 0.28, rgb <1.00, 1.00, 1.00>, 0.0, ase2) // #62 
atom(<  0.47, -14.21,  -8.98>, 0.28, rgb <1.00, 1.00, 1.00>, 0.0, ase2) // #63 
atom(<  3.59, -14.09,  -8.47>, 0.28, rgb <1.00, 1.00, 1.00>, 0.0, ase2) // #64 
atom(< -5.07,  -3.52,  -6.50>, 0.28, rgb <1.00, 1.00, 1.00>, 0.0, ase2) // #65 
atom(< -2.46,  -2.37,  -9.02>, 0.28, rgb <1.00, 1.00, 1.00>, 0.0, ase2) // #66 
atom(< -1.39,  -2.91,  -5.90>, 0.28, rgb <1.00, 1.00, 1.00>, 0.0, ase2) // #67 
atom(<  0.71,  -1.18,  -7.76>, 0.28, rgb <1.00, 1.00, 1.00>, 0.0, ase2) // #68 
atom(< -4.98,  -1.76,  -8.39>, 0.28, rgb <1.00, 1.00, 1.00>, 0.0, ase2) // #69 
atom(< -4.30,  -1.54,  -3.21>, 0.28, rgb <1.00, 1.00, 1.00>, 0.0, ase2) // #70 
atom(< -0.51,  -1.64,  -6.88>, 0.28, rgb <1.00, 1.00, 1.00>, 0.0, ase2) // #71 
atom(<  1.57,  -3.19,  -2.84>, 0.28, rgb <1.00, 1.00, 1.00>, 0.0, ase2) // #72 
atom(<  3.22,  -1.52,  -6.54>, 0.28, rgb <1.00, 1.00, 1.00>, 0.0, ase2) // #73 
atom(< -9.03, -12.71, -15.74>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #74 
atom(< -6.05, -12.45, -16.52>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #75 
atom(< -6.84, -13.03, -13.60>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #76 
atom(< -3.64, -11.08, -15.17>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #77 
atom(< -4.36, -11.76, -12.11>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #78 
atom(< -1.38, -11.44, -13.02>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #79 
atom(< -9.10,  -9.82, -15.24>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #80 
atom(< -6.95, -10.22, -13.00>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #81 
atom(< -6.02,  -9.57, -15.92>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #82 
atom(< -3.69,  -8.06, -14.49>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #83 
atom(< -4.43,  -8.84, -11.70>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #84 
atom(< -1.49,  -8.52, -12.41>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #85 
atom(< -9.24,  -7.03, -14.53>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #86 
atom(< -7.03,  -7.27, -12.47>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #87 
atom(< -6.31,  -6.62, -15.41>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #88 
atom(< -3.53,  -5.18, -13.97>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #89 
atom(< -4.52,  -5.84, -11.14>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #90 
atom(< -1.56,  -5.52, -11.84>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #91 
atom(< -9.21,  -4.26, -14.12>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #92 
atom(< -7.03,  -4.37, -12.00>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #93 
atom(< -6.40,  -3.70, -14.93>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #94 
atom(< -4.83,  -2.48, -12.97>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #95 
atom(< -4.79,  -3.09, -10.32>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #96 
atom(< -1.58,  -2.66, -11.23>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #97 
atom(< -4.18, -13.42, -11.51>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #98 
atom(< -4.46, -12.60, -15.74>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #99 
atom(< -7.30, -12.74, -15.29>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #100 
atom(< -0.59, -12.80, -13.82>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #101 
atom(< -3.15, -11.51, -13.42>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #102 
atom(< -9.89, -11.39, -14.94>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #103 
atom(< -6.30, -10.78, -17.15>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #104 
atom(< -6.06, -11.69, -12.70>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #105 
atom(< -4.28, -10.40, -10.91>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #106 
atom(< -4.56,  -9.59, -14.94>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #107 
atom(< -7.37,  -9.99, -14.74>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #108 
atom(< -0.60,  -9.86, -13.18>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #109 
atom(< -3.15,  -8.47, -12.82>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #110 
atom(<-10.05,  -8.59, -14.41>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #111 
atom(< -6.64,  -8.02, -16.60>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #112 
atom(< -6.00,  -8.77, -12.52>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #113 
atom(< -4.24,  -7.50, -10.49>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #114 
atom(< -4.57,  -6.63, -14.55>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #115 
atom(< -7.51,  -7.10, -14.18>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #116 
atom(< -0.71,  -6.88, -12.65>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #117 
atom(< -3.24,  -5.62, -12.34>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #118 
atom(< -9.98,  -5.77, -13.60>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #119 
atom(< -6.66,  -5.02, -16.09>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #120 
atom(< -6.04,  -5.80, -12.14>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #121 
atom(< -7.57,  -4.02, -13.68>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #122 
atom(< -4.35,  -4.70,  -9.84>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #123 
atom(< -4.70,  -3.69, -14.23>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #124 
atom(< -0.61,  -3.90, -12.04>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #125 
atom(< -3.33,  -2.76, -11.89>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #126 
atom(<-10.14,  -2.91, -13.40>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #127 
atom(< -6.54,  -2.04, -15.71>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #128 
atom(< -6.03,  -2.72, -11.66>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #129 
atom(< -4.80,  -0.93, -13.58>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #130 
atom(< -4.45,  -1.60,  -9.28>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #131 
atom(< -0.78,  -1.05, -11.50>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #132 
atom(< -9.70, -13.98, -15.44>, 0.28, rgb <1.00, 1.00, 1.00>, 0.0, ase2) // #133 
atom(< -6.10, -13.55, -17.58>, 0.28, rgb <1.00, 1.00, 1.00>, 0.0, ase2) // #134 
atom(< -5.92, -14.26, -13.50>, 0.28, rgb <1.00, 1.00, 1.00>, 0.0, ase2) // #135 
atom(< -3.31, -13.80, -11.93>, 0.28, rgb <1.00, 1.00, 1.00>, 0.0, ase2) // #136 
atom(< -4.60, -13.23, -14.91>, 0.28, rgb <1.00, 1.00, 1.00>, 0.0, ase2) // #137 
atom(< -1.27, -13.33, -14.39>, 0.28, rgb <1.00, 1.00, 1.00>, 0.0, ase2) // #138 
atom(< -9.56,  -2.13, -13.04>, 0.28, rgb <1.00, 1.00, 1.00>, 0.0, ase2) // #139 
atom(< -6.06,  -1.22, -15.17>, 0.28, rgb <1.00, 1.00, 1.00>, 0.0, ase2) // #140 
atom(< -6.69,  -1.97, -11.39>, 0.28, rgb <1.00, 1.00, 1.00>, 0.0, ase2) // #141 
atom(< -3.81,  -0.72, -13.81>, 0.28, rgb <1.00, 1.00, 1.00>, 0.0, ase2) // #142 
atom(< -8.98,  -0.04, -15.41>, 0.28, rgb <1.00, 1.00, 1.00>, 0.0, ase2) // #143 
atom(< -8.59,  -1.01, -10.24>, 0.28, rgb <1.00, 1.00, 1.00>, 0.0, ase2) // #144 
atom(< -5.10,  -0.14, -12.97>, 0.28, rgb <1.00, 1.00, 1.00>, 0.0, ase2) // #145 
atom(< -4.70,  -0.71,  -9.72>, 0.28, rgb <1.00, 1.00, 1.00>, 0.0, ase2) // #146 
atom(< -0.55,  -1.05, -12.53>, 0.28, rgb <1.00, 1.00, 1.00>, 0.0, ase2) // #147 
atom(<  3.63, -13.10, -10.90>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #148 
atom(<  6.59, -12.86, -11.69>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #149 
atom(<  5.82, -13.40,  -8.75>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #150 
atom(<  8.98, -11.58, -10.24>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #151 
atom(<  8.32, -12.08,  -7.30>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #152 
atom(< 11.20, -11.92,  -8.11>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #153 
atom(<  3.45, -10.19, -10.37>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #154 
atom(<  5.74, -10.47,  -8.24>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #155 
atom(<  6.42,  -9.93, -11.15>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #156 
atom(<  8.74,  -8.60,  -9.98>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #157 
atom(<  8.18,  -9.14,  -6.81>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #158 
atom(< 11.29,  -8.86,  -7.74>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #159 
atom(<  3.32,  -7.31,  -9.98>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #160 
atom(<  5.38,  -7.74,  -7.73>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #161 
atom(<  6.25,  -7.01, -10.61>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #162 
atom(<  8.62,  -5.63,  -9.87>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #163 
atom(<  7.91,  -6.50,  -6.43>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #164 
atom(< 10.97,  -5.95,  -7.12>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #165 
atom(<  3.17,  -4.31,  -9.43>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #166 
atom(<  5.42,  -4.77,  -7.19>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #167 
atom(<  6.23,  -4.20, -10.08>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #168 
atom(<  8.68,  -2.93, -10.31>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #169 
atom(<  8.23,  -3.84,  -7.50>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #170 
atom(< 10.09,  -2.57,  -6.20>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #171 
atom(<  8.46, -13.75,  -6.67>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #172 
atom(<  8.18, -13.02, -10.90>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #173 
atom(<  5.35, -13.14, -10.45>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #174 
atom(< 12.05, -13.19,  -8.99>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #175 
atom(<  9.52, -11.65,  -8.53>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #176 
atom(<  2.68, -11.78, -10.22>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #177 
atom(<  6.66, -11.23, -12.38>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #178 
atom(<  6.85, -11.93,  -8.39>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #179 
atom(<  8.35, -10.73,  -6.09>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #180 
atom(<  7.92, -10.13, -10.20>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #181 
atom(<  5.17, -10.21,  -9.93>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #182 
atom(< 11.93, -10.24,  -8.59>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #183 
atom(<  9.41,  -9.06,  -8.17>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #184 
atom(<  2.57,  -8.85,  -9.58>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #185 
atom(<  6.13,  -8.33, -11.84>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #186 
atom(<  6.63,  -9.04,  -7.71>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #187 
atom(<  8.43,  -7.86,  -5.59>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #188 
atom(<  7.83,  -7.15,  -9.89>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #189 
atom(<  5.00,  -7.40,  -9.48>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #190 
atom(< 11.81,  -7.19,  -7.95>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #191 
atom(<  9.13,  -6.23,  -7.71>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #192 
atom(<  2.36,  -6.02,  -9.21>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #193 
atom(<  6.31,  -5.40, -11.35>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #194 
atom(<  6.29,  -6.27,  -7.44>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #195 
atom(<  4.88,  -4.42,  -8.84>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #196 
atom(<  7.66,  -4.97,  -5.17>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #197 
atom(<  8.56,  -4.24, -11.68>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #198 
atom(< 11.59,  -4.20,  -7.37>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #199 
atom(<  7.98,  -4.12,  -9.33>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #200 
atom(<  2.29,  -3.20,  -8.39>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #201 
atom(<  5.53,  -2.73, -10.86>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #202 
atom(<  6.49,  -3.34,  -6.87>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #203 
atom(<  7.80,  -1.52, -10.67>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #204 
atom(<  5.69,  -8.92,  -5.05>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #205 
atom(< 10.25,  -0.85,  -6.03>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #206 
atom(<  2.97, -14.39, -10.68>, 0.28, rgb <1.00, 1.00, 1.00>, 0.0, ase2) // #207 
atom(<  6.44, -13.82, -12.84>, 0.28, rgb <1.00, 1.00, 1.00>, 0.0, ase2) // #208 
atom(<  6.71, -14.53,  -8.68>, 0.28, rgb <1.00, 1.00, 1.00>, 0.0, ase2) // #209 
atom(<  9.34, -14.17,  -7.07>, 0.28, rgb <1.00, 1.00, 1.00>, 0.0, ase2) // #210 
atom(<  8.33, -13.89, -10.33>, 0.28, rgb <1.00, 1.00, 1.00>, 0.0, ase2) // #211 
atom(< 11.44, -13.94,  -9.43>, 0.28, rgb <1.00, 1.00, 1.00>, 0.0, ase2) // #212 
atom(<  2.76,  -2.33,  -8.09>, 0.28, rgb <1.00, 1.00, 1.00>, 0.0, ase2) // #213 
atom(<  4.54,  -2.70, -10.56>, 0.28, rgb <1.00, 1.00, 1.00>, 0.0, ase2) // #214 
atom(<  5.95,  -2.56,  -7.22>, 0.28, rgb <1.00, 1.00, 1.00>, 0.0, ase2) // #215 
atom(<  8.00,  -0.63, -11.14>, 0.28, rgb <1.00, 1.00, 1.00>, 0.0, ase2) // #216 
atom(<  3.75,  -0.96, -10.17>, 0.28, rgb <1.00, 1.00, 1.00>, 0.0, ase2) // #217 
atom(<  4.48,  -1.55,  -5.44>, 0.28, rgb <1.00, 1.00, 1.00>, 0.0, ase2) // #218 
atom(<  6.87,  -1.22, -10.28>, 0.28, rgb <1.00, 1.00, 1.00>, 0.0, ase2) // #219 
atom(<  5.62,  -4.89,  -2.88>, 0.28, rgb <1.00, 1.00, 1.00>, 0.0, ase2) // #220 
atom(< 10.18,  -0.46,  -6.99>, 0.28, rgb <1.00, 1.00, 1.00>, 0.0, ase2) // #221 
atom(< -1.31, -12.19, -16.90>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #222 
atom(<  1.66, -11.94, -17.69>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #223 
atom(<  0.89, -12.50, -14.76>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #224 
atom(<  4.03, -10.65, -16.33>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #225 
atom(<  3.32, -11.19, -13.37>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #226 
atom(<  6.25, -11.07, -14.13>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #227 
atom(< -1.40,  -9.19, -16.28>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #228 
atom(<  0.83,  -9.43, -14.11>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #229 
atom(<  1.50,  -8.93, -17.15>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #230 
atom(<  3.94,  -7.72, -15.75>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #231 
atom(<  3.21,  -8.15, -12.80>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #232 
atom(<  6.16,  -8.11, -13.57>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #233 
atom(< -1.50,  -6.32, -15.88>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #234 
atom(<  0.74,  -6.64, -13.62>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #235 
atom(<  1.33,  -6.10, -16.72>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #236 
atom(<  3.87,  -4.84, -15.23>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #237 
atom(<  3.18,  -5.22, -12.16>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #238 
atom(<  6.11,  -5.13, -13.13>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #239 
atom(< -1.42,  -3.52, -15.58>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #240 
atom(<  0.62,  -3.78, -13.26>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #241 
atom(<  1.48,  -3.14, -16.29>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #242 
atom(<  3.81,  -2.01, -14.77>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #243 
atom(<  2.81,  -2.36, -11.89>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #244 
atom(<  5.93,  -2.28, -12.59>, 0.99, rgb <0.94, 0.78, 0.62>, 0.0, ase2) // #245 
atom(<  3.52, -12.85, -12.66>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #246 
atom(<  3.26, -12.09, -16.91>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #247 
atom(<  0.42, -12.22, -16.45>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #248 
atom(<  7.11, -12.31, -15.00>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #249 
atom(<  4.54, -10.86, -14.62>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #250 
atom(< -1.96, -10.85, -15.83>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #251 
atom(<  1.36, -10.29, -18.35>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #252 
atom(<  1.71, -11.00, -14.15>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #253 
atom(<  3.32,  -9.83, -12.15>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #254 
atom(<  3.02,  -9.16, -16.28>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #255 
atom(<  0.32,  -9.09, -15.79>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #256 
atom(<  6.87,  -9.43, -14.57>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #257 
atom(<  4.47,  -8.15, -14.05>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #258 
atom(< -2.10,  -7.87, -15.27>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #259 
atom(<  1.35,  -7.34, -17.99>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #260 
atom(<  1.77,  -7.98, -13.79>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #261 
atom(<  3.35,  -6.85, -11.63>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #262 
atom(<  2.98,  -6.23, -15.96>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #263 
atom(<  0.20,  -6.35, -15.37>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #264 
atom(<  6.92,  -6.51, -13.85>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #265 
atom(<  4.38,  -5.25, -13.54>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #266 
atom(< -2.10,  -4.89, -14.82>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #267 
atom(<  1.25,  -4.49, -17.42>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #268 
atom(<  1.74,  -5.23, -13.22>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #269 
atom(<  0.31,  -3.73, -15.00>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #270 
atom(<  3.12,  -3.88, -11.10>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #271 
atom(<  2.93,  -3.39, -15.39>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #272 
atom(<  6.87,  -3.54, -13.44>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #273 
atom(<  4.17,  -2.18, -13.07>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #274 
atom(< -1.56,  -1.59, -14.38>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #275 
atom(<  1.05,  -1.42, -16.68>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #276 
atom(<  1.31,  -2.20, -12.82>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #277 
atom(<  3.13,  -0.36, -15.09>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #278 
atom(<  2.81,  -1.18, -10.53>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #279 
atom(<  6.32,  -0.65, -13.20>, 0.63, rgb <0.18, 0.31, 0.97>, 0.0, ase2) // #280 
atom(< -2.03, -13.39, -16.49>, 0.28, rgb <1.00, 1.00, 1.00>, 0.0, ase2) // #281 
atom(<  1.48, -12.94, -18.78>, 0.28, rgb <1.00, 1.00, 1.00>, 0.0, ase2) // #282 
atom(<  1.87, -13.64, -14.74>, 0.28, rgb <1.00, 1.00, 1.00>, 0.0, ase2) // #283 
atom(<  4.37, -13.27, -13.11>, 0.28, rgb <1.00, 1.00, 1.00>, 0.0, ase2) // #284 
atom(<  3.19, -12.86, -16.24>, 0.28, rgb <1.00, 1.00, 1.00>, 0.0, ase2) // #285 
atom(<  6.46, -12.95, -15.48>, 0.28, rgb <1.00, 1.00, 1.00>, 0.0, ase2) // #286 
atom(< -0.50,  -1.30, -14.08>, 0.28, rgb <1.00, 1.00, 1.00>, 0.0, ase2) // #287 
atom(<  0.79,  -0.86, -15.82>, 0.28, rgb <1.00, 1.00, 1.00>, 0.0, ase2) // #288 
atom(<  1.40,  -1.78, -13.79>, 0.28, rgb <1.00, 1.00, 1.00>, 0.0, ase2) // #289 
atom(<  2.13,  -0.32, -14.77>, 0.28, rgb <1.00, 1.00, 1.00>, 0.0, ase2) // #290 
atom(< -4.05,  -7.90, -17.73>, 0.28, rgb <1.00, 1.00, 1.00>, 0.0, ase2) // #291 
atom(< -1.63,  -0.39, -11.40>, 0.28, rgb <1.00, 1.00, 1.00>, 0.0, ase2) // #292 
atom(<  3.05,  -0.16, -16.13>, 0.28, rgb <1.00, 1.00, 1.00>, 0.0, ase2) // #293 
atom(<  2.37,  -1.73,  -9.77>, 0.28, rgb <1.00, 1.00, 1.00>, 0.0, ase2) // #294 
atom(<  5.42,  -0.16, -13.39>, 0.28, rgb <1.00, 1.00, 1.00>, 0.0, ase2) // #295 
atom(<  9.14,  -3.27,  -6.79>, 0.94, rgb <0.50, 0.82, 0.89>, 0.0, ase2) // #296 
