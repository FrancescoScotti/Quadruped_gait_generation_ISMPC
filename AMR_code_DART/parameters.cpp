#pragma once

#include <math.h>

// Definable parameters
// ********************

// Times
const double mpcTimeStep = 0.01; //0.05;
const double controlTimeStep = 0.01;
const double singleSupportDuration = 0.35; //0.3
const double doubleSupportDuration = 0.1; // 0.2
const double predictionTime = 1.0;

// Walk parameters
const double stepHeight = 0.033; //0.03
const double comTargetHeight = 0.69;
const double kSwingFoot = 0.05; 

// Constraints
const double thetaMax = 0.30;
const double footConstraintSquareWidth = 0.09;
const double deltaXMax = 0.25;
const double deltaYIn = 0.15;
const double deltaYOut = 0.28;

// Cost function weights for horizontal QP
const double qZ = 1;
const double qF = 1000000;
// Cost function weights for vertical QP
const double q_force = 1.0;
const double q_position = 1000000000000.0;

// Kinematic control
const double IKerrorGain = 1.0; 

// Used in the code
// ****************
const double mass_hrp4 = 50.0; 
const double g = 9.81; 
const double eta = sqrt(g/comTargetHeight);
const int N = round(predictionTime/mpcTimeStep);
const int S = round(singleSupportDuration/mpcTimeStep);
const int F = round(doubleSupportDuration/mpcTimeStep);
const int M = 2; 


