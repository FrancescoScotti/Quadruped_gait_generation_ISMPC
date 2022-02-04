#pragma once

#include <Eigen/Core>
#include "qpOASES/qpOASES.hpp"
#include "types.hpp"
#include "parameters.cpp"
#include "utils.cpp"
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <math.h>       /* isnan, sqrt */
//#include "FootstepPlan.hpp"
#include <chrono>

class MPCSolver{
    public:
    MPCSolver(const Eigen::MatrixXd& ftsp_and_timings);
    ~MPCSolver();

    // Compute the next desired state starting from the current state
    State solve(State current, WalkState walkState, const Eigen::MatrixXd& ftsp_and_timings);

    // some stuff
    int itr;
    int fsCount, old_fsCount, adaptation_memo, ds_samples, ct;

    double xz_dot, yz_dot;

    private:
    // Matrices for prediction
    Eigen::VectorXd p;
    Eigen::MatrixXd P;

    // Matrices for cost function
    Eigen::MatrixXd costFunctionH;
    Eigen::VectorXd costFunctionF;

    // Matrices for stability constraint
    Eigen::MatrixXd Aeq;
    Eigen::VectorXd beq;

    //Matrices for balance constraint
    Eigen::MatrixXd AZmp;
    Eigen::VectorXd bZmpMax;
    Eigen::VectorXd bZmpMin;

    // Matrices for kinematic constraints
    Eigen::MatrixXd AFootsteps;
    Eigen::VectorXd bFootstepsMax;
    Eigen::VectorXd bFootstepsMin;

    // Matrices for the stacked constraints
    Eigen::MatrixXd AConstraint;
    Eigen::VectorXd bConstraintMax;
    Eigen::VectorXd bConstraintMin;

    // Matrices for Z QP and constraints
    Eigen::MatrixXd costFunctionH_z, costFunctionH_xy;
    Eigen::VectorXd costFunctionF_z, costFunctionF_x, costFunctionF_y;
    Eigen::MatrixXd Aeq_z;
    Eigen::VectorXd beq_z;
    Eigen::MatrixXd Aineq_z;
    Eigen::VectorXd bMax_z;
    Eigen::VectorXd bMin_z;
    Eigen::MatrixXd Af_z;
    Eigen::VectorXd fMax_z;
    Eigen::VectorXd fMin_z;
    Eigen::MatrixXd S_bar_z, S_bar_z_v;
    Eigen::MatrixXd T_bar_z, T_bar_z_v;
    Eigen::MatrixXd S_bar_g_z, S_bar_g_z_v;
    Eigen::MatrixXd T_bar_g_z, T_bar_g_z_v;
    Eigen::VectorXd lambda, etas;
    Eigen::MatrixXd deltas;
    Eigen::VectorXd z_trajectory;
    Eigen::VectorXd f_trajectory;
    Eigen::MatrixXd A_virtual_ZMP;
    Eigen::VectorXd b_virtual_ZMP_x, b_virtual_ZMP_y;
    Eigen::MatrixXd S_p, T_p;
    Eigen::MatrixXd S_v, T_v;
    Eigen::MatrixXd phi_input, phi_state, C_sc;

    //Matrices State update
    Eigen::Matrix2d A_xy;    
    Eigen::MatrixXd B_xy;    
    Eigen::MatrixXd A_z;    
    Eigen::MatrixXd B_z, Bg_z, C_p, C_v;
    Eigen::MatrixXd C_z;
    
    // Midpoint of ZMP constraint
    Eigen::MatrixXd ftsp_midpoint; 

    // useful parameters
    double h_des = comTargetHeight; //0.69
    double h_ub = 0.015;
    double h_lb = 0.025;
    double total_flight_phase = 0.02;
    double a_f;

};
