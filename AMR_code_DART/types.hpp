#pragma once

#include <Eigen/Core>
#include "utils.cpp"

// Contains the state of the LIP robot
struct State {
    Eigen::Vector3d comPos;
    Eigen::Vector3d comVel;
    Eigen::Vector3d comAcc;
    Eigen::Vector3d zmpPos;

    Eigen::Vector3d leftBackFootPos;
    Eigen::Vector3d leftBackFootVel;
    Eigen::Vector3d leftBackFootAcc;

    Eigen::Vector3d rightBackFootPos;
    Eigen::Vector3d rightBackFootVel;
    Eigen::Vector3d rightBackFootAcc;

    Eigen::Vector3d leftFrontFootPos;
    Eigen::Vector3d leftFrontFootVel;
    Eigen::Vector3d leftFrontFootAcc;

    Eigen::Vector3d rightFrontFootPos;
    Eigen::Vector3d rightFrontFootVel;
    Eigen::Vector3d rightFrontFootAcc;

    Eigen::Vector3d torsoOrient;
    Eigen::Vector3d leftBackFootOrient;
    Eigen::Vector3d rightBackFootOrient;
    Eigen::Vector3d leftFrontFootOrient;
    Eigen::Vector3d rightFrontFootOrient;


    inline Eigen::VectorXd getComPose() {
	Eigen::VectorXd comPose(6);
        comPose << torsoOrient, comPos;
        return comPose;
    }

    inline Eigen::VectorXd getSupportFootPose(bool supportFoot) {
	Eigen::VectorXd sfPose(6);
        if (supportFoot == 0) sfPose << leftBackFootOrient, leftBackFootPos;
        else sfPose << rightBackFootOrient, rightBackFootPos;
        return sfPose;
    }

    inline Eigen::VectorXd getFrontSwingFootPose(bool supportFoot) {
	Eigen::VectorXd sfPose(6);
        if (supportFoot == 1) sfPose << rightFrontFootOrient, rightFrontFootPos;
        else sfPose << leftFrontFootOrient, leftFrontFootPos;
        return sfPose;
    }

    inline Eigen::VectorXd getBackSwingFootPose(bool supportFoot) {
	Eigen::VectorXd sfPose(6);
        if (supportFoot == 1) sfPose << leftBackFootOrient, leftBackFootPos;
        else sfPose << rightBackFootOrient, rightBackFootPos;
        return sfPose;
    }

    inline Eigen::VectorXd getRelComPose(bool supportFoot) {
	return vvRel(getComPose(), getSupportFootPose(supportFoot));
    }

    inline Eigen::VectorXd getRelFrontSwingFootPose(bool supportFoot) {
	return vvRel(getFrontSwingFootPose(supportFoot), getSupportFootPose(supportFoot));
    }

    inline Eigen::VectorXd getRelBackSwingFootPose(bool supportFoot) {
	return vvRel(getBackSwingFootPose(supportFoot), getSupportFootPose(supportFoot));
    }

};

struct WalkState {
    bool supportFoot;
    double simulationTime;
    int mpcIter, controlIter, footstepCounter, indInitial;
};

struct Vref {
    Vref(double _x, double _y, double _omega) : x(_x), y(_y), omega(_omega) {}

    double x = 0;
    double y = 0;
    double omega = 0;
}; 
