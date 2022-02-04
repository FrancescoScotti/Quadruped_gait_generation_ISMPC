#pragma once

#include "./Eigen/Core"
#include "./dart/dart.hpp"
#include "MPCSolver.hpp"
#include "Utility.hpp"
#include  <fstream>
#include "utils.cpp"
#include "types.hpp"
//#include "FootstepPlan.hpp"
#include "StateFiltering.hpp"
#include <memory>
#include <random>

class Controller
{
public:
  Controller(dart::dynamics::SkeletonPtr _robot, dart::simulation::WorldPtr _world);
  virtual ~Controller();

  Eigen::Vector3d getZmpFromExternalForces();

  void update();

  Eigen::MatrixXd getTorsoAndSwfJacobian();
  Eigen::MatrixXd getTorsoAndSwfJacobian_frontRight();
  Eigen::MatrixXd getTorsoAndSwfJacobian_frontLeft();
  Eigen::MatrixXd getTorsoAndSwfJacobian_backRight();
  Eigen::MatrixXd getTorsoAndSwfJacobian_backLeft();
  Eigen::MatrixXd getTorsoAndSwfJacobianDeriv();

  Eigen::VectorXd getJointVelocitiesQp(State current, State desired);
  Eigen::VectorXd getJointVelocitiesQpAcceleration(State current, State desired);

  void setInitialConfiguration();
  void ArmSwing();
  void AnkleRegulation();
  void FixWaistToChestJoints();

  Eigen::VectorXd getJointVelocitiesStacked(Eigen::VectorXd, Eigen::VectorXd,
		  Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd);
  Eigen::VectorXd getJointVelocitiesStacked_frontRight(Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd);
  Eigen::VectorXd getJointVelocitiesStacked_frontLeft(Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd);

  Eigen::VectorXd getJointVelocitiesStacked_backRight(Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd);
  Eigen::VectorXd getJointVelocitiesStacked_backLeft(Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd);

  Eigen::Vector3d getRPY(dart::dynamics::BodyNode*, dart::dynamics::BodyNode*);
  Eigen::Vector3d getRPY(dart::dynamics::BodyNode*);

  void storeData();

  virtual void keyboard(unsigned char _key, int _x, int _y);

private:

  dart::dynamics::SkeletonPtr mRobot;

  dart::dynamics::BodyNode* mTorso;

  dart::simulation::WorldPtr mWorld;

  bool balancePoint;
  bool TORSO = false;
  bool COM = true;

  Eigen::MatrixXd ftsp_and_time;

  MPCSolver* solver;

  dart::dynamics::BodyNode* mBackLeftFoot;  
  dart::dynamics::BodyNode* mFrontLeftFoot;
  dart::dynamics::BodyNode* mBackRightFoot;
  dart::dynamics::BodyNode* mFrontRightFoot; 
  dart::dynamics::BodyNode* mSupportFoot;
  dart::dynamics::BodyNode* mSupportFoot_front;
  dart::dynamics::BodyNode* mSwingFoot_front;  
  dart::dynamics::BodyNode* mSwingFoot_back;
  dart::dynamics::BodyNode* mBase;

  State desired;
  State current;
  WalkState walkState;
  
  // our decoaretions 
  
  int file_size_com;
  int file_size_footsteps;
  
  Eigen::VectorXd x_des;
  Eigen::VectorXd y_des;
  Eigen::VectorXd z_des;

  Eigen::VectorXd xd_des;
  Eigen::VectorXd yd_des;
  Eigen::VectorXd zd_des;

  Eigen::VectorXd xFootDesired_rf;
  Eigen::VectorXd xFootDesired_lf;
  Eigen::VectorXd xFootDesired_rr;
  Eigen::VectorXd xFootDesired_lr;
  
  Eigen::VectorXd yFootDesired_rf;
  Eigen::VectorXd yFootDesired_lf;
  Eigen::VectorXd yFootDesired_rr;
  Eigen::VectorXd yFootDesired_lr;
  
  Eigen::VectorXd zFootDesired_rf;
  Eigen::VectorXd zFootDesired_lf;
  Eigen::VectorXd zFootDesired_rr;
  Eigen::VectorXd zFootDesired_lr;
  

  Eigen::Vector3d previousCoM, currentCoM;
  void createForceLine();
  const Eigen::Vector4d DefaultForceLineColor
    = Eigen::Vector4d(1.0, 0.0, 0.0, 1.0);
  dart::dynamics::LineSegmentShapePtr mForceLine;

};
