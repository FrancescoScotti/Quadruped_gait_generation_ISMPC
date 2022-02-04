#include "Controller.hpp"
#include "parameters.cpp"

using namespace std;
using namespace dart::dynamics;
using namespace dart::simulation;

static std::ofstream fout1(realpath("../data/currentFoot.txt", NULL), std::ofstream::trunc);
static std::ofstream fout2(realpath("../data/currentCom.txt", NULL), std::ofstream::trunc);
static std::ofstream fout3(realpath("../data/desiredFoot.txt", NULL), std::ofstream::trunc);
static std::ofstream fout4(realpath("../data/desiredCom.txt", NULL), std::ofstream::trunc);
static std::ofstream fout5(realpath("../data/currentComAcc.txt", NULL), std::ofstream::trunc);
static std::ofstream fout6(realpath("../data/desiredComAcc.txt", NULL), std::ofstream::trunc);
static std::ofstream fout7(realpath("../data/jacobianComAcc.txt", NULL), std::ofstream::trunc);
static std::ofstream fout8(realpath("../data/measuredZmp.txt", NULL), std::ofstream::trunc);
static std::ofstream fout9(realpath("../data/desiredZmp.txt", NULL), std::ofstream::trunc);


static std::ofstream foutDebug(realpath("../data/debug.txt", NULL), std::ofstream::trunc);

    

Controller::Controller(dart::dynamics::SkeletonPtr _robot, dart::simulation::WorldPtr _world)
	: mRobot(_robot), mWorld(_world)
{
    // some useful pointers to robot limbs
    /*
    mLeftFoot = mRobot->getBodyNode("LF_HIP");
    mRightFoot = mRobot->getBodyNode("LF_HIP");
    mBase = mRobot->getBodyNode("base");
    mTorso = mRobot->getBodyNode("base");
    /**/

    // INITIAL CONFIGURATION - top view of the robot in current configuration

    //     --  ||1)           (2||  --
    //         ||               || 
    //     --  ||3)           (4||  --

    //     <-------------  x axis
    //     |
    //     V               y axis

    // 1)
    mFrontRightFoot = mRobot->getBodyNode("RF_FOOT");

    // 3)
    mFrontLeftFoot = mRobot->getBodyNode("LF_FOOT");

    // 2)
    mBackRightFoot = mRobot->getBodyNode("RH_FOOT");

    // 4)
    mBackLeftFoot = mRobot->getBodyNode("LH_FOOT");


    ///////////////////////////////////////////////


    mBase = mRobot->getBodyNode("base");
    mTorso = mRobot->getBodyNode("base_inertia");


    // Initialize walk state
    walkState.mpcIter = 0;
    walkState.controlIter = 0;
    walkState.footstepCounter = 0;
    walkState.supportFoot = true;

    balancePoint = TORSO;

    Eigen::Vector3d comInitialPosition;
    if (balancePoint == TORSO) {
        comInitialPosition = mBase->getCOM();
    } else {
        comInitialPosition = mRobot->getCOM();
    }

    // load desired footsteps

    std::ifstream inFile;
    std::string myText;

  
    //Read the footsteps
    
    // Desired footsteps and timings in world frame

    int N_footsteps = 40;
    ftsp_and_time = Eigen::MatrixXd::Zero(N_footsteps,4);
    int start = 1;
    for (int i = start; i < N_footsteps; i++) {
         ftsp_and_time(i,0) = (i-start)*0.2; //0.09  0.1  0.13 0.15  0.18
         ftsp_and_time(i,1) = pow(-(double)1,(i-start))*0.08;  //0.08
         ftsp_and_time(i,2) = 0.0;
         ftsp_and_time(i,3) = (double)((double)mpcTimeStep/(double)controlTimeStep)*(S+F)*i;  // time: BAD IMPLEMENTATION
    }
 

    bool firstSupportFootIsLeft = true;
    Eigen::VectorXd leftFootPose(6), rightFootPose(6);

    
    // Instantiate MPC solver
    const Eigen::MatrixXd& ftsp_and_time_ref = ftsp_and_time;
    solver = new MPCSolver(ftsp_and_time_ref);
    current.comPos = balancePoint == TORSO ? mBase->getCOM() : mRobot->getCOM();
    current.comVel = balancePoint == TORSO ? mBase->getCOMLinearVelocity() : mRobot->getCOMLinearVelocity();
    current.comAcc = balancePoint == TORSO ? mBase->getCOMLinearAcceleration() : mRobot->getCOMLinearAcceleration();
    desired.comPos = Eigen::Vector3d(current.comPos(0), current.comPos(1), comTargetHeight);
    desired.comVel = Eigen::Vector3d::Zero();
    desired.zmpPos = Eigen::Vector3d(current.comPos(0), current.comPos(1),0.0);

    desired.leftBackFootPos = mBackLeftFoot->getCOM();
    desired.rightBackFootPos = mBackRightFoot->getCOM();
    desired.leftBackFootOrient = getRPY(mBackLeftFoot);
    desired.rightBackFootOrient = getRPY(mBackRightFoot);

    desired.leftFrontFootPos = mFrontLeftFoot->getCOM();
    desired.rightFrontFootPos = mFrontRightFoot->getCOM();
    desired.leftFrontFootOrient = getRPY(mFrontLeftFoot);
    desired.rightFrontFootOrient = getRPY(mFrontRightFoot);

    desired.torsoOrient = Eigen::Vector3d(0.0, 0.0, 0.0);

    desired.comAcc = eta * eta * (desired.comPos - desired.zmpPos);

    previousCoM = current.comPos;
    
    ///// Compute or load the desired CoM velocity and position ///////////////
    
    int file_size_footsteps = 2000;

	xFootDesired_rf = Eigen::VectorXd::Zero(file_size_footsteps);
    xFootDesired_lf = Eigen::VectorXd::Zero(file_size_footsteps);
    xFootDesired_rr = Eigen::VectorXd::Zero(file_size_footsteps);
    xFootDesired_lr = Eigen::VectorXd::Zero(file_size_footsteps);
    yFootDesired_rf = Eigen::VectorXd::Zero(file_size_footsteps);
    yFootDesired_lf = Eigen::VectorXd::Zero(file_size_footsteps);
    yFootDesired_rr = Eigen::VectorXd::Zero(file_size_footsteps);
    yFootDesired_lr = Eigen::VectorXd::Zero(file_size_footsteps);
    zFootDesired_rf = Eigen::VectorXd::Zero(file_size_footsteps);
    zFootDesired_lf = Eigen::VectorXd::Zero(file_size_footsteps);
    zFootDesired_rr = Eigen::VectorXd::Zero(file_size_footsteps);
    zFootDesired_lr = Eigen::VectorXd::Zero(file_size_footsteps);
    int i = 0;
    //inFile.open("../MATLAB/trotting/phipi4/foot_rl_trot_phipi4.txt");
    inFile.open("../MATLAB/walking/phi0_10cm_50/foot_rl_walk_phi0.txt");
    float foot_lr_x,foot_lr_y,foot_lr_z;
    while (getline (inFile, myText)) {
        // Output the text from the file
        if (3 == std::sscanf(myText.c_str(), "%f %f %f\n", &foot_lr_x, &foot_lr_y, &foot_lr_z))
        {
            xFootDesired_lr(i) = foot_lr_x;
            yFootDesired_lr(i) = foot_lr_y;
            zFootDesired_lr(i) = foot_lr_z;
            i++;
        }
    }
    inFile.close();

    i = 0;
    //inFile.open("../MATLAB/trotting/phipi4/foot_rr_trot_phipi4.txt");
    inFile.open("../MATLAB/walking/phi0_10cm_50/foot_rr_walk_phi0.txt");
    float foot_rr_x,foot_rr_y,foot_rr_z;
    while (getline (inFile, myText)) {
        // Output the text from the file
        if (3 == std::sscanf(myText.c_str(), "%f %f %f\n", &foot_rr_x, &foot_rr_y, &foot_rr_z))
        {
            xFootDesired_rr(i) = foot_rr_x;
            yFootDesired_rr(i) = foot_rr_y;
            zFootDesired_rr(i) = foot_rr_z;
            i++;
        }
    }
    inFile.close();


    i = 0;
    //inFile.open("../MATLAB/trotting/phipi4/foot_fl_trot_phipi4.txt");
    inFile.open("../MATLAB/walking/phi0_10cm_50/foot_fl_walk_phi0.txt");
    float foot_lf_x,foot_lf_y,foot_lf_z;
    while (getline (inFile, myText)) {
        // Output the text from the file
        if (3 == std::sscanf(myText.c_str(), "%f %f %f\n", &foot_lf_x, &foot_lf_y, &foot_lf_z))
        {
            xFootDesired_lf(i) = foot_lf_x;
            yFootDesired_lf(i) = foot_lf_y;
            zFootDesired_lf(i) = foot_lf_z;
            i++;
        }
    }
    inFile.close();


    i = 0;
    //inFile.open("../MATLAB/trotting/phipi4/foot_fr_trot_phipi4.txt");
    inFile.open("../MATLAB/walking/phi0_10cm_50/foot_fr_walk_phi0.txt");
    float foot_rf_x,foot_rf_y,foot_rf_z;
    while (getline (inFile, myText)) {
        // Output the text from the file
        if (3 == std::sscanf(myText.c_str(), "%f %f %f\n", &foot_rf_x, &foot_rf_y, &foot_rf_z))
        {
            xFootDesired_rf(i) = foot_rf_x;
            yFootDesired_rf(i) = foot_rf_y;
            zFootDesired_rf(i) = foot_rf_z;
            i++;
        }
    }
    inFile.close();





    i = 0;
    //inFile.open("../MATLAB/trotting/phipi4/ComTrajectory_trot_phipi4.txt");
    inFile.open("../MATLAB/walking/phi0_10cm_50/ComTrajectory_walk_phi0.txt");
    if (!inFile) {
          std::cout << "Unable to open File" << std::endl;
          exit(1);
    }

    int file_size_com = 0;
    while (getline(inFile, myText)) {
         file_size_com = file_size_com + 1;
    }
    inFile.close();
    float x;
    float y;
    float z;
    x_des = Eigen::VectorXd::Zero(file_size_com);
    y_des = Eigen::VectorXd::Zero(file_size_com);
    z_des = Eigen::VectorXd::Zero(file_size_com);
    //inFile.open("../MATLAB/trotting/phipi4/ComTrajectory_trot_phipi4.txt");
    inFile.open("../MATLAB/walking/phi0_10cm_50/ComTrajectory_walk_phi0.txt");
    while (getline (inFile, myText)) {
        // Output the text from the file
        if (3 == std::sscanf(myText.c_str(), "%f %f %f\n", &x, &y, &z))
        {
            x_des(i) = x;
            y_des(i) = y;
            z_des(i) = z;
            i++;
        }
    }
    inFile.close();


    i = 0;
    //inFile.open("../MATLAB/trotting/phipi4/ComVelocity_trot_phipi4.txt");
    inFile.open("../MATLAB/walking/phi0_10cm_50/ComVelocity_walk_phi0.txt");
    if (!inFile) {
          std::cout << "Unable to open File" << std::endl;
          exit(1);
    }

    file_size_com = 0;
    while (getline(inFile, myText)) {
         file_size_com = file_size_com + 1;
    }
    inFile.close();
    float xd;
    float yd;
    float zd;
    xd_des = Eigen::VectorXd::Zero(file_size_com);
    yd_des = Eigen::VectorXd::Zero(file_size_com);
    zd_des = Eigen::VectorXd::Zero(file_size_com);
    //inFile.open("../MATLAB/trotting/phipi4/ComVelocity_trot_phipi4.txt");
    inFile.open("../MATLAB/walking/phi0_10cm_50/ComVelocity_walk_phi0.txt");
    while (getline (inFile, myText)) {
        // Output the text from the file
        if (3 == std::sscanf(myText.c_str(), "%f %f %f\n", &xd, &yd, &zd))
        {
            xd_des(i) = xd;
            yd_des(i) = yd;
            zd_des(i) = zd;
            i++;
        }
    }
    inFile.close();

}

Controller::~Controller() {}



void Controller::update() {

    if (((int) walkState.simulationTime) % 5 == 0 ) {

        previousCoM = current.comPos;

    }

    if (walkState.simulationTime >= ftsp_and_time(walkState.footstepCounter,3) -1  && false) {
    	walkState.controlIter = 0;
    	walkState.mpcIter = 0;
        walkState.footstepCounter++;
	walkState.supportFoot = !walkState.supportFoot;
        //std::cout << "Iteration " << walkState.controlIter << " Footstep " << walkState.footstepCounter << std::endl;
        std::cout <<  mWorld->getSimFrames() << std::endl;
    }


    // STOP EXECUTION AFTER 10 SECONDS
    if (mWorld->getSimFrames() == 3000) exit(1);
    // simulation time
    walkState.simulationTime = mWorld->getSimFrames();
    // Retrieve current position and orientation of useful quantities in the world frame
    current.comPos = balancePoint == TORSO ? mBase->getCOM() : mRobot->getCOM();
    current.comVel = balancePoint == TORSO ? mBase->getCOMLinearVelocity() : mRobot->getCOMLinearVelocity();
    current.comAcc = balancePoint == TORSO ? mBase->getCOMLinearAcceleration() : mRobot->getCOMLinearAcceleration();
    /*
    current.comPos = mTorso->getCOM();
    current.comVel = mTorso->getCOMLinearVelocity();
    current.comAcc = mTorso->getCOMLinearAcceleration();
    /**/
    current.zmpPos = current.comPos - current.comAcc / (eta*eta);

    current.leftBackFootPos = mBackLeftFoot->getCOM();
    current.rightBackFootPos = mBackRightFoot->getCOM();
    current.leftBackFootOrient = getRPY(mBackLeftFoot);
    current.rightBackFootOrient = getRPY(mBackRightFoot);

    current.leftFrontFootPos = mFrontLeftFoot->getCOM();
    current.rightFrontFootPos = mFrontRightFoot->getCOM();
    current.leftFrontFootOrient = getRPY(mFrontLeftFoot);
    current.rightFrontFootOrient = getRPY(mFrontRightFoot);

    current.torsoOrient = getRPY(mTorso);

    if (walkState.supportFoot == 0 ) {
    //std::cout << "current.rightFootPos " << current.rightFootPos(0) << std::endl;
    //std::cout << "desired.rightFootPos " << desired.rightFootPos(0) << std::endl;
    } else {
    //std::cout << "current.leftFootPos " << current.leftFootPos(0) << std::endl;
    //std::cout << "desired.leftFootPos " << desired.leftFootPos(0) << std::endl;
    }
 
    // Filtered Desired State for quasi-closed loop behavior


/**/
    // Compute the new desired state using MPC: create a constant reference to the candidate footsteps and timings to be passed in the MPCSolver
    //const Eigen::MatrixXd& ftsp_and_time_ref = ftsp_and_time;
    //desired = solver->solve(desired, walkState, ftsp_and_time_ref);

    //desired.comPos << 0.08*cos(2*3.14*0.01*(mWorld->getSimFrames())/2.0),0*0.06*sin(2*3.14*0.01*(mWorld->getSimFrames())/2.0),0.0;
 
    walkState.supportFoot = 1;

/*
    // Compute inverse kinematics
     Eigen::VectorXd qDot =  getJointVelocitiesStacked(desired.getRelComPose(walkState.supportFoot),  current.getRelComPose(walkState.supportFoot),
  	 	 desired.getRelFrontSwingFootPose(walkState.supportFoot),  current.getRelFrontSwingFootPose(walkState.supportFoot), desired.getRelBackSwingFootPose(walkState.supportFoot),  current.getRelBackSwingFootPose(walkState.supportFoot), desired.comVel);
/**/

/*
    Eigen::VectorXd desired_Com_pose = Eigen::VectorXd::Zero(6);
    Eigen::VectorXd current_Com_pose = Eigen::VectorXd::Zero(6);
    Eigen::VectorXd desired_FrontSwingFootPose = Eigen::VectorXd::Zero(6);
    Eigen::VectorXd current_FrontSwingFootPose = Eigen::VectorXd::Zero(6);
    Eigen::VectorXd desired_BackSwingFootPose = Eigen::VectorXd::Zero(6);
    Eigen::VectorXd current_BackSwingFootPose = Eigen::VectorXd::Zero(6);


    current_Com_pose << current.torsoOrient, current.comPos;
    desired_Com_pose << desired.torsoOrient, desired.comPos;

    current_FrontSwingFootPose << current.rightFrontFootOrient, current.rightFrontFootPos;
    desired_FrontSwingFootPose << desired.rightFrontFootOrient, desired.rightFrontFootPos;

    current_BackSwingFootPose << current.leftBackFootOrient, current.leftBackFootPos;
    desired_BackSwingFootPose << desired.leftBackFootOrient, desired.leftBackFootPos;

    Eigen::VectorXd qDot =  getJointVelocitiesStacked(desired_Com_pose,  current_Com_pose,
  	 	 desired_FrontSwingFootPose,  current_FrontSwingFootPose, desired_BackSwingFootPose,  current_BackSwingFootPose, desired.comVel);

/**/

    // Change frame for inverse kinematics
    //Eigen::VectorXd qDot =  getJointVelocitiesStacked_worldframe(desired.getComPose(),  current.getComPose(),desired.getSwingFootPose(walkState.supportFoot),  current.getSwingFootPose(walkState.supportFoot),  desired.comVel);

    // THEY MUST BE 12 i.e. 4 DOF for each joint!!!
    //qDot = Eigen::VectorXd::Zero(12);

    

    desired.comPos(0) = x_des(mWorld->getSimFrames());
    desired.comPos(1) = y_des(mWorld->getSimFrames());
    desired.comPos(2) = z_des(mWorld->getSimFrames());

/*   TEST MOVE LEG
    desired.comPos(0) = -0.1;
    desired.comPos(1) = 0.0 ;
    desired.comPos(2) = 0.56;
/**/

    Eigen::VectorXd desired_com_vel = Eigen::VectorXd::Zero(6);
    desired.comVel(0) = xd_des(mWorld->getSimFrames());
    desired.comVel(1) = yd_des(mWorld->getSimFrames());
    desired.comVel(2) = zd_des(mWorld->getSimFrames());
    
    desired_com_vel << Eigen::Vector3d::Zero(3), desired.comVel;
    
    // Compute or load the desired velocity and position for each leg
    // given by hand
    desired.rightFrontFootPos << xFootDesired_rf(mWorld->getSimFrames()),yFootDesired_rf(mWorld->getSimFrames()),zFootDesired_rf(mWorld->getSimFrames());
    desired.leftFrontFootPos <<  xFootDesired_lf(mWorld->getSimFrames()),yFootDesired_lf(mWorld->getSimFrames()),zFootDesired_lf(mWorld->getSimFrames());
    desired.rightBackFootPos <<  xFootDesired_rr(mWorld->getSimFrames()),yFootDesired_rr(mWorld->getSimFrames()),zFootDesired_rr(mWorld->getSimFrames());
    desired.leftBackFootPos <<  xFootDesired_lr(mWorld->getSimFrames()),yFootDesired_lr(mWorld->getSimFrames()),zFootDesired_lr(mWorld->getSimFrames());

    // desired velocities legs w.r.t. CoM
    Eigen::VectorXd desired_vel_swing = Eigen::VectorXd::Zero(6);
    desired_vel_swing = -desired_com_vel + desired_vel_swing; // move velocity frame (....)

    /// Task vectors for each leg ///

    Eigen::VectorXd desired_rightFrontFootPose = Eigen::VectorXd::Zero(6);
    Eigen::VectorXd current_rightFrontFootPose = Eigen::VectorXd::Zero(6);
    // CoM and Swing foot positions in CoM Frame
    //desired.rightFrontFootPos(2) = 0.05 + 0.05*sin(2.0*3.14*0.01*(mWorld->getSimFrames())/3.0); TEST MOVE LEG
    current_rightFrontFootPose << current.rightFrontFootOrient, -current.comPos + current.rightFrontFootPos;
    desired_rightFrontFootPose << desired.rightFrontFootOrient, -desired.comPos + desired.rightFrontFootPos; //-desired.comPos + desired.rightFrontFootPos;

    Eigen::VectorXd desired_leftFrontFootPose = Eigen::VectorXd::Zero(6);
    Eigen::VectorXd current_leftFrontFootPose = Eigen::VectorXd::Zero(6);
    // CoM and Swing foot positions in CoM Frame
    current_leftFrontFootPose << current.leftFrontFootOrient, -current.comPos + current.leftFrontFootPos;
    desired_leftFrontFootPose << desired.leftFrontFootOrient, -desired.comPos + desired.leftFrontFootPos; //-desired.comPos + desired.rightFrontFootPos;

    Eigen::VectorXd desired_rightBackFootPose = Eigen::VectorXd::Zero(6);
    Eigen::VectorXd current_rightBackFootPose = Eigen::VectorXd::Zero(6);
    // CoM and Swing foot positions in CoM Frame
    current_rightBackFootPose << current.rightBackFootOrient, -current.comPos + current.rightBackFootPos;
    desired_rightBackFootPose << desired.rightBackFootOrient, -desired.comPos + desired.rightBackFootPos; //-desired.comPos + desired.rightFrontFootPos;


    Eigen::VectorXd desired_leftBackFootPose = Eigen::VectorXd::Zero(6);
    Eigen::VectorXd current_leftBackFootPose = Eigen::VectorXd::Zero(6);
    // CoM and Swing foot positions in CoM Frame
    current_leftBackFootPose << current.leftBackFootOrient, -current.comPos + current.leftBackFootPos;
    desired_leftBackFootPose << desired.leftBackFootOrient, -desired.comPos + desired.leftBackFootPos; //-desired.comPos + desired.rightFrontFootPos;

    // ........  //

    // THIS IS THE WAY, maybe...

    // desired CoM is given in the world frame
    // desired CoM Velocity is given in the world frame

    // desired Swing is given in the world frame
    // desired Swing Velocity is diven in the world frame





    // Conceptual height problem (?)   --- it moves the leg but not the CoM

    // -- maybe the other legs are too fixed and the foot simply slips away



    //Eigen::VectorXd qDot =  getJointVelocitiesStacked_frontRight(desired_rightFrontFootPose,  current_rightFrontFootPose, desired_vel_swing);

    Eigen::VectorXd qDot = getJointVelocitiesStacked_frontRight(desired_rightFrontFootPose,  current_rightFrontFootPose, desired_vel_swing) +
                           getJointVelocitiesStacked_frontLeft(desired_leftFrontFootPose,  current_leftFrontFootPose, desired_vel_swing) + 
                           getJointVelocitiesStacked_backRight(desired_rightBackFootPose,  current_rightBackFootPose, desired_vel_swing) +  
                           getJointVelocitiesStacked_backLeft(desired_leftBackFootPose,  current_leftBackFootPose, desired_vel_swing);


    //std::cout << qDot << std::endl;   
    /*std::cout << "desired.comPos(0) " << desired.comPos(0) << " current.comPos(0) " << current.comPos(0) << " frame " <<  (mWorld->getSimFrames())  << std::endl;
    std::cout << "desired.comPos(1) " << desired.comPos(1) << " current.comPos(1) " << current.comPos(1) << " frame " <<  (mWorld->getSimFrames())  << std::endl;*/
    std::cout << "error_x" << desired.comPos(0)-current.comPos(0) << " frame " <<  (mWorld->getSimFrames())  << std::endl;
    std::cout << "error_y" << desired.comPos(1)-current.comPos(1) << " frame " <<  (mWorld->getSimFrames())  << std::endl;

    //int stop = getchar();

    // ........  //



    // Set the velocity of the floating base to zero
    for (int i = 0; i < 6; ++i){
        mRobot->setCommand(i, 0.0);
    }
    // Set the velocity of each joint from the computed inverse kinematics
    for (int i = 0; i < 12; ++i){
        mRobot->setCommand(i+6,qDot(i)); //velocity joint control
    }

/**/
    // Store the results in files (for plotting)
    //storeData();
    // Arm swing


    // Update the iteration counters, if a step is finished reset and change support foot
    ++walkState.controlIter;
    walkState.mpcIter = floor(walkState.controlIter*controlTimeStep/mpcTimeStep);
    // When a footstep has been completed

    //int stop = getchar();
    if (mWorld->getSimFrames() >= 1) createForceLine();  

}



void Controller::AnkleRegulation() {

}



Eigen::MatrixXd Controller::getTorsoAndSwfJacobian() { //FIXME this is still expressed in relative frame
    if (walkState.supportFoot == true) {  //walkState.supportFoot == true
        mSupportFoot = mBackRightFoot;
        mSupportFoot_front = mFrontLeftFoot;
        // FOR NOW: Swing LF_FOOT and LH_FOOT !!!!!!
        mSwingFoot_front = mFrontRightFoot;  // DEFINE GAIT !!!!!
        mSwingFoot_back = mBackLeftFoot;  // DEFINE GAIT !!!!!
    } else {
        mSupportFoot = mRobot->getBodyNode("LF_FOOT");
        mSupportFoot_front = mRobot->getBodyNode("RH_FOOT");

        mSwingFoot_front = mRobot->getBodyNode("LH_FOOT");  // DEFINE GAIT !!!!!
        mSwingFoot_back = mRobot->getBodyNode("RF_FOOT");  // DEFINE GAIT !!!!!
    }

    Eigen::MatrixXd Jacobian_supportToBase;
    Eigen::MatrixXd Jacobian_supportToBase_front;

    /*

    if (balancePoint == TORSO) {
	Jacobian_supportToBase =  mRobot->getJacobian(mTorso,mSupportFoot) - mRobot->getJacobian(mSupportFoot,mSupportFoot); //(mBase,mSupportFoot)
    } else {
	Jacobian_supportToBase =  (mRobot->getCOMJacobian(mSupportFoot) - mRobot->getJacobian(mSupportFoot,mSupportFoot));
    }

    for (unsigned int i=0; i<12; i++){ 
    // WHY THIS?????
    if (i>=4)  Jacobian_supportToBase.col(i).setZero();
    }

    Eigen::MatrixXd Jacobian_SupportToSwing_front =  mRobot->getJacobian(mSwingFoot_front,mSupportFoot) - mRobot->getJacobian(mSupportFoot,mSupportFoot);
    Eigen::MatrixXd Jacobian_SupportToSwing_back =  mRobot->getJacobian(mSwingFoot_back,mSupportFoot) - mRobot->getJacobian(mSupportFoot,mSupportFoot);

    /**/


    if (balancePoint == TORSO) {
	Jacobian_supportToBase =  mRobot->getJacobian(mSupportFoot,mTorso) - mRobot->getJacobian(mTorso,mTorso); //(mBase,mSupportFoot)
	Jacobian_supportToBase_front =  mRobot->getJacobian(mSupportFoot_front,mTorso) - mRobot->getJacobian(mTorso,mTorso); //(mBase,mSupportFoot)
    } else {
	Jacobian_supportToBase =  (mRobot->getCOMJacobian(mSupportFoot) - mRobot->getJacobian(mSupportFoot,mSupportFoot));
	Jacobian_supportToBase =  mRobot->getJacobian(mSupportFoot,mTorso) - mRobot->getJacobian(mTorso,mTorso); //(mBase,mSupportFoot)
	Jacobian_supportToBase_front =  mRobot->getJacobian(mSupportFoot_front,mTorso) - mRobot->getJacobian(mTorso,mTorso); //(mBase,mSupportFoot)
    }

    for (unsigned int i=0; i<12; i++){ 
    // WHY THIS?????
    //if (i>=4)  Jacobian_supportToBase.col(i).setZero();
    }

    Eigen::MatrixXd Jacobian_SupportToSwing_front =  mRobot->getJacobian(mSwingFoot_front,mTorso) - mRobot->getJacobian(mTorso,mTorso);
    Eigen::MatrixXd Jacobian_SupportToSwing_back =  mRobot->getJacobian(mSwingFoot_back,mTorso) - mRobot->getJacobian(mTorso,mTorso);

    /*
    std::cout << mRobot->getJacobian(mSupportFoot_front,mTorso) - mRobot->getJacobian(mTorso,mTorso) << std::endl;
    std::cout << "here" << std::endl;
    std::cout << Jacobian_supportToBase.rows() << " " << Jacobian_supportToBase.cols() << std::endl;
    std::cout << Jacobian_supportToBase_front.rows() << " " << Jacobian_supportToBase_front.cols() << std::endl;
    std::cout << Jacobian_SupportToSwing_front.rows() << " " << Jacobian_SupportToSwing_front.cols() << std::endl;
    std::cout << Jacobian_SupportToSwing_back.rows() << " " << Jacobian_SupportToSwing_back.cols() << std::endl;
    /**/
  
    Eigen::MatrixXd Jacobian_tot_(24, 18);

    Jacobian_tot_ << Jacobian_supportToBase, Jacobian_supportToBase_front, Jacobian_SupportToSwing_front, Jacobian_SupportToSwing_back;


    Eigen::MatrixXd Jacobian_tot(24, 12);

    // Remove the floating base columns
    Jacobian_tot = Jacobian_tot_.block<24,12>(0, 6);

    return Jacobian_tot;
}

Eigen::MatrixXd Controller::getTorsoAndSwfJacobian_frontRight() { //FIXME this is still expressed in relative frame

    mSwingFoot_front = mFrontRightFoot;  // DEFINE GAIT !!!!!
    Eigen::MatrixXd Jacobian_SupportToSwing_front =  mRobot->getJacobian(mSwingFoot_front,mTorso) - mRobot->getJacobian(mTorso,mTorso);

    Eigen::MatrixXd Jacobian_tot_(6, 18);

    Jacobian_tot_ << Jacobian_SupportToSwing_front;

    Eigen::MatrixXd Jacobian_tot(6, 12);

    // Remove the floating base columns
    Jacobian_tot = Jacobian_tot_.block<6,12>(0, 6);

    return Jacobian_tot;

}

Eigen::MatrixXd Controller::getTorsoAndSwfJacobian_frontLeft() { //FIXME this is still expressed in relative frame

    mSwingFoot_front = mFrontLeftFoot;  // DEFINE GAIT !!!!!
    Eigen::MatrixXd Jacobian_SupportToSwing_front =  mRobot->getJacobian(mSwingFoot_front,mTorso) - mRobot->getJacobian(mTorso,mTorso);

    Eigen::MatrixXd Jacobian_tot_(6, 18);

    Jacobian_tot_ << Jacobian_SupportToSwing_front;

    Eigen::MatrixXd Jacobian_tot(6, 12);

    // Remove the floating base columns
    Jacobian_tot = Jacobian_tot_.block<6,12>(0, 6);

    return Jacobian_tot;

}

Eigen::MatrixXd Controller::getTorsoAndSwfJacobian_backRight() { //FIXME this is still expressed in relative frame

    mSwingFoot_front = mBackRightFoot;  // DEFINE GAIT !!!!!
    Eigen::MatrixXd Jacobian_SupportToSwing_front =  mRobot->getJacobian(mSwingFoot_front,mTorso) - mRobot->getJacobian(mTorso,mTorso);

    Eigen::MatrixXd Jacobian_tot_(6, 18);

    Jacobian_tot_ << Jacobian_SupportToSwing_front;

    Eigen::MatrixXd Jacobian_tot(6, 12);

    // Remove the floating base columns
    Jacobian_tot = Jacobian_tot_.block<6,12>(0, 6);

    return Jacobian_tot;

}

Eigen::MatrixXd Controller::getTorsoAndSwfJacobian_backLeft() { //FIXME this is still expressed in relative frame

    mSwingFoot_front = mBackLeftFoot;  // DEFINE GAIT !!!!!
    Eigen::MatrixXd Jacobian_SupportToSwing_front =  mRobot->getJacobian(mSwingFoot_front,mTorso) - mRobot->getJacobian(mTorso,mTorso);

    Eigen::MatrixXd Jacobian_tot_(6, 18);

    Jacobian_tot_ << Jacobian_SupportToSwing_front;

    Eigen::MatrixXd Jacobian_tot(6, 12);

    // Remove the floating base columns
    Jacobian_tot = Jacobian_tot_.block<6,12>(0, 6);

    return Jacobian_tot;

}








Eigen::MatrixXd Controller::getTorsoAndSwfJacobianDeriv() {
/**
   if (walkState.supportFoot == true) {
        mSupportFoot = mRobot->getBodyNode("r_sole");
        mSwingFoot = mRobot->getBodyNode("l_sole");
    } else {
        mSupportFoot = mRobot->getBodyNode("l_sole");
        mSwingFoot = mRobot->getBodyNode("r_sole");
    }

    Eigen::MatrixXd Jacobian_supportToBaseDeriv;

    if (balancePoint == TORSO) {
	Jacobian_supportToBaseDeriv =  mRobot->getJacobianSpatialDeriv(mTorso,mSupportFoot) - mRobot->getJacobianSpatialDeriv(mSupportFoot,mSupportFoot); //(mBase,mSupportFoot)
    } else {
	Jacobian_supportToBaseDeriv =  (mRobot->getJacobianSpatialDeriv(mSupportFoot) - mRobot->getJacobianSpatialDeriv(mSupportFoot,mSupportFoot));
    }

    for (unsigned int i=0; i<44; i++){ 
    if (i!=5 || i!=6)  Jacobian_supportToBaseDeriv.col(i).setZero();
    }

    Eigen::MatrixXd Jacobian_SupportToSwingDeriv =  mRobot->getJacobianClassicDeriv(mSwingFoot,mSupportFoot) - mRobot->getJacobianClassicDeriv(mSupportFoot,mSupportFoot);

    Eigen::MatrixXd Jacobian_tot_(12, 56);

    Jacobian_tot_ << Jacobian_supportToBaseDeriv, Jacobian_SupportToSwingDeriv;

    Eigen::MatrixXd Jacobian_tot(12, 50);

    // Remove the floating base columns
    Jacobian_tot = Jacobian_tot_.block<12,50>(0, 6);

    return Jacobian_tot;
/**/ 

}

Eigen::VectorXd Controller::getJointVelocitiesQp(State current, State desired) {
    int nVariables = 50;

    double jointVelocitiesGain = 0.00001;//0.00001;
    Eigen::MatrixXd taskGain = Eigen::MatrixXd::Identity(12,12);

    // Torso Orientation
    taskGain(0,0) = 1;
    taskGain(1,1) = 1;
    taskGain(2,2) = 1;

    // CoM Position
    taskGain(3,3) = 5;
    taskGain(4,4) = 5;
    taskGain(5,5) = 1;


    // Swing Foot Orientation
    taskGain(6,6) = 5;
    taskGain(7,7) = 5;
    taskGain(8,8) = 1;

    // Swing Foot Position
    taskGain(9,9) = 5;
    taskGain(10,10) = 5;
    taskGain(11,11) = 5;

    // Construct stack of error, wrap to Pi angular errors
    Eigen::VectorXd desired_pos(12);
    //desired_pos << desired.getRelComPose(walkState.supportFoot), desired.getRelSwingFootPose(walkState.supportFoot);

    Eigen::VectorXd actual_pos(12);
    //actual_pos << current.getRelComPose(walkState.supportFoot), current.getRelSwingFootPose(walkState.supportFoot);


    Eigen::VectorXd errorStack = actual_pos - desired_pos;
    errorStack(0) = wrapToPi(errorStack(0));
    errorStack(1) = wrapToPi(errorStack(1));
    errorStack(2) = wrapToPi(errorStack(2));
    errorStack(6) = wrapToPi(errorStack(6));
    errorStack(7) = wrapToPi(errorStack(7));
    errorStack(8) = wrapToPi(errorStack(8));

    // Cost Function
    Eigen::MatrixXd Jacobian_tot = getTorsoAndSwfJacobian();
    Eigen::MatrixXd costFunctionH = mWorld->getTimeStep()*mWorld->getTimeStep()*Jacobian_tot.transpose()*taskGain*Jacobian_tot +
			jointVelocitiesGain*Eigen::MatrixXd::Identity(nVariables,nVariables);

    Eigen::VectorXd costFunctionF = mWorld->getTimeStep()*Jacobian_tot.transpose() * taskGain * IKerrorGain * errorStack;

    // Constraint RHipYawPitch and LHipYawPitch to be at the same angle
    // not working on HRP4
    Eigen::MatrixXd AHip = Eigen::MatrixXd::Zero(1,nVariables);
    AHip(0,52-6) = -mWorld->getTimeStep();
    AHip(0,46-6) = mWorld->getTimeStep();
    Eigen::VectorXd bHip(1);
    bHip(0) = mRobot->getPosition(54-6) - mRobot->getPosition(48-6);

    // Solve the QP
    Eigen::VectorXd solution = solveQP(costFunctionH, costFunctionF, 0*AHip, 0*bHip, 0*bHip);

    return solution;
}

Eigen::VectorXd Controller::getJointVelocitiesStacked(Eigen::VectorXd desider_pos_base, Eigen::VectorXd actPosBase,
		Eigen::VectorXd desider_pos_SwingFoot_front, Eigen::VectorXd actPosSwingFoot_front, Eigen::VectorXd desider_pos_SwingFoot_back, Eigen::VectorXd actPosSwingFoot_back, Eigen::VectorXd desider_com_vel){

        Eigen::VectorXd ComVref = Eigen::VectorXd::Zero(24);
        ComVref(3) = desider_com_vel(0);
        ComVref(4) = desider_com_vel(1);
        ComVref(5) = desider_com_vel(2);
     
        // in order: virtual_support_foot, opposite_to_front_foot, SwingFoot_front, SwingFoot_back;

        // virtual_support_foot, opposite_to_front_foot  equally contribute 

	Eigen::VectorXd desired_pos(24);
	desired_pos << -desider_pos_base, -desider_pos_base, desider_pos_SwingFoot_front-desider_pos_base, desider_pos_SwingFoot_back-desider_pos_base;

	// Assemble actual positions and orientations
	Eigen::VectorXd actual_pos(24);
	actual_pos << -actPosBase, -actPosBase, actPosSwingFoot_front-actPosBase, actPosSwingFoot_back-actPosBase;

        // running 0.18 foostep x
        //desired_pos(1) = desired_pos(1) + 0.01;
        //desired_pos(3) = desired_pos(3) + 0.04; //0.05
        // running 0.2 foostep x


	// Get the proper jacobian and pseudoinvert it
	Eigen::MatrixXd Jacobian_tot = getTorsoAndSwfJacobian();
	Eigen::MatrixXd PseudoJacobian_tot = (Jacobian_tot.transpose())*(Jacobian_tot*Jacobian_tot.transpose()+0.001*Eigen::MatrixXd::Identity(24,24)).inverse();

	Eigen::MatrixXd _taskGain = Eigen::MatrixXd::Identity(24,24);




        //// VERY GOOD
	// Torso Orientation
	_taskGain(0,0) = 1;//0.1
	_taskGain(1,1) = 1;//0.1;
	_taskGain(2,2) = 1;//0.1;

	// CoM Position
	_taskGain(3,3) = 50;  //0.1  10
	_taskGain(4,4) = 50;
	_taskGain(5,5) = 50; //1 IT WAS ONE BEFORE !!!!

	// Front Swing Foot Orientation
	_taskGain(6,6) = 1;
	_taskGain(7,7) = 1;
	_taskGain(8,8) = 1;

	// Front Swing Foot Position
	_taskGain(9,9) = 50;
	_taskGain(10,10) = 50;
	_taskGain(11,11) = 50;

	// Back Swing Foot Orientation
	_taskGain(12,12) = 2;
	_taskGain(13,13) = 6;
	_taskGain(14,14) = 1;

	// Back Swing Foot Position
	_taskGain(15,15) = 5;
	_taskGain(16,16) = 5;
	_taskGain(17,17) = 5;

	// Back Swing Foot Orientation
	_taskGain(18,18) = 2;
	_taskGain(19,19) = 6;
	_taskGain(20,20) = 1;

	// Back Swing Foot Position
	_taskGain(21,21) = 5;
	_taskGain(22,22) = 5;
	_taskGain(23,23) = 5;


        double ikGain = 7;

	Eigen::VectorXd qDot(12);
	qDot = PseudoJacobian_tot*(ComVref+ikGain*_taskGain*(desired_pos - actual_pos));

	return qDot;
}


Eigen::VectorXd Controller::getJointVelocitiesStacked_frontRight(Eigen::VectorXd desider_pos_SwingFoot, Eigen::VectorXd actPosSwingFoot, Eigen::VectorXd desider_com_vel){

        Eigen::VectorXd ComVref = Eigen::VectorXd::Zero(6);
        ComVref<<0.0,0.0,0.0,desider_com_vel(0), desider_com_vel(1), desider_com_vel(5);
 
	Eigen::VectorXd desired_pos(6);
	desired_pos << desider_pos_SwingFoot;

	// Assemble actual positions and orientations
	Eigen::VectorXd actual_pos(6);
	actual_pos << actPosSwingFoot;

	// Get the proper jacobian and pseudoinvert it
	Eigen::MatrixXd Jacobian_tot = getTorsoAndSwfJacobian_frontRight();
	Eigen::MatrixXd PseudoJacobian_tot = (Jacobian_tot.transpose())*(Jacobian_tot*Jacobian_tot.transpose()+0.001*Eigen::MatrixXd::Identity(6,6)).inverse();

	Eigen::MatrixXd _taskGain = Eigen::MatrixXd::Identity(6,6);

        //// VERY GOOD
	// Torso Orientation
	_taskGain(0,0) = 1;
	_taskGain(1,1) = 1;//0.001;
	_taskGain(2,2) = 1;//0;

	// CoM Position
	_taskGain(3,3) = 5;  //0.1  10
	_taskGain(4,4) = 5;
	_taskGain(5,5) = 5; //1 IT WAS ONE BEFORE !!!!

        double ikGain = 10;
	Eigen::VectorXd qDot(12); 
	qDot = PseudoJacobian_tot*(ComVref+ikGain*_taskGain*(desired_pos - actual_pos));

	return qDot;
}


Eigen::VectorXd Controller::getJointVelocitiesStacked_frontLeft(Eigen::VectorXd desider_pos_SwingFoot, Eigen::VectorXd actPosSwingFoot, Eigen::VectorXd desider_com_vel){

        Eigen::VectorXd ComVref = Eigen::VectorXd::Zero(6);
        ComVref<<0.0,0.0,0.0,desider_com_vel(0), desider_com_vel(1), desider_com_vel(5);
 
	Eigen::VectorXd desired_pos(6);
	desired_pos << desider_pos_SwingFoot;

	// Assemble actual positions and orientations
	Eigen::VectorXd actual_pos(6);
	actual_pos << actPosSwingFoot;

	// Get the proper jacobian and pseudoinvert it
	Eigen::MatrixXd Jacobian_tot = getTorsoAndSwfJacobian_frontLeft();
	Eigen::MatrixXd PseudoJacobian_tot = (Jacobian_tot.transpose())*(Jacobian_tot*Jacobian_tot.transpose()+0.001*Eigen::MatrixXd::Identity(6,6)).inverse();

	Eigen::MatrixXd _taskGain = Eigen::MatrixXd::Identity(6,6);

        //// VERY GOOD
	// Torso Orientation
	_taskGain(0,0) = 1;
	_taskGain(1,1) = 1;//0.001;
	_taskGain(2,2) = 1;//0;

	// CoM Position
	_taskGain(3,3) = 5;  //0.1  10
	_taskGain(4,4) = 5;
	_taskGain(5,5) = 5; //1 IT WAS ONE BEFORE !!!!

        double ikGain = 10;
	Eigen::VectorXd qDot(12); 
	qDot = PseudoJacobian_tot*(ComVref+ikGain*_taskGain*(desired_pos - actual_pos));

	return qDot;
}

Eigen::VectorXd Controller::getJointVelocitiesStacked_backRight(Eigen::VectorXd desider_pos_SwingFoot, Eigen::VectorXd actPosSwingFoot, Eigen::VectorXd desider_com_vel){

        Eigen::VectorXd ComVref = Eigen::VectorXd::Zero(6);
        ComVref<<0.0,0.0,0.0,desider_com_vel(0), desider_com_vel(1), desider_com_vel(5);
 
	Eigen::VectorXd desired_pos(6);
	desired_pos << desider_pos_SwingFoot;

	// Assemble actual positions and orientations
	Eigen::VectorXd actual_pos(6);
	actual_pos << actPosSwingFoot;

	// Get the proper jacobian and pseudoinvert it
	Eigen::MatrixXd Jacobian_tot = getTorsoAndSwfJacobian_backRight();
	Eigen::MatrixXd PseudoJacobian_tot = (Jacobian_tot.transpose())*(Jacobian_tot*Jacobian_tot.transpose()+0.001*Eigen::MatrixXd::Identity(6,6)).inverse();

	Eigen::MatrixXd _taskGain = Eigen::MatrixXd::Identity(6,6);

        //// VERY GOOD
	// Torso Orientation
	_taskGain(0,0) = 1;
	_taskGain(1,1) = 1;//0.001;
	_taskGain(2,2) = 1;//0;

	// CoM Position
	_taskGain(3,3) = 5;  //0.1  10
	_taskGain(4,4) = 5;
	_taskGain(5,5) = 5; //1 IT WAS ONE BEFORE !!!!

        double ikGain = 10;
	Eigen::VectorXd qDot(12); 
	qDot = PseudoJacobian_tot*(ComVref+ikGain*_taskGain*(desired_pos - actual_pos));

	return qDot;
}


Eigen::VectorXd Controller::getJointVelocitiesStacked_backLeft(Eigen::VectorXd desider_pos_SwingFoot, Eigen::VectorXd actPosSwingFoot, Eigen::VectorXd desider_com_vel){

        Eigen::VectorXd ComVref = Eigen::VectorXd::Zero(6);
        ComVref<<0.0,0.0,0.0,desider_com_vel(0), desider_com_vel(1), desider_com_vel(5);
 
	Eigen::VectorXd desired_pos(6);
	desired_pos << desider_pos_SwingFoot;

	// Assemble actual positions and orientations
	Eigen::VectorXd actual_pos(6);
	actual_pos << actPosSwingFoot;

	// Get the proper jacobian and pseudoinvert it
	Eigen::MatrixXd Jacobian_tot = getTorsoAndSwfJacobian_backLeft();
	Eigen::MatrixXd PseudoJacobian_tot = (Jacobian_tot.transpose())*(Jacobian_tot*Jacobian_tot.transpose()+0.001*Eigen::MatrixXd::Identity(6,6)).inverse();

	Eigen::MatrixXd _taskGain = Eigen::MatrixXd::Identity(6,6);

        //// VERY GOOD
	// Torso Orientation
	_taskGain(0,0) = 1;
	_taskGain(1,1) = 1;//0.001;
	_taskGain(2,2) = 1;//0;

	// CoM Position
	_taskGain(3,3) = 5;  //0.1  10
	_taskGain(4,4) = 5;
	_taskGain(5,5) = 5; //1 IT WAS ONE BEFORE !!!!

        double ikGain = 10;
	Eigen::VectorXd qDot(12); 
	qDot = PseudoJacobian_tot*(ComVref+ikGain*_taskGain*(desired_pos - actual_pos));

	return qDot;
}




























Eigen::VectorXd Controller::getJointVelocitiesQpAcceleration(State current, State desired) {
    int nVariables = 24;

    // Fetch Jacobian and Jacobian derivative
    Eigen::MatrixXd jac = getTorsoAndSwfJacobian();
    Eigen::MatrixXd jacDot = getTorsoAndSwfJacobianDeriv();

    Eigen::VectorXd q(24);
    Eigen::VectorXd qDot(24);
    Eigen::VectorXd qDoubleDot(24);

    for (int i = 0; i < 24; ++i){
        q(i) = mRobot->getPosition(i+6);
        qDot(i) = mRobot->getVelocity(i+6);
        qDoubleDot(i) = mRobot->getAcceleration(i+6);
    }

    // Construct stack of error, wrap to Pi angular errors
    Eigen::VectorXd desiredPos(12);
    //desiredPos << desired.getRelComPose(walkState.supportFoot), desired.getRelSwingFootPose(walkState.supportFoot);
    Eigen::VectorXd actual_pos(12);
    //actual_pos << current.getRelComPose(walkState.supportFoot), current.getRelSwingFootPose(walkState.supportFoot);

    Eigen::VectorXd errorStack = actual_pos - desiredPos;
    errorStack(0) = wrapToPi(errorStack(0));
    errorStack(1) = wrapToPi(errorStack(1));
    errorStack(2) = wrapToPi(errorStack(2));
    errorStack(6) = wrapToPi(errorStack(6));
    errorStack(7) = wrapToPi(errorStack(7));
    errorStack(8) = wrapToPi(errorStack(8));

    // Get desired velocity and acceleration
    Eigen::VectorXd desiredVel(12);
    desiredVel << Eigen::VectorXd::Zero(6), Eigen::VectorXd::Zero(6);
    Eigen::VectorXd desiredAcc(12);
    desiredAcc << Eigen::VectorXd::Zero(6), Eigen::VectorXd::Zero(6);

    // Cost Function
    double delta = mWorld->getTimeStep();
    double alpha = 0;
    double beta = 0;
    double gamma = 100000;
    
    Eigen::MatrixXd costFunctionH = alpha * (jac + jacDot * delta).transpose() * (jac + jacDot * delta) +
                                    beta  * (jac * delta).transpose() * (jac * delta) +
                                    gamma * (jac * delta * delta).transpose() * (jac * delta * delta) +
                                    Eigen::MatrixXd::Identity(nVariables, nVariables);

    Eigen::VectorXd costFunctionF = alpha * (jac + jacDot * delta).transpose() * (jacDot * qDot - desiredAcc) +
                                    beta  * (jac * delta).transpose() * (jac * qDot - desiredVel + delta * jac * qDoubleDot) +
                                    gamma * (jac * delta * delta).transpose() * (errorStack + jac * qDot);

    // Dummy constraint
    Eigen::MatrixXd AHip = Eigen::MatrixXd::Zero(1,nVariables);
    Eigen::VectorXd bHip(1);
    bHip(0) = 0;

    // Solve the QP
    Eigen::VectorXd solution = solveQP(costFunctionH, costFunctionF, AHip, bHip, bHip);
    return solution;
}

Eigen::Vector3d Controller::getRPY(dart::dynamics::BodyNode* body, dart::dynamics::BodyNode* referenceFrame) {
    Eigen::MatrixXd rotMatrix = body->getTransform(referenceFrame).rotation();

    Eigen::Vector3d RPY;
    RPY << atan2(rotMatrix(2,1),rotMatrix(2,2)),
        atan2(-rotMatrix(2,0),sqrt(rotMatrix(2,1)*rotMatrix(2,1)+rotMatrix(2,2)*rotMatrix(2,2))),
        atan2(rotMatrix(1,0),rotMatrix(0,0));

    return RPY;
}

Eigen::Vector3d Controller::getRPY(dart::dynamics::BodyNode* body) {
    Eigen::MatrixXd rotMatrix = body->getTransform().rotation();

    Eigen::Vector3d RPY;
    RPY << atan2(rotMatrix(2,1),rotMatrix(2,2)),
        atan2(-rotMatrix(2,0),sqrt(rotMatrix(2,1)*rotMatrix(2,1)+rotMatrix(2,2)*rotMatrix(2,2))),
        atan2(rotMatrix(1,0),rotMatrix(0,0));

    return RPY;
}

Eigen::Vector3d Controller::getZmpFromExternalForces()
{
    Eigen::Vector3d zmp_v;
    bool left_contact = false;
    bool right_contact = false;
/*
    Eigen::Vector3d left_cop;
    if(abs(mLeftFoot->getConstraintImpulse()[5]) > 0.01){
        left_cop << -mLeftFoot->getConstraintImpulse()(1)/mLeftFoot->getConstraintImpulse()(5), mLeftFoot->getConstraintImpulse()(0)/mLeftFoot->getConstraintImpulse()(5), 0.0;
        Eigen::Matrix3d iRotation = mLeftFoot->getWorldTransform().rotation();
        Eigen::Vector3d iTransl   = mLeftFoot->getWorldTransform().translation();
        left_cop = iTransl + iRotation*left_cop;
        left_contact = true;
    }

    Eigen::Vector3d right_cop;
    if(abs(mRightFoot->getConstraintImpulse()[5]) > 0.01){
        right_cop << -mRightFoot->getConstraintImpulse()(1)/mRightFoot->getConstraintImpulse()(5), mRightFoot->getConstraintImpulse()(0)/mRightFoot->getConstraintImpulse()(5), 0.0;
        Eigen::Matrix3d iRotation = mRightFoot->getWorldTransform().rotation();
        Eigen::Vector3d iTransl   = mRightFoot->getWorldTransform().translation();
        right_cop = iTransl + iRotation*right_cop;
        right_contact = true;
    }

    if(left_contact && right_contact){
        zmp_v << (left_cop(0)*mLeftFoot->getConstraintImpulse()[5] + right_cop(0)*mRightFoot->getConstraintImpulse()[5])/(mLeftFoot->getConstraintImpulse()[5] + mRightFoot->getConstraintImpulse()[5]),
                 (left_cop(1)*mLeftFoot->getConstraintImpulse()[5] + right_cop(1)*mRightFoot->getConstraintImpulse()[5])/(mLeftFoot->getConstraintImpulse()[5] + mRightFoot->getConstraintImpulse()[5]),
		 0.0;
    }else if(left_contact){
        zmp_v << left_cop(0), left_cop(1), 0.0;
    }else if(right_contact){
        zmp_v << right_cop(0), right_cop(1), 0.0;
    }else{
        // No contact detected
        zmp_v << 0.0, 0.0, 0.0;
    }
/**/
    return zmp_v;
}

void Controller::keyboard(unsigned char /*_key*/, int /*_x*/, int /*_y*/) {}

void Controller::setInitialConfiguration() {
  Eigen::VectorXd q = mRobot->getPositions();


    // Floating Base
    q[0] = 0.0;
    q[1] = 0.0;
    q[2] = 0.0;
    q[3] = 0.44;
    q[4] = 0.00;
    q[5] = 0.605;

    mRobot->setPositions(q);

/*
    // Right Leg
    q[44] = 0.0;            // hip yaw
    q[45] = 3*M_PI/180;    // hip roll
    q[46] = -25*M_PI/180;   // hip pitch
    q[47] = 50*M_PI/180;   // knee pitch
    q[48] = -30*M_PI/180; // ankle pitch
    q[49] = -4*M_PI/180;   // ankle roll         
    // Left Leg
    q[50] = 0.0;            // hip yaw
    q[51] = -3*M_PI/180;    // hip roll
    q[52] = -25*M_PI/180;   // hip pitch
    q[53] = 50*M_PI/180;   // knee pitch
    q[54] = -30*M_PI/180; // ankle pitch
    q[55] = 4*M_PI/180;   // ankle roll       



    // Additional arm position setting
    mRobot->setPosition(mRobot->getDof("R_SHOULDER_P")->getIndexInSkeleton(), (4)*M_PI/180 );
    mRobot->setPosition(mRobot->getDof("R_SHOULDER_R")->getIndexInSkeleton(), -8*M_PI/180  );
    mRobot->setPosition(mRobot->getDof("R_SHOULDER_Y")->getIndexInSkeleton(), 0 );

    mRobot->setPosition(mRobot->getDof("R_ELBOW_P")->getIndexInSkeleton(), -25*M_PI/180 );
    mRobot->setPosition(mRobot->getDof("R_ELBOW_P")->getIndexInSkeleton(), -25*M_PI/180 );

    mRobot->setPosition(mRobot->getDof("L_SHOULDER_P")->getIndexInSkeleton(), (4)*M_PI/180  );
    mRobot->setPosition(mRobot->getDof("L_SHOULDER_R")->getIndexInSkeleton(), 8*M_PI/180  );
    mRobot->setPosition(mRobot->getDof("L_SHOULDER_Y")->getIndexInSkeleton(), 0 );

    mRobot->setPosition(mRobot->getDof("L_ELBOW_P")->getIndexInSkeleton(), -25*M_PI/180 ); 
    mRobot->setPosition(mRobot->getDof("L_ELBOW_P")->getIndexInSkeleton(), -25*M_PI/180 );

    // Set initial chest position 
    mRobot->setPosition(mRobot->getDof("CHEST_P")->getIndexInSkeleton(), 0.02 );
    mRobot->setPosition(mRobot->getDof("CHEST_Y")->getIndexInSkeleton(), 0.02 );
/**/ 
}


void Controller::ArmSwing() {

// TO DO: add variable period to the swing trajecotry
mRobot->setPosition(mRobot->getDof("R_SHOULDER_P")->getIndexInSkeleton(), (4-2.5*sin(2*M_PI*0.01*(mWorld->getSimFrames())/0.9))*M_PI/180 );   // 2.5 amplitude
mRobot->setPosition(mRobot->getDof("L_SHOULDER_P")->getIndexInSkeleton(), (4+2.5*sin(2*M_PI*0.01*(mWorld->getSimFrames())/0.9))*M_PI/180 );

}

void Controller::FixWaistToChestJoints() {


mRobot->setPosition(mRobot->getDof("CHEST_P")->getIndexInSkeleton(), 0.0 );
mRobot->setPosition(mRobot->getDof("CHEST_Y")->getIndexInSkeleton(), 0.0 );


}


void Controller::storeData() {

     Eigen::VectorXd COMPOS_meas  = Eigen::VectorXd::Zero(3);
     COMPOS_meas = mRobot->getCOM();
     Eigen::VectorXd COMVEL_meas  = Eigen::VectorXd::Zero(3);
     COMVEL_meas = mTorso->getCOMLinearVelocity();
     Eigen::VectorXd FOOT_meas  = Eigen::VectorXd::Zero(3);
     FOOT_meas = mSupportFoot->getCOM();

     Eigen::VectorXd ZMPPOS_meas_cop =  Eigen::VectorXd::Zero(3);
     ZMPPOS_meas_cop = COMPOS_meas - mRobot->getCOMLinearAcceleration()/(eta*eta);
     ZMPPOS_meas_cop = getZmpFromExternalForces();
     

     ofstream myfile;

     myfile.open ("./Data/x_RF.txt",ios::app);
     myfile << mSupportFoot->getCOM()(0) <<endl; 
     myfile.close();
     myfile.open ("./Data/y_RF.txt",ios::app);
     myfile << mSupportFoot->getCOM()(1) <<endl; 
     myfile.close();

     myfile.open ("./Data/x_m.txt",ios::app);
     myfile << COMPOS_meas(0) <<endl; 
     myfile.close();
     myfile.open ("./Data/y_m.txt",ios::app);
     myfile << COMPOS_meas(1) <<endl; 
     myfile.close();
     myfile.open ("./Data/xz_m_cop.txt",ios::app);
     myfile << ZMPPOS_meas_cop(0) <<endl; //current.zmpPos
     myfile.close();
     myfile.open ("./Data/yz_m_cop.txt",ios::app);
     myfile << ZMPPOS_meas_cop(1) <<endl; 
     myfile.close();
     myfile.open ("./Data/x.txt",ios::app);
     myfile << desired.comPos(0) <<endl; 
     myfile.close();
     myfile.open ("./Data/y.txt",ios::app);
     myfile << desired.comPos(1) <<endl; 
     myfile.close();
     myfile.open ("./Data/xz.txt",ios::app);
     myfile << desired.zmpPos(0) <<endl; 
     myfile.close();
     myfile.open ("./Data/yz.txt",ios::app);
     myfile << desired.zmpPos(1) <<endl; 
     myfile.close();

}


void Controller::createForceLine()
  {

    dart::dynamics::SimpleFramePtr lineFrame
        = std::make_shared<dart::dynamics::SimpleFrame>(
            dart::dynamics::Frame::World());

    mForceLine = std::make_shared<dart::dynamics::LineSegmentShape>( current.comPos, previousCoM, 3.0);
    mForceLine->addDataVariance(dart::dynamics::Shape::DYNAMIC_VERTICES);
    mForceLine->addConnection(0,1) ;		


    lineFrame->setShape(mForceLine);
    lineFrame->createVisualAspect();
    lineFrame->getVisualAspect()->setColor(DefaultForceLineColor);

    mWorld->addSimpleFrame(lineFrame);

  }
