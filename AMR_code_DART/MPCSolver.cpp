#include "MPCSolver.hpp"

static std::ofstream foutDebug(realpath("../data/debug.txt", NULL), std::ofstream::trunc);

MPCSolver::MPCSolver(const Eigen::MatrixXd& ftsp_and_timings) {

     // read .txt files to store footsteps - first store the file size
     std::ifstream inFile; inFile.open("../vertical_motion/z.txt");
     if (!inFile) {
          std::cout << "Unable to open File" << std::endl;
          exit(1);
     }
     int file_size = 0; double buffer;
     while (inFile >> buffer) {
         file_size = file_size + 1;
     }
     inFile.close();
     z_trajectory = Eigen::VectorXd::Zero(file_size); 
     f_trajectory = Eigen::VectorXd::Zero(file_size); 
     std::ifstream inFile_; inFile.open("../vertical_motion/z.txt"); inFile_.open("../vertical_motion/f.txt");
     if (!inFile) {
          std::cout << "Unable to open File" << std::endl;
          exit(1);
     }
     for (int i = 0; i < file_size; i++) {
          inFile >> z_trajectory(i);
          inFile_ >> f_trajectory(i);
     }
     inFile.close(); inFile_.close();

    // HORIZONTAL QP

    // Matrices for cost function (unique QP for XY)
    costFunctionH = Eigen::MatrixXd::Zero(2*(N+M),2*(N+M));
    costFunctionF = Eigen::VectorXd::Zero(2*(N+M));

    // Matrices for stability constraint
    Aeq = Eigen::MatrixXd::Zero(2,(N*2)+(M*2));
    beq = Eigen::VectorXd::Zero(2);

    // Matrices for ZMP constraint
    AZmp = Eigen::MatrixXd::Zero(2*N,2*(N+M));
    bZmpMax = Eigen::VectorXd::Zero(2*N);
    bZmpMin = Eigen::VectorXd::Zero(2*N);

    // Matrices for kinematic constraint
    AFootsteps = Eigen::MatrixXd::Zero(2*M,2*(M+N));
    bFootstepsMax = Eigen::VectorXd::Zero(2*M);
    bFootstepsMin = Eigen::VectorXd::Zero(2*M);

    // Matrices for all constraints stacked
    AConstraint = Eigen::MatrixXd::Zero(2*(N+M)+2,2*(N+M));
    bConstraintMin = Eigen::VectorXd::Zero(2*(N+M)+2);
    bConstraintMax = Eigen::VectorXd::Zero(2*(N+M)+2);

    // Matrices for ZMP prediction (might be useful)
    p = Eigen::VectorXd::Ones(N);
    P = Eigen::MatrixXd::Ones(N,N)*mpcTimeStep;
    for(int i = 0; i < N; ++i) {
    	for(int j = 0; j < N; ++j) {
            if (j > i) P(i, j)=0;
        }
    }

    // TO DO
    double ch = cosh(eta*mpcTimeStep);
    double sh = sinh(eta*mpcTimeStep);
    Eigen::MatrixXd A_upd = Eigen::MatrixXd::Zero(3,3);
    Eigen::VectorXd B_upd = Eigen::VectorXd::Zero(3);
    A_upd<<ch,sh/eta,1-ch,eta*sh,ch,-eta*sh,0,0,1;
    B_upd<<mpcTimeStep-sh/eta,1-ch,mpcTimeStep;

    // Fs sequence midpoint for anticipative tail
    // TO DO:
/*
    x_midpoint.resize(500);
    y_midpoint.resize(500);
    int midpointFootstepIndex = 0;
    for (int i = 0; i < 500; i++) {
        int timeTillFootstepEnd = footstepPlan->getFootstepEndTiming(midpointFootstepIndex) - i;

        if (timeTillFootstepEnd <= 0) { // is footstep finished?
            midpointFootstepIndex++;
            timeTillFootstepEnd = footstepPlan->getFootstepEndTiming(midpointFootstepIndex) - i;
        }

        if (timeTillFootstepEnd > ds_samples) { // are we in single support?
            x_midpoint(i) = footstepPlan->getFootstepPosition(midpointFootstepIndex)(0);
            y_midpoint(i) = footstepPlan->getFootstepPosition(midpointFootstepIndex)(1);
        } else { // we are in double support
            x_midpoint(i) = footstepPlan->getFootstepPosition(midpointFootstepIndex)(0) * (double)timeTillFootstepEnd / ds_samples
                            + footstepPlan->getFootstepPosition(midpointFootstepIndex + 1)(0) * (1 - (double)timeTillFootstepEnd / ds_samples);
            y_midpoint(i) = footstepPlan->getFootstepPosition(midpointFootstepIndex)(1) * (double)timeTillFootstepEnd / ds_samples
                            + footstepPlan->getFootstepPosition(midpointFootstepIndex + 1)(1) * (1 - (double)timeTillFootstepEnd / ds_samples);
        }
    }
/**/
     old_fsCount = 0;
     ct = 0;

     xz_dot = 0.0;
     yz_dot = 0.0;

    // VERTICAL QP


    // Matrices for cost function (vertical QP)
    costFunctionH_z = Eigen::MatrixXd::Zero(N,N);
    costFunctionF_z = Eigen::VectorXd::Zero(N);

    // Matrices for equality constraint (vertical QP) - NOT USED IN THE WALKING GAIT


    // Matrices for vertical CoM constraints
    Aineq_z = Eigen::MatrixXd::Zero(N,N);
    bMax_z = Eigen::VectorXd::Zero(N);
    bMin_z = Eigen::VectorXd::Zero(N);
    // Normal force bounds 0<f<inf
    Af_z = Eigen::MatrixXd::Zero(N,N);
    fMax_z = Eigen::VectorXd::Zero(N);
    fMin_z = Eigen::VectorXd::Zero(N);

    // Build here the constraints for z
    A_z = Eigen::MatrixXd::Zero(2,2);
    A_z << 1.0, mpcTimeStep, 0.0, 1.0;
    B_z = Eigen::MatrixXd::Zero(2,1); B_z << 0.0,  mpcTimeStep/ mass_hrp4;
    Bg_z = Eigen::MatrixXd::Zero(2,1); Bg_z << 0.0, - mpcTimeStep;
    C_z = Eigen::MatrixXd::Zero(1,2); C_z << 1.0, 0.0;
    C_p = Eigen::MatrixXd::Zero(1,2); C_p << 1.0, 0.0;
    C_v = Eigen::MatrixXd::Zero(1,2); C_v << 0.0, 1.0;    
    S_bar_z = Eigen::MatrixXd::Zero(N,N);
    T_bar_z = Eigen::MatrixXd::Zero(N,2);    
    S_bar_g_z = Eigen::MatrixXd::Zero(N,N);
    T_bar_g_z = Eigen::MatrixXd::Zero(N,2);
    S_bar_z_v = Eigen::MatrixXd::Zero(N,N);
    T_bar_z_v = Eigen::MatrixXd::Zero(N,2);    
    S_bar_g_z_v = Eigen::MatrixXd::Zero(N,N);
    T_bar_g_z_v = Eigen::MatrixXd::Zero(N,2);
    S_p = Eigen::MatrixXd::Zero(N,N);
    S_v = Eigen::MatrixXd::Zero(N,N);
    T_p = Eigen::MatrixXd::Zero(N,2);
    T_v = Eigen::MatrixXd::Zero(N,2);

    Eigen::MatrixXd& A_z_ = A_z;
    for (int k = 0; k < N; k++) {
        T_bar_z.block(k,0,1,2) = C_z*(matrixPower(A_z_,k+1));
        T_bar_z_v.block(k,0,1,2) = C_v*(matrixPower(A_z_,k+1));
        for (int j = 0; j < k; j++) {
             S_bar_z.block(k,j,1,1) = C_z*(matrixPower(A_z_,k-j))*B_z;
             S_bar_z_v.block(k,j,1,1) = C_v*(matrixPower(A_z_,k-j))*B_z;
             S_bar_g_z.block(k,j,1,1) = C_z*(matrixPower(A_z_,k-j))*Bg_z;
             S_bar_g_z_v.block(k,j,1,1) = C_v*(matrixPower(A_z_,k-j))*Bg_z;
        }
    }
    T_bar_g_z = S_bar_g_z*p*g;
    T_bar_g_z_v = S_bar_g_z_v*p*g;
    costFunctionH_z = 0.01*S_bar_z_v.transpose()*S_bar_z_v + (Eigen::MatrixXd::Identity(N,N)) ;
    Aineq_z << S_bar_z;
    bMax_z.block(0,0,N,1) = p*10000.0;
    bMin_z.block(0,0,N,1) = p*0.0;
    // initialize lambda vector
    lambda = p*eta*eta;
    etas = p*eta;


    // build midpoint sequence
    ftsp_midpoint = Eigen::MatrixXd::Zero(ftsp_and_timings.rows()*(S+F),3);
    Eigen::MatrixXd double_support_transition = Eigen::MatrixXd::Zero(F,1);
    Eigen::MatrixXd p_dst = Eigen::MatrixXd::Ones(F,1);
    for (int i = 0; i < F; i++) double_support_transition(i,0) = (double)i/(double)F;
    for (int i = 0; i < ftsp_and_timings.rows() - 1 ; i++) {
        // single support phases
        ftsp_midpoint.block(i*(S+F),0,S,1) = ftsp_and_timings(i,0)*Eigen::MatrixXd::Ones(S,1);
        ftsp_midpoint.block(i*(S+F),1,S,1) = ftsp_and_timings(i,1)*Eigen::MatrixXd::Ones(S,1);
        ftsp_midpoint.block(i*(S+F),2,S,1) = ftsp_and_timings(i,2)*Eigen::MatrixXd::Ones(S,1);
        // double support phases
        ftsp_midpoint.block(i*(S+F)+S,0,F,1) = ftsp_and_timings(i,0)*p_dst + (ftsp_and_timings(i+1,0)-ftsp_and_timings(i,0))*double_support_transition;
        ftsp_midpoint.block(i*(S+F)+S,1,F,1) =  ftsp_and_timings(i,1)*p_dst + (ftsp_and_timings(i+1,1)-ftsp_and_timings(i,1))*double_support_transition;
        ftsp_midpoint.block(i*(S+F)+S,2,F,1) =  ftsp_and_timings(i,2)*p_dst + (ftsp_and_timings(i+1,2)-ftsp_and_timings(i,2))*double_support_transition;
    }
    
    // useful multiplier
    deltas = Eigen::MatrixXd::Zero(1,N);
    for (int i = 0; i<N; i++) deltas(0,i) = exp(-mpcTimeStep*eta*i);

    // cost function xy
    costFunctionH_xy = (Eigen::MatrixXd::Identity(N,N));

    // matrices and parameters fo horizontal QP
    A_xy = Eigen::MatrixXd::Zero(2,2);
    B_xy = Eigen::MatrixXd::Zero(2,1); 
    A_virtual_ZMP = Eigen::MatrixXd::Zero(F,N);
    b_virtual_ZMP_x = Eigen::VectorXd::Zero(F);
    b_virtual_ZMP_y = Eigen::VectorXd::Zero(F);
    a_f = (mpcTimeStep + (1/eta)*(1-exp(eta*mpcTimeStep)))/(1 - exp(eta*mpcTimeStep));
    C_sc = Eigen::MatrixXd::Zero(1,2);
    phi_state = Eigen::MatrixXd::Identity(2,2);
    phi_input = Eigen::MatrixXd::Zero(2,N);

}

MPCSolver::~MPCSolver() {}

State MPCSolver::solve(State current, WalkState walkState, const Eigen::MatrixXd& ftsp_and_timings) {

    itr = walkState.mpcIter;
    fsCount = walkState.footstepCounter;

    // Output struct
    State next = current;

    auto start = std::chrono::high_resolution_clock::now();

    if ((walkState.controlIter % (int)(100*mpcTimeStep)) == 0) {


    ///////// STAGE ONE: VERTICAL QP ///////// 


    Eigen::MatrixXd I_c = Eigen::MatrixXd::Identity(N,N);
    Eigen::MatrixXd Aeq_z; Eigen::VectorXd beq_z;

    if (walkState.mpcIter < S){
       Aeq_z = Eigen::MatrixXd::Zero(F, N);
       beq_z = Eigen::VectorXd::Zero(F);
    } else {
       Aeq_z = Eigen::MatrixXd::Zero(S+F-walkState.mpcIter, N);
       beq_z = Eigen::VectorXd::Zero(S+F-walkState.mpcIter);
    }

    for (int i = 0; i < N ; i++) {
        if (walkState.mpcIter < S){
           if (i >= S && i < S+F) { 
              I_c(i,i) = 0.0;
              Aeq_z(i-S,i-walkState.mpcIter) = 1.0;
           }
        } else {
           if (i < S+F-walkState.mpcIter){
              I_c(i,i) = 0.0;
              Aeq_z(i,i) = 1.0;
           }
        }
    }
/*
    std::cout << " " << std::endl;
    std::cout << Aeq_z << std::endl;
    std::cout << " " << std::endl;
/**/

    Eigen::MatrixXd state_z = Eigen::MatrixXd::Zero(2,1);
    state_z << current.comPos(2), current.comVel(2);
    double q_p, q_u, q_v;
    q_p = 1005000.0;
    q_u  = 0.01;
    q_v = 100.0;    
    

    costFunctionH_z = q_p*S_bar_z.transpose()*S_bar_z + q_v*S_bar_z_v.transpose()*S_bar_z_v + q_u*(Eigen::MatrixXd::Identity(N,N)) ;
    costFunctionF_z << q_p*S_bar_z.transpose()*(T_bar_z*state_z + T_bar_g_z - p*h_des - ftsp_midpoint.block((int)(walkState.simulationTime/(mpcTimeStep/controlTimeStep)),2,N,1)) + q_v*S_bar_z_v.transpose()*(T_bar_z_v*state_z + T_bar_g_z_v) + q_u*(-p*mass_hrp4*g);


    double is_running = 0.0;
    if (walkState.footstepCounter>1) is_running = 1.0;
    int iter_step = walkState.mpcIter;



     // Eigen::VectorXd decisionVariables_z = solveQP_hpipm_z(costFunctionH_z, costFunctionF_z, test*Aineq_z.block(0,0,N,N), test*bMin_z.block(0,0,N,1), test*bMax_z.block(0,0,N,1)); 
Eigen::VectorXd decisionVariables_z = solveQP_hpipm_z(costFunctionH_z, costFunctionF_z, Aineq_z.block(0,0,N,N), bMin_z.block(0,0,N,1),bMax_z.block(0,0,N,1), is_running*Aeq_z, is_running*beq_z); 
   

    //decisionVariables_z = Eigen::VectorXd::Ones(N)*mass_hrp4*g;

    Eigen::Vector2d nextStateZ = A_z*state_z + B_z*decisionVariables_z(0) + Bg_z*g;
    next.comPos(2) = nextStateZ(0);
    next.comVel(2) = nextStateZ(1);
    if (isnan(next.comPos(2))) next.comPos(2) = h_des;
    if (isnan(next.comVel(2))) next.comVel(2) = 0.0;

/**/
    std::cout << "mpcIter " << walkState.mpcIter << std::endl;
    std::cout << "CoM z " << next.comPos(2) << std::endl;
    std::cout << "f " << decisionVariables_z(0) << std::endl;


    // Reset constraint matrices
    AZmp.setZero();
    AFootsteps.setZero();

    ///////// STAGE TWO: COMPUTE PIECE-WISE CONSTANT LAMBDA ///////// 
    //Eigen::VectorXd z_acceleration = (1.0/mass_hrp4)*decisionVariables_z - p*g;
    //Eigen::VectorXd z_position = S_bar_z*decisionVariables_z+T_bar_z*state_z+T_bar_g_z;

    // LOADED TRAJECTORY: STAGE ONE IS PERFORMED OFF-LINE!

    Eigen::VectorXd z_acceleration = (1.0/mass_hrp4)*decisionVariables_z - p*g;
    Eigen::VectorXd z_position = S_bar_z*decisionVariables_z+T_bar_z*state_z+T_bar_g_z;

    //Eigen::VectorXd z_acceleration = (1.0/mass_hrp4)*f_trajectory.block((int)(walkState.simulationTime/(mpcTimeStep/controlTimeStep)),0,N,1) - p*g;
    //Eigen::VectorXd z_position = z_trajectory.block((int)(walkState.simulationTime/(mpcTimeStep/controlTimeStep)),0,N,1);

    //next.comPos(2) = z_trajectory((int)(walkState.simulationTime/(mpcTimeStep/controlTimeStep)));
    //next.comVel(2) = (z_trajectory((int)(walkState.simulationTime/(mpcTimeStep/controlTimeStep))+1)-z_trajectory((int)(walkState.simulationTime/(mpcTimeStep/controlTimeStep))))/controlTimeStep;

    for (int j = 0; j< N; j++){ 
        lambda(j) = (g + z_acceleration(j))/z_position(j);
        etas(j) = sqrt(lambda(j));
        if (lambda(j) <= 0.0) etas(j) = eta;
    }
    std::cout << "lambda "<< lambda(0) << std::endl;

    //////// STAGE THREE: HORIZONTAL QP ////////

    Eigen::MatrixXd state_x = Eigen::MatrixXd::Zero(2,1);
    state_x << current.comPos(0), current.comVel(0);
    Eigen::MatrixXd state_y = Eigen::MatrixXd::Zero(2,1);
    state_y << current.comPos(1), current.comVel(1);

    Eigen::VectorXd decisionVariables_x = Eigen::VectorXd::Zero(N);
    Eigen::VectorXd decisionVariables_y = Eigen::VectorXd::Zero(N);

    if (lambda(0) > 2.0) {  

    // ZMP constraints
    Eigen::VectorXd Zmin_x, Zmax_x, Zmin_y, Zmax_y;
    Eigen::MatrixXd A_ineq_xy = Eigen::MatrixXd::Identity(N,N);
 
    if (walkState.footstepCounter > 1) {
    Zmin_x = ftsp_midpoint.block((int)(walkState.simulationTime/(mpcTimeStep/controlTimeStep)),0,N,1) - p*footConstraintSquareWidth/2;
    Zmax_x = ftsp_midpoint.block((int)(walkState.simulationTime/(mpcTimeStep/controlTimeStep)),0,N,1) + p*footConstraintSquareWidth/2;
    Zmin_y = ftsp_midpoint.block((int)(walkState.simulationTime/(mpcTimeStep/controlTimeStep)),1,N,1) - p*footConstraintSquareWidth/2;
    Zmax_y = ftsp_midpoint.block((int)(walkState.simulationTime/(mpcTimeStep/controlTimeStep)),1,N,1) + p*footConstraintSquareWidth/2;
    } else {
    Zmin_x = ftsp_midpoint.block((int)(walkState.simulationTime/(mpcTimeStep/controlTimeStep)),0,N,1) - p*1;
    Zmax_x = ftsp_midpoint.block((int)(walkState.simulationTime/(mpcTimeStep/controlTimeStep)),0,N,1) + p*1;
    Zmin_y = ftsp_midpoint.block((int)(walkState.simulationTime/(mpcTimeStep/controlTimeStep)),1,N,1) - p*1;
    Zmax_y = ftsp_midpoint.block((int)(walkState.simulationTime/(mpcTimeStep/controlTimeStep)),1,N,1) + p*1;   
    }


    // Stability constraints 
    Eigen::VectorXd beq_x, beq_y; 
    Eigen::MatrixXd Aeq_x, Aeq_y; 
    Aeq_x = Eigen::MatrixXd::Zero(1,N);
    Aeq_y = Eigen::MatrixXd::Zero(1,N);
    beq_x = Eigen::VectorXd::Zero(1,1);
    beq_y = Eigen::VectorXd::Zero(1,1);  

    phi_state = Eigen::MatrixXd::Identity(2,2);
    
    for (int i = 0; i < N; i++) {

        if (lambda(i) < 2.0) {
            A_xy << 1.0, mpcTimeStep, 0.0, 1.0;
            B_xy << 0.0, 0.0;
        } else {
          double ch = cosh(sqrt(lambda(i))*mpcTimeStep);
          double sh = sinh(sqrt(lambda(i))*mpcTimeStep);
          A_xy << ch, sh/sqrt(lambda(i)), sqrt(lambda(i))*sh, ch;
          B_xy << 1-ch, -sqrt(lambda(i))*sh;
        }
        phi_state = A_xy*phi_state;
        phi_input.block(0,i,2,1) = B_xy;

        for (int j = i+1; j < N; j++) {
          double ch = cosh(sqrt(lambda(j))*mpcTimeStep);
          double sh = sinh(sqrt(lambda(j))*mpcTimeStep);
          A_xy << ch, sh/sqrt(lambda(j)), sqrt(lambda(j))*sh, ch;
          if (lambda(j) < 2.0) A_xy << 1, mpcTimeStep, 0, 1;
          phi_input.block(0,i,2,1) = A_xy*phi_input.block(0,i,2,1);
        }
        
     }

    double eta_sc = sqrt(lambda(N-1));
    eta_sc = eta;
    C_sc << 1.0, 1.0/eta_sc;
    Aeq_x = C_sc*phi_input;
    Aeq_y = C_sc*phi_input;

    Eigen::VectorXd buffer_eq_constr = -C_sc*phi_state*state_x + eta_sc*mpcTimeStep*deltas*ftsp_midpoint.block((int)(walkState.simulationTime/(mpcTimeStep/controlTimeStep))+N,0,N,1);
    beq_x(0) = buffer_eq_constr(0);
    buffer_eq_constr = -C_sc*phi_state*state_y + eta_sc*mpcTimeStep*deltas*ftsp_midpoint.block((int)(walkState.simulationTime/(mpcTimeStep/controlTimeStep))+N,1,N,1);
    beq_y(0) = buffer_eq_constr(0);


    // cost function xy
    costFunctionF_x = - ftsp_midpoint.block((int)(walkState.simulationTime/(mpcTimeStep/controlTimeStep)),0,N,1);
    costFunctionF_y = - ftsp_midpoint.block((int)(walkState.simulationTime/(mpcTimeStep/controlTimeStep)),1,N,1);



    // solve QP

    decisionVariables_x = solveQP_hpipm_xy_piecewiseconstantZMP(costFunctionH_xy, costFunctionF_x, Aeq_x, beq_x, A_ineq_xy, Zmin_x, Zmax_x);
    decisionVariables_y = solveQP_hpipm_xy_piecewiseconstantZMP(costFunctionH_xy, costFunctionF_y, Aeq_y, beq_y, A_ineq_xy, Zmin_y, Zmax_y);

    }

    //double zmp_x_input = 0.0; //decisionVariables_x(0);
    //double zmp_y_input = 0.0;//decisionVariables_y(0);
    double zmp_x_input = decisionVariables_x(0);
    double zmp_y_input = decisionVariables_y(0);

    // State integration
    if (lambda(0) < 2.0) {
       A_xy << 1.0, mpcTimeStep, 0.0, 1.0;
       B_xy << 0.0, 0.0;
    } else {
      double ch = cosh(sqrt(lambda(0))*mpcTimeStep);
      double sh = sinh(sqrt(lambda(0))*mpcTimeStep);
      A_xy << ch, sh/sqrt(lambda(0)), sqrt(lambda(0))*sh, ch;
      B_xy << 1.0-ch, -sqrt(lambda(0))*sh;
    }
    

    state_x = A_xy*state_x + B_xy*zmp_x_input;
    state_y = A_xy*state_y + B_xy*zmp_y_input;
    next.comPos(0) = state_x(0);
    next.comVel(0) = state_x(1);
    next.comPos(1) = state_y(0);
    next.comVel(1) = state_y(1);


    std::cout << "zmp x " << zmp_x_input << std::endl;
    std::cout << "next.comPos(0) " << next.comPos(0) << std::endl;
    std::cout << "next.comVel(0) " << next.comVel(0) << std::endl;


    }



    auto finish = std::chrono::high_resolution_clock::now();
    auto interval = std::chrono::time_point_cast<std::chrono::microseconds>(finish) - std::chrono::time_point_cast<std::chrono::microseconds>(start);
    //std::cout << "Elapsed time: " << interval.count() << std::endl;


    // Swing foot trajectory generator
    Eigen::Vector4d footstepPredicted;
    footstepPredicted << ftsp_and_timings(fsCount+0,0), ftsp_and_timings(fsCount+0,1), 0.0, 0.0;


    // Update swing foot position
    
    // If it's the first step, or we are in double support, keep the foot on the ground
    double timeSinceFootstepStart = (F + walkState.controlIter) * controlTimeStep;
    double singleSupportEnd = ((double)(S+F+F)*(mpcTimeStep/controlTimeStep)) * controlTimeStep;
    double swingFootHeight = -(4*stepHeight/pow(singleSupportEnd,2)) * timeSinceFootstepStart * (timeSinceFootstepStart - singleSupportEnd);
    int samplesTillNextFootstep = S+F - walkState.controlIter;

    if (walkState.footstepCounter <= 1) swingFootHeight = 0.0;

    total_flight_phase = 0.03;

/*

    // If support foot is left, move right foot
    if (walkState.supportFoot == 1) {


	next.rightFootPos += kSwingFoot * (footstepPredicted.head(3) - current.rightFootPos);
        

        next.rightFootPos(2) = swingFootHeight;
        if (samplesTillNextFootstep < F && walkState.footstepCounter > 1) next.leftFootPos(2) = total_flight_phase;
        else next.leftFootPos(2) = 0.0;
        next.rightFootVel = Eigen::Vector3d::Zero();
        next.rightFootAcc = Eigen::Vector3d::Zero();
        next.rightFootOrient(2) += 5 * timeSinceFootstepStart * kSwingFoot * wrapToPi(footstepPredicted(3) - current.rightFootOrient(2));
        //if (samplesTillNextFootstep > F + 5 && samplesTillNextFootstep < S + F - 5) next.rightFootOrient(1) = 0.03;


    } else {


	next.leftFootPos += kSwingFoot * (footstepPredicted.head(3) - current.leftFootPos);
       

        next.leftFootPos(2) = swingFootHeight;
        if (samplesTillNextFootstep < F && walkState.footstepCounter > 1) next.rightFootPos(2) = total_flight_phase;
        else next.rightFootPos(2) = 0.0;
        next.leftFootVel = Eigen::Vector3d::Zero();
        next.leftFootAcc = Eigen::Vector3d::Zero();
        next.leftFootOrient(2) += 5 * timeSinceFootstepStart * kSwingFoot * wrapToPi(footstepPredicted(3) - current.leftFootOrient(2));

        //if (samplesTillNextFootstep > F + 5 && samplesTillNextFootstep < S + F - 5) next.leftFootOrient(1) = 0.03;

    }

    //std::cout << "swingFootHeight " << swingFootHeight << std::endl;

    next.torsoOrient(2) = wrapToPi((next.leftFootOrient(2) + next.rightFootOrient(2)) / 2.0);
/**/

     
     
    // Return next robot state

    return next;
}



