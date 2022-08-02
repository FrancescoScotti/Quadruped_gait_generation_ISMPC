% Compute the gait generation with ZMP and CoM trajectories
% using ISMPC.

init_quadruped2;

% Proposed foot steps, from "init_quadruped.m"
fs_plan = center;

%parameters 

disp_L = (disp_o+disp_i)/2;


animated_plots = true;
sim_time = 20; % seconds

NF = N_gait; % NF
foot_size = 0.02; % define the size of the foot, just for the plots
centroid_size = foot_size;
cont=1;
step_duration = 50;
delta = 0.01;
for i=0:step_duration:sim_time/delta + 320 
    fs_timing(cont) = i;
    cont = cont + 1;
end


% general MPC parameters
C = 100;
P = 200;
F = 3;
mpcTimeStep = 0.01; %sampling time
sim_duration = floor(sim_time/mpcTimeStep); % in discrete steps
eta = sqrt(9.8/height);

dsSamples_array = [30];
for i=2:1:N_gait
    if mod(i,2) == 0
        dsSamples_array = [dsSamples_array,30];
    else
        dsSamples_array = [dsSamples_array,30];
    end
end
dsSamples = 30;
% ssSamples = 12;

% ZMP constraint parameters
wx = centroid_size;
wy = centroid_size;
initial_wx = centroid_size;
initial_wy = centroid_size; %/2 + 2*disp_B;



% state initialization
x = disp_C/2;
xd = 0.0;
xz = disp_C/2;
y = 0.0;
yd = 0.0;
yz = 0.0;
current_xfs = fs_plan(1,1);
current_yfs = fs_plan(1,2);
fsCounter = 1;

% state update matrices
ch = cosh(eta*mpcTimeStep);
sh = sinh(eta*mpcTimeStep);
A_upd = [ch, sh/eta, 1-ch; eta*sh, ch, -eta*sh; 0, 0, 1];
B_upd = [mpcTimeStep-sh/eta; 1-ch; mpcTimeStep];

% stored strajectories
x_store = x;
xz_store = xz;
y_store = y;
yz_store = yz;
xfs_store = current_xfs;
yfs_store = current_yfs;

% ZMP constraint constant matrices
Pzmp = tril(ones(C,C)) * mpcTimeStep;
pzmp = ones(C,1);


% compute the centerline for the entire footstep plan
% cl_x = linspace(fs_plan(1,1), fs_plan(1+1,1), dsSamples)';
% cl_y = linspace(fs_plan(1,2), fs_plan(1+1,2), dsSamples)';
cl_x = fs_plan(1,1) * ones(fs_timing(2)-fs_timing(1)-dsSamples, 1);
cl_y = fs_plan(1,2) * ones(fs_timing(2)-fs_timing(1)-dsSamples, 1);
cl_x = [cl_x;linspace(fs_plan(1,1),fs_plan(2,1),dsSamples)'];
cl_y = [cl_y;linspace(fs_plan(1,2),fs_plan(2,2),dsSamples)'];
cl_fscounter = 1;
for i = 2:NF-1
    cl_x = [cl_x; fs_plan(i,1) * ones(fs_timing(2)-dsSamples-fs_timing(1), 1); ...
        linspace(fs_plan(i,1), fs_plan(i+1,1), dsSamples)'];
    cl_y = [cl_y; fs_plan(i,2) * ones(fs_timing(2)-dsSamples-fs_timing(1), 1); ...
        linspace(fs_plan(i,2), fs_plan(i+1,2), dsSamples)'];
end


% useful initializations
predicted_xfs(1) = 0.0;
predicted_yfs(1) = 0.0;
ct = 0;

%handlers 
 h = [zeros(4)];
 p = [zeros(2)];
 
 
 foot_plan1 = foot_plan;
 
 counter = 1; %da 1 a 8 è il counter per i plot della camminata tra triangolo e poligono 
 
 

file1 = fopen('ComTrajectory_walk_phipi4.txt', 'a+');
file2 = fopen('ComVelocity_walk_phipi4.txt', 'a+');

% file per ogni zampa 
file_fl = fopen('foot_fl_walk_phipi4.txt', 'a+');
file_fr = fopen('foot_fr_walk_phipi4.txt', 'a+');
file_rl = fopen('foot_rl_walk_phipi4.txt', 'a+');
file_rr = fopen('foot_rr_walk_phipi4.txt', 'a+');

for j = 1:sim_duration

    % specific dsSample
    dsSamples = dsSamples_array(fsCounter);
    
    % impulsive disturbance
    % do not care for now
    if fsCounter == 2 && ct >= 21 && ct < 36
       bang_x = mpcTimeStep*0.3*0;   
       bang_y = -mpcTimeStep*0.5*0;  %guess m/s da quantificare in N
    else
       bang_x = 0; 
       bang_y = 0;       
    end   
    
    % apply the disturbance     
    xd = xd + bang_x;
    yd = yd + bang_y;
    
    % persistent disturbance, if any   
    w_x(j) = 0;
    w_y(j) = 0.0;
       
    % useful information
    samplesLeftTillNextFootsteps = fs_timing(fsCounter+1) - j;      
      
    %% ZMP constraint
    
    % compute mapping timestep -> predicted footstep number
    mapping = zeros(C,F+1);
    predictedFootstep = 0;
    
    for i = 1:C
        if j + i >= fs_timing(fsCounter + predictedFootstep + 1)
            predictedFootstep = predictedFootstep + 1;
        end
        samplesRemaining = fs_timing(fsCounter + predictedFootstep+1)-(j+i);
        if samplesRemaining > dsSamples  
           mapping(i,predictedFootstep+1) = 1;
           
        else
            mapping(i,predictedFootstep+1) = samplesRemaining/dsSamples;
            mapping(i,predictedFootstep+2) = 1 - samplesRemaining/dsSamples;
        end
    end
    
    A_zmpconstr = [Pzmp, -mapping(:,2:end), zeros(C, C+F); ...
        zeros(C, C+F), Pzmp, -mapping(:,2:end); ...
        -Pzmp, mapping(:,2:end), zeros(C, C+F); ...
        zeros(C, C+F), -Pzmp, mapping(:,2:end)];
    
    b_zmpconstr = [pzmp * (- xz + wx/2) + mapping(:,1) * current_xfs; ...
        pzmp * (- yz + wy/2) + mapping(:,1) * current_yfs; ...
        - pzmp * (- xz - wx/2) - mapping(:,1) * current_xfs; ...
        - pzmp * (- yz - wy/2) - mapping(:,1) * current_yfs];
  

    %% kinematic constraint
    
    
    difference_matrix = eye(F);
    for i = 2:F
        difference_matrix(i,i-1) = -1;
    end
                                            % +++ dx_z, x_f^j, dy_z, y_f^j
    A_kinconstr_component = blkdiag(zeros(C,C), difference_matrix, zeros(C,C), difference_matrix);
    A_kinconstr = [A_kinconstr_component; - A_kinconstr_component];
    b_kinconstr = zeros(4*(C+F), 1);
    for i = 1:F
        b_kinconstr(C+i) = disp_forw;% wkx/2;
        b_kinconstr(3*C+2*F+i) = disp_forw;% wkx/2;
        
        % differentiate between left and right footstep
        if mod(fsCounter+i, 2) == 0
            b_kinconstr(2*C+F+i) = disp_L/2 + disp_L/2;%ell + wky/2;
            b_kinconstr(4*C+3*F+i) = (disp_L/2 + disp_L/2);%-(ell - wky/2);
        else
            b_kinconstr(2*C+F+i) =disp_L/2 + disp_L/2; %-ell + wky/2;
            b_kinconstr(4*C+3*F+i) =-(-disp_L/2 - disp_L/2); %-(-ell - wky/2);
        end
        
    end
    if fsCounter==1 
        b_kinconstr(C+1) = disp_forw_dummy;% wkx/2;
        b_kinconstr(3*C+2*F+1) =disp_forw_dummy;% wkx/2;
        b_kinconstr(2*C+F+1) = disp_L/2 + disp_L/2;%ell + wky/2;
        b_kinconstr(4*C+3*F+1) = (disp_L/2 + disp_L/2);%-(ell - wky/2);
    end
        

    
    % add the contribution from the current footstep
    b_kinconstr(C+1) = b_kinconstr(C+1) + current_xfs;
    b_kinconstr(2*C+F+1) = b_kinconstr(2*C+F+1) + current_yfs;
    b_kinconstr(3*C+2*F+1) = b_kinconstr(3*C+2*F+1) - current_xfs;
    b_kinconstr(4*C+3*F+1) = b_kinconstr(4*C+3*F+1) - current_yfs;
    
    %% stability constraint (numerical anticipative)
    
    % compute anticipative tail
    anticipative_x = exp(-eta*mpcTimeStep*(C+1:P)) * (1-exp(-eta*mpcTimeStep)) *  (cl_x(j+C+1:j+P) - xfs_store(fsCounter)) + ...
        exp(-eta*mpcTimeStep*P) * (cl_x(P) - xfs_store(fsCounter));
    anticipative_y = exp(-eta*mpcTimeStep*(C+1:P)) * (1-exp(-eta*mpcTimeStep)) *  (cl_y(j+C+1:j+P) - yfs_store(fsCounter)) + ...
        exp(-eta*mpcTimeStep*P) * (cl_y(P)  - yfs_store(fsCounter));
     
    lambda = exp(- eta * mpcTimeStep);
    b_exp = exp(- eta * mpcTimeStep * (0:C-1))';
    
    A_stabconstr = [(1/eta)*(1-lambda)/(1-lambda^C) * b_exp' - mpcTimeStep * ones(1,C)*exp(- eta * mpcTimeStep * C), ...
        zeros(1,F+C+F); ...
        zeros(1,C+F), ...
        (1/eta)*(1-lambda)/(1-lambda^C) * b_exp' - mpcTimeStep * ones(1,C)*exp(- eta * mpcTimeStep * C), ...
        zeros(1,F)];
    
    b_stabconstr = [x + xd/eta - xz - anticipative_x ; ...
        y + yd/eta - yz - anticipative_y ];
    
% % STABILITY CONSTRAINT SCRITTE A MODO NOSTRO 
%     % derivata cl_x,cl_y (leo)
%     for i=C:1:P-1
%         cl_xd(j+i)=(cl_x(j+i+1)-cl_x(j+i))/mpcTimeStep;
%         cl_yd(j+i) = (cl_y(j+i+1)-cl_y(j+i))/mpcTimeStep;
%     end
%     
%     % compute anticipative tail 
%     anticipative_x = exp(-eta*mpcTimeStep*(C:P-1)) * (1-exp(-eta*mpcTimeStep)) *  (cl_xd(j+C:j+P-1)');
%     anticipative_y = exp(-eta*mpcTimeStep*(C:P-1)) * (1-exp(-eta*mpcTimeStep)) *  (cl_yd(j+C:j+P-1)');
% 
%     lambda = exp(- eta * mpcTimeStep);
%     b_exp = exp(- eta * mpcTimeStep * (0:C-1))';
% 
%     A_stabconstr = [(1/eta)*(1-lambda) * b_exp', ...
%         zeros(1,F+C+F); ...
%         zeros(1,C+F), ...
%         (1/eta) * (1-lambda) * b_exp', ...
%         zeros(1,F)];
%   
%     b_stabconstr = [x + xd/eta - xz - anticipative_x ; ...
%         y + yd/eta - yz - anticipative_y ];
 

    % cost function
    
    Qzdot = 1;
    Qfootsteps = 1000000000; %100000;
    H_half = blkdiag(Qzdot * eye(C,C), Qfootsteps * eye(F,F));
    H = blkdiag(H_half, H_half);
    
    f = [zeros(C,1); - Qfootsteps * fs_plan(fsCounter+1:fsCounter+F,1); ...
         zeros(C,1); - Qfootsteps * fs_plan(fsCounter+1:fsCounter+F,2)];
    
    %% solve QP
    
    AdsConstr = zeros(2,size(H,1));
    AsdConstr(1,C+1) = (samplesLeftTillNextFootsteps < dsSamples  && fsCounter > 1);
    AsdConstr(1,C+F+C+1) = (samplesLeftTillNextFootsteps < dsSamples  && fsCounter > 1);
    
    bdsConstr = (samplesLeftTillNextFootsteps < dsSamples && fsCounter > 1)*[predicted_xfs(1)-current_xfs;predicted_yfs(1)-current_yfs];
    
    options = optimset('Display','off');
    decisionVariables = quadprog(H, f, [A_zmpconstr;A_kinconstr], [b_zmpconstr;b_kinconstr], [A_stabconstr; 0*AdsConstr], [b_stabconstr; 0*bdsConstr],...
       [],[],[],options);

    predicted_xzd = decisionVariables(1:C);
    predicted_xfs = decisionVariables(C+1:C+F);
    predicted_yzd = decisionVariables(C+F+1:C+F+C);
    predicted_yfs = decisionVariables(C+F+C+1:C+F+C+F);
    
    % update state and store
    
    predstate_x = [x; xd; xz];
    predstate_y = [y; yd; yz];
    
    % compute prediction
    
    up_to = 2;
    
    for i = 1:up_to
        predstate_x = A_upd * predstate_x + B_upd * predicted_xzd(i) + [0;mpcTimeStep;0]*w_x(j);
        predicted_x(i) = predstate_x(1);
        predicted_xd(i) = predstate_x(2);
        predicted_xz(i) = predstate_x(3);
        
        predstate_y = A_upd * predstate_y + B_upd * predicted_yzd(i) + [0;mpcTimeStep;0]*w_y(j);
        predicted_y(i) = predstate_y(1);
        predicted_yd(i) = predstate_y(2);
        predicted_yz(i) = predstate_y(3);
    end
    
    x = predicted_x(1);
    xd = predicted_xd(1);
    xz = predicted_xz(1);
    
    y = predicted_y(1);
    yd = predicted_yd(1) ;
    yz = predicted_yz(1);
    
    % store trajectories
    x_store(end+1) = x;
    xd_store(j) = xd;    
    xz_store(end+1) = xz;
    
    y_store(end+1) = y;
    yd_store(j) = yd;        
    yz_store(end+1) = yz;
    

    %% SECOND QUAD_PROG
 
H1 = eye(2);
% J = (X_f-x_f)^2+(Y_f-y_f)^2
% i passi sono in foot_plan

 if counter == 2 && fsCounter <= 4 % dummy_case
 
    [x_free,y_free,changed] = compute_one_feet_walk(counter,[predicted_xfs(1),predicted_yfs(1)],[foot_plan(fsCounter,1),foot_plan(fsCounter,2),foot_plan(fsCounter,5),foot_plan(fsCounter,6),foot_plan(fsCounter,3),foot_plan(fsCounter,4)],[foot_plan(fsCounter+1,7),foot_plan(fsCounter+1,8)]);
    if changed ==1
        for l=1:1:8
            foot_plan(fsCounter+l,7) = x_free;
            foot_plan(fsCounter+l,8) = y_free;
        end
    changed = 0;
    end
    
    f = [-foot_plan(fsCounter+1,7);-foot_plan(fsCounter+1,8)];
    
    
    A = [0 1 ;
        0 -1 ;
         1 0 ];
     
    b = [disp_o_dummy+foot_plan(fsCounter,8);
        -foot_plan(fsCounter,8)+disp_i_dummy;
        foot_plan(fsCounter,7)+disp_forw_dummy];
    
 elseif counter == 2 && fsCounter > 4 % not_dummy_case
 
    [x_free,y_free,changed] = compute_one_feet_walk(counter,[predicted_xfs(1),predicted_yfs(1)],[foot_plan(fsCounter,1),foot_plan(fsCounter,2),foot_plan(fsCounter,5),foot_plan(fsCounter,6),foot_plan(fsCounter,3),foot_plan(fsCounter,4)],[foot_plan(fsCounter+1,7),foot_plan(fsCounter+1,8)]);
    if changed ==1
        for l=1:1:8
            foot_plan(fsCounter+l,7) = x_free;
            foot_plan(fsCounter+l,8) = y_free;
        end
    changed = 0;
    end
    
    f = [-foot_plan(fsCounter+1,7);-foot_plan(fsCounter+1,8)];
    
    
    A = [0 1 ;
        0 -1 ;
         1 0 ];
     
    b = [disp_o+foot_plan(fsCounter,8);
        -foot_plan(fsCounter,8)+disp_i;
        foot_plan(fsCounter,7)+disp_forw];
    

    

 elseif counter == 4 && fsCounter <= 4 % dummy_case
    [x_free,y_free,changed] = compute_one_feet_walk(counter,[predicted_xfs(1),predicted_yfs(1)],[foot_plan(fsCounter,1),foot_plan(fsCounter,2), foot_plan(fsCounter,5),foot_plan(fsCounter,6), foot_plan(fsCounter,7),foot_plan(fsCounter,8)],[foot_plan(fsCounter+1,3),foot_plan(fsCounter+1,4)]);
    if changed ==1
        for l=1:1:8
            foot_plan(fsCounter+l,3) = x_free;
            foot_plan(fsCounter+l,4) = y_free;
        end
        changed = 0;
    end
    
    f = [-foot_plan(fsCounter+1,3);-foot_plan(fsCounter+1,4)];


    
     A = [0 1 ;
          0 -1 ;
          1 0 ];
             
    b = [disp_i_dummy+foot_plan(fsCounter,4);
        -foot_plan(fsCounter,4)+disp_o_dummy;
        foot_plan(fsCounter,3)+disp_forw_dummy];

 elseif counter == 4 && fsCounter > 4 % not_dummy_case
    
    [x_free,y_free,changed] = compute_one_feet_walk(counter,[predicted_xfs(1),predicted_yfs(1)],[foot_plan(fsCounter,1),foot_plan(fsCounter,2), foot_plan(fsCounter,5),foot_plan(fsCounter,6), foot_plan(fsCounter,7),foot_plan(fsCounter,8)],[foot_plan(fsCounter+1,3),foot_plan(fsCounter+1,4)]);
    if changed ==1
        for l=1:1:8
            foot_plan(fsCounter+l,3) = x_free;
            foot_plan(fsCounter+l,4) = y_free;
        end
        changed = 0;
    end
    
    f = [-foot_plan(fsCounter+1,3);-foot_plan(fsCounter+1,4)];


    
     A = [0 1 ;
          0 -1 ;
          1 0 ];
             
    b = [disp_i+foot_plan(fsCounter,4);
        -foot_plan(fsCounter,4)+disp_o;
        foot_plan(fsCounter,3)+disp_forw];

 elseif counter == 6
    
    [x_free,y_free,changed] = compute_one_feet_walk(counter,[predicted_xfs(1),predicted_yfs(1)],[foot_plan(fsCounter,3),foot_plan(fsCounter,4), foot_plan(fsCounter,7),foot_plan(fsCounter,8),foot_plan(fsCounter,1),foot_plan(fsCounter,2)],[foot_plan(fsCounter+1,5),foot_plan(fsCounter+1,6)]);
    if changed ==1
        for l=1:1:8
            foot_plan(fsCounter+l,5) = x_free;
            foot_plan(fsCounter+l,6) = y_free;
        end
        changed=0;
    end

    f = [-foot_plan(fsCounter+1,5);-foot_plan(fsCounter+1,6)]; 
    
    A = [0 1 ;
         0 -1 ;
         1 0 ];
    b = [disp_i+foot_plan(fsCounter,6);
        -foot_plan(fsCounter,6)+disp_o;
        foot_plan(fsCounter,5)+disp_forw];
    
  elseif counter == 8
      
    [x_free,y_free,changed] = compute_one_feet_walk(counter,[predicted_xfs(1),predicted_yfs(1)],[foot_plan(fsCounter,3),foot_plan(fsCounter,4), foot_plan(fsCounter,7),foot_plan(fsCounter,8), foot_plan(fsCounter,5),foot_plan(fsCounter,6)],[foot_plan(fsCounter+1,1),foot_plan(fsCounter+1,2)]);
    if changed ==1
        for l=1:1:8
            foot_plan(fsCounter+l,1) = x_free;
            foot_plan(fsCounter+l,2) = y_free;
        end
        changed=0;
    end

    f = [-foot_plan(fsCounter+1,1);-foot_plan(fsCounter+1,2)];
    
    A = [0 1 ;
         0 -1 ;
         1 0 ];
     
    b = [disp_o+foot_plan(fsCounter,2);
        -foot_plan(fsCounter,2)+disp_i;
        foot_plan(fsCounter,1)+disp_forw];
end


options = optimset('Display','off');
    
if counter==2
    X = quadprog(H1,f,A,b,[],[],[],[],[],options);
    for l=1:1:8
        foot_plan(fsCounter+l,7)= X(1);
        foot_plan(fsCounter+l,8)= X(2);
    end

elseif counter == 4
    X = quadprog(H1,f,A,b,[],[],[],[],[],options);
    for l=1:1:8
        foot_plan(fsCounter+l,3)= X(1);
        foot_plan(fsCounter+l,4)= X(2);
    end

elseif counter == 6
    X = quadprog(H1,f,A,b,[],[],[],[],[],options);
    for l=1:1:8
        foot_plan(fsCounter+l,5)= X(1);
        foot_plan(fsCounter+l,6)= X(2);
    end

elseif counter == 8
    X = quadprog(H1,f,A,b,[],[],[],[],[],options);
    for l=1:1:8
        foot_plan(fsCounter+l,1)= X(1);
        foot_plan(fsCounter+1,2)= X(2);
    end
end
 
 


 % write CoM on file1.txt (N.B. eliminare ogni volta prima di salvare)
hey = [x_store(j),y_store(j),height];

fprintf(file1, '%d %d %d\n', hey);

% write CoM Velocity on file2.txt (N.B. eliminare ogni volta prima di salvare)
aux = [xd_store(j),yd_store(j),0];

fprintf(file2, '%d %d %d\n', aux);   
    % update some things
    
    ct = ct + 1;
    
    if j+1 >= fs_timing(fsCounter + 1)
        disp('updating footstep counter')
%         j
%         fsCounter
        fsCounter = fsCounter + 1
        counter = counter+1;
        current_xfs = predicted_xfs(1);
        current_yfs = predicted_yfs(1);
        xfs_store(end+1) = predicted_xfs(1);
        yfs_store(end+1) = predicted_yfs(1);
        
        if (fsCounter >= 2 )
        
        % update footstep plam
        fs_plan(:,1) = fs_plan(:,1) + (predicted_xfs(1)-fs_plan(fsCounter,1));
        fs_plan(:,2) = fs_plan(:,2) + (predicted_yfs(1)-fs_plan(fsCounter,2));
        
      
        % compute the centerline for the entire footstep plan
        cl_x = fs_plan(1,1) * ones(fs_timing(2)-fs_timing(1), 1);
        cl_y = fs_plan(1,2) * ones(fs_timing(2)-fs_timing(1), 1);
        cl_fscounter = 1;
        for i = 2:NF-1
            cl_x = [cl_x; fs_plan(i,1) * ones(fs_timing(2)-dsSamples-fs_timing(1), 1); ...
            linspace(fs_plan(i,1), fs_plan(i+1,1), dsSamples)'];
            cl_y = [cl_y; fs_plan(i,2) * ones(fs_timing(2)-dsSamples-fs_timing(1), 1); ...
            linspace(fs_plan(i,2), fs_plan(i+1,2), dsSamples)'];
        end
        
        
        end        
        
        a_count = 0;
        ct = 0;
    end


end
% write feet on file.txt (N.B. eliminare ogni volta prima di salvare)

conteggio = 1;
for i = 1:1:sim_duration/step_duration
    for k = 1:1:step_duration
        
        % sto con quattro zampe poggiate 
        if mod(conteggio,2) == 1
            fprintf(file_fl, '%d %d %d\n',[foot_plan(i,7),foot_plan(i,8), 0.0]);
            fprintf(file_rr, '%d %d %d\n',[foot_plan(i,3),foot_plan(i,4), 0.0]);
            fprintf(file_fr, '%d %d %d\n',[foot_plan(i,5),foot_plan(i,6), 0.0]);
            fprintf(file_rl, '%d %d %d\n',[foot_plan(i,1),foot_plan(i,2), 0.0]);
            
        
        % primo triangolo 1-2-3
        elseif conteggio == 2
            fprintf(file_fl, '%d %d %d\n',[foot_plan(i,7)+(foot_plan(i+1,7)-foot_plan(i,7))/step_duration*k,foot_plan(i,8)+(foot_plan(i+1,8)-foot_plan(i,8))/step_duration*k, -0.000032*k^2+0.0016*k]);
            fprintf(file_rr, '%d %d %d\n',[foot_plan(i,3),foot_plan(i,4), 0.0]);
            fprintf(file_fr, '%d %d %d\n',[foot_plan(i,5),foot_plan(i,6), 0.0]);
            fprintf(file_rl, '%d %d %d\n',[foot_plan(i,1),foot_plan(i,2), 0.0]);
            
            
        % secondo triangolo 1-3-4
        elseif conteggio == 4
            fprintf(file_fl, '%d %d %d\n',[foot_plan(i,7),foot_plan(i,8), 0.0]);
            fprintf(file_rr, '%d %d %d\n',[foot_plan(i,3)+(foot_plan(i+1,3)-foot_plan(i,3))/step_duration*k,foot_plan(i,4)+(foot_plan(i+1,4)-foot_plan(i,4))/step_duration*k, -0.000032*k^2+0.0016*k]);
            fprintf(file_fr, '%d %d %d\n',[foot_plan(i,5),foot_plan(i,6), 0.0]);
            fprintf(file_rl, '%d %d %d\n',[foot_plan(i,1),foot_plan(i,2), 0.0]);
           
        
        % terzo triangolo 1-2-4
        elseif conteggio == 6
            fprintf(file_fl, '%d %d %d\n',[foot_plan(i,7),foot_plan(i,8), 0.0]);
            fprintf(file_rr, '%d %d %d\n',[foot_plan(i,3),foot_plan(i,4), 0.0]);
            fprintf(file_fr, '%d %d %d\n',[foot_plan(i,5)+(foot_plan(i+1,5)-foot_plan(i,5))/step_duration*k,foot_plan(i,6)+(foot_plan(i+1,6)-foot_plan(i,6))/step_duration*k, -0.000032*k^2+0.0016*k]);
            fprintf(file_rl, '%d %d %d\n',[foot_plan(i,1),foot_plan(i,2), 0.0]);
    
        
        % quarto triangolo 2-3-4
        elseif conteggio == 8
            fprintf(file_fl, '%d %d %d\n',[foot_plan(i,7),foot_plan(i,8), 0.0]);
            fprintf(file_rr, '%d %d %d\n',[foot_plan(i,3),foot_plan(i,4), 0.0]);
            fprintf(file_fr, '%d %d %d\n',[foot_plan(i,5),foot_plan(i,6), 0.0]);
            fprintf(file_rl, '%d %d %d\n',[foot_plan(i,1)+(foot_plan(i+1,1)-foot_plan(i,1))/step_duration*k,foot_plan(i,2)+(foot_plan(i+1,2)-foot_plan(i,2))/step_duration*k, -0.000032*k^2+0.0016*k]);
            
        end
    end
    if conteggio == 8
        conteggio = 1;
    else 
        conteggio = conteggio + 1;
    end 
    
end          
fclose(file1);
fclose(file2);
fclose(file_rl);
fclose(file_rr);
fclose(file_fl);
fclose(file_fr);



