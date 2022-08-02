%clear; clc; close all; 


%% decision parameters  
N_gait = 100; % NF 
   
disp_A = 0.1; % passo sulle x

phi = pi/4; %lateral motion 


%% robot paramenters 
m = 30.5;        
foot_size = 0.02;% define the size of the foot, just for the plots

height = 0.56;
disp_B = 0.259394;    % 0.2 (nostra iniziale) distanza sulle y dal asse centrale 
disp_C = 0.88;      % 1 (nostra iniziale) lunghezza corpo 
%     ^ Y
%     |  
%     |   
%    [1)..A..>[1')----(4]
%     |                |
%     B                |
% -0--|-------C--------|------> X
%    -B                |
%     |                |
%    [2)--------------(3]..A..>(3']
%     |
%     |
%% admissible region of footplacement ,
disp_i = 0.4;
disp_o = 0.4;
disp_i_dummy = disp_i/2;
disp_o_dummy = disp_o/2;
disp_forw = 0.5;
disp_forw_dummy = disp_forw/2;
%     --------^--------------------------
%     |       |                          |
%     |       |                          |
%     |    disp_o                        |
%     |       |                          |
%     |       |                          |
%     |-------X----------forw----------->|
%     |       |                          |
%     |    disp_i                        |
%     |       |                          |
%     --------v---------------------------
%%
g = 9.81;
disp_vertical = min (disp_i,disp_o);
disp_vertical_dummy = disp_vertical/2;


x_passo = disp_A*cos(phi);
y_passo = disp_A*sin(phi);
x_passo_dummy = disp_A*cos(phi)/2;
y_passo_dummy = disp_A*sin(phi)/2;

counter = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% constraints on the feasibility of the gait
if counter==1
    if y_passo_dummy > disp_vertical_dummy 
        if phi>atan(disp_vertical_dummy/disp_forw_dummy) %oltre il vertice
            y_passo_dummy = disp_vertical_dummy;
            x_passo_dummy = disp_vertical_dummy*cos(phi)/sin(phi);
        else
            x_passo_dummy = disp_forw_dummy;
            y_passo_dummy = disp_forw_dummy*sin(phi)/cos(phi);
        end
        
    elseif x_passo_dummy > disp_forw_dummy
        if phi>atan(disp_vertical_dummy/disp_forw_dummy) %oltre il vertice
            y_passo_dummy = disp_vertical_dummy;
            x_passo_dummy = disp_vertical_dummy*cos(phi)/sin(phi);
        else
            x_passo_dummy = disp_forw_dummy;
            y_passo_dummy = disp_forw_dummy*sin(phi)/cos(phi);
        end
    end
end 
counter = 2;

if y_passo > disp_vertical && counter>=2
    if phi>atan(disp_vertical/disp_forw) %oltre il vertice
        y_passo = disp_vertical;
        x_passo = disp_vertical*cos(phi)/sin(phi);
    else
        x_passo = disp_forw;
        y_passo = disp_forw*sin(phi)/cos(phi);
    end
        
elseif x_passo > disp_forw
    if phi>atan(disp_vertical/disp_forw) %oltre il vertice
        y_passo = disp_vertical;
        x_passo = disp_vertical*cos(phi)/sin(phi);
    else
        x_passo = disp_forw;
        y_passo = disp_forw*sin(phi)/cos(phi);
    end
end

    
   

%%dummy
foot_des_back_left = [zeros(N_gait,1),disp_B*ones(N_gait,1)];
foot_des_back_right = [zeros(N_gait,1),-disp_B*ones(N_gait,1)];
foot_des_front_left = [disp_C*ones(N_gait,1),disp_B*ones(N_gait,1)];
foot_des_front_right = [disp_C*ones(N_gait,1),-disp_B*ones(N_gait,1)];

foot_des_front_left(3,1) = disp_C + x_passo_dummy;
foot_des_front_left(4,1) = foot_des_front_left(3,1);
foot_des_front_left(4,1) = foot_des_front_left(3,1);
foot_des_front_left(5,1) = foot_des_front_left(3,1);


foot_des_back_right(2,1) = foot_des_back_right(1,1);
foot_des_back_right(3,1) = foot_des_back_right(1,1);
foot_des_back_right(4,1) = foot_des_back_right(3,1);
foot_des_back_right(5,1) = foot_des_back_right(4,1)+x_passo_dummy;




foot_des_front_left(3,2) = disp_B + y_passo_dummy;
foot_des_front_left(4,2) = foot_des_front_left(3,2);
foot_des_front_left(5,2) = foot_des_front_left(3,2);


foot_des_back_right(2,2) = foot_des_back_right(1,2);
foot_des_back_right(3,2) = foot_des_back_right(1,2);
foot_des_back_right(4,2) = foot_des_back_right(3,2);
foot_des_back_right(5,2) = foot_des_back_right(4,2)+y_passo_dummy;

% Da finire la gait cambiata..i.e. piede 4 parte per primo 
for j=6:8:N_gait
        
        foot_des_front_right(j,1) = foot_des_front_right(j-1,1);
        foot_des_front_right(j+1,1) = foot_des_front_right(j,1) + x_passo;
        foot_des_front_right(j+2,1) = foot_des_front_right(j+1,1);
        foot_des_front_right(j+3,1) = foot_des_front_right(j+1,1);
        foot_des_front_right(j+4,1) = foot_des_front_right(j+1,1);
        foot_des_front_right(j+5,1) = foot_des_front_right(j+1,1);
        foot_des_front_right(j+6,1) = foot_des_front_right(j+1,1);
        foot_des_front_right(j+7,1) = foot_des_front_right(j+1,1);
        
        foot_des_back_left(j,1) = foot_des_back_left(j-1,1);
        foot_des_back_left(j+1,1) = foot_des_back_left(j,1);
        foot_des_back_left(j+2,1) = foot_des_back_left(j,1);
        foot_des_back_left(j+3,1) = foot_des_back_left(j+2,1)+x_passo;
        foot_des_back_left(j+4,1) = foot_des_back_left(j+3,1);
        foot_des_back_left(j+5,1) = foot_des_back_left(j+3,1);
        foot_des_back_left(j+6,1) = foot_des_back_left(j+3,1);
        foot_des_back_left(j+7,1) = foot_des_back_left(j+3,1);
        
        
        foot_des_front_left(j,1) = foot_des_front_left(j-1,1);
        foot_des_front_left(j+1,1) = foot_des_front_left(j,1);
        foot_des_front_left(j+2,1) = foot_des_front_left(j,1);
        foot_des_front_left(j+3,1) = foot_des_front_left(j,1);
        foot_des_front_left(j+4,1) = foot_des_front_left(j,1);
        foot_des_front_left(j+5,1) = foot_des_front_left(j+4,1)+x_passo;
        foot_des_front_left(j+6,1) = foot_des_front_left(j+5,1);
        foot_des_front_left(j+7,1) = foot_des_front_left(j+5,1);
        
        foot_des_back_right(j,1) = foot_des_back_right(j-1,1);
        foot_des_back_right(j+1,1) = foot_des_back_right(j,1);
        foot_des_back_right(j+2,1) = foot_des_back_right(j,1);
        foot_des_back_right(j+3,1) = foot_des_back_right(j,1);
        foot_des_back_right(j+4,1) = foot_des_back_right(j,1);
        foot_des_back_right(j+5,1) = foot_des_back_right(j,1);
        foot_des_back_right(j+6,1) = foot_des_back_right(j,1);
        foot_des_back_right(j+7,1) = foot_des_back_right(j+6,1)+x_passo;
        
        
        
        
        foot_des_front_right(j,2) = foot_des_front_right(j-1,2);
        foot_des_front_right(j+1,2) = foot_des_front_right(j,2) + y_passo;
        foot_des_front_right(j+2,2) = foot_des_front_right(j+1,2);
        foot_des_front_right(j+3,2) = foot_des_front_right(j+1,2);
        foot_des_front_right(j+4,2) = foot_des_front_right(j+1,2);
        foot_des_front_right(j+5,2) = foot_des_front_right(j+1,2);
        foot_des_front_right(j+6,2) = foot_des_front_right(j+1,2);
        foot_des_front_right(j+7,2) = foot_des_front_right(j+1,2);
        
        foot_des_back_left(j,2) = foot_des_back_left(j-1,2);
        foot_des_back_left(j+1,2) = foot_des_back_left(j,2);
        foot_des_back_left(j+2,2) = foot_des_back_left(j,2);
        foot_des_back_left(j+3,2) = foot_des_back_left(j+2,2)+y_passo;
        foot_des_back_left(j+4,2) = foot_des_back_left(j+3,2);
        foot_des_back_left(j+5,2) = foot_des_back_left(j+3,2);
        foot_des_back_left(j+6,2) = foot_des_back_left(j+3,2);
        foot_des_back_left(j+7,2) = foot_des_back_left(j+3,2);
        
        
        foot_des_front_left(j,2) = foot_des_front_left(j-1,2);
        foot_des_front_left(j+1,2) = foot_des_front_left(j,2);
        foot_des_front_left(j+2,2) = foot_des_front_left(j,2);
        foot_des_front_left(j+3,2) = foot_des_front_left(j,2);
        foot_des_front_left(j+4,2) = foot_des_front_left(j,2);
        foot_des_front_left(j+5,2) = foot_des_front_left(j+4,2)+y_passo;
        foot_des_front_left(j+6,2) = foot_des_front_left(j+5,2);
        foot_des_front_left(j+7,2) = foot_des_front_left(j+5,2);
        
        foot_des_back_right(j,2) = foot_des_back_right(j-1,2);
        foot_des_back_right(j+1,2) = foot_des_back_right(j,2);
        foot_des_back_right(j+2,2) = foot_des_back_right(j,2);
        foot_des_back_right(j+3,2) = foot_des_back_right(j,2);
        foot_des_back_right(j+4,2) = foot_des_back_right(j,2);
        foot_des_back_right(j+5,2) = foot_des_back_right(j,2);
        foot_des_back_right(j+6,2) = foot_des_back_right(j,2);
        foot_des_back_right(j+7,2) = foot_des_back_right(j+6,2)+y_passo;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Possible rotation of the step direction 
theta = 0; %diagonal motion 
R_z_theta = [cos(theta),-sin(theta); sin(theta),cos(theta)];

for i = 1 : N_gait
    foot_des_back_right(i,:) = (R_z_theta*(foot_des_back_right(i,:)'))';
    foot_des_back_left(i,:)  = (R_z_theta*(foot_des_back_left(i,:)'))';
    foot_des_front_right(i,:) = (R_z_theta*(foot_des_front_right(i,:)'))';
    foot_des_front_left(i,:)  = (R_z_theta*(foot_des_front_left(i,:)'))';
end

% all the planned feet
foot_plan = [foot_des_back_left ,foot_des_back_right,foot_des_front_right,foot_des_front_left];  

   
% Compute the center  of  the support poligon
center = zeros(N_gait,2);

center(1,1) = disp_C/2; % first zmp center

count = 1;
for j = 1:8:N_gait-4 
    syms x y real;
    % quadruplo supporto
    for k = 0:2:6
%         fig(j+k) = polyshape([foot_plan(j+k,1),foot_plan(j+k,3),foot_plan(j+k,5),foot_plan(j+k,7)],[foot_plan(j+k,2),foot_plan(j+k,4),foot_plan(j+k,6),foot_plan(j+k,8)]);
%         [x,y] = centroid(fig(j+k));
%         center(j+k,:) = [x,y];
        diagonal_1 = polyfit([foot_plan(j+k,1),foot_plan(j+k,5)],[foot_plan(j+k,2),foot_plan(j+k,6)],1);
        equation_1 = diagonal_1(1)*x + diagonal_1(2);

        diagonal_2 = polyfit([foot_plan(j+k,3),foot_plan(j+k,7)],[foot_plan(j+k,4),foot_plan(j+k,8)],1);
        equation_2 = diagonal_2(1)*x + diagonal_2(2);

        temp = solve([y-equation_1;y-equation_2],x,y);
        center(j+k,:) = [temp.x,temp.y];
    end
    
    % primo triangolo
%     fig(j+1) = polyshape([foot_plan(j+1,1),foot_plan(j+1,3),foot_plan(j+1,5)],[foot_plan(j+1,2),foot_plan(j+1,4),foot_plan(j+1,6)]);
%     [x,y] = centroid(fig(j));
%     center(j+1,:) = [x,y]; % ultima prova
    center(j+1,:) = [center(j,1),center(j,2)];
    
    
    % secondo triangolo
%     fig(j+3) = polyshape([foot_plan(j+3,1),foot_plan(j+3,5),foot_plan(j+3,7)],[foot_plan(j+3,2),foot_plan(j+3,6),foot_plan(j+3,8)]);
%     [x,y] = centroid(fig(j+2));
%     center(j+3,:) = [x,y]; % ultima prova
    center(j+3,:) = [center(j+2,1),center(j+2,2)];
    
    % terzo triangolo
%     fig(j+5) = polyshape([foot_plan(j+5,1),foot_plan(j+5,3),foot_plan(j+5,7)],[foot_plan(j+5,2),foot_plan(j+5,4),foot_plan(j+5,8)]);
%     [x,y] = centroid(fig(j+4));
%     center(j+5,:) = [x,y]; % ultima prova
    center(j+5,:) = [center(j+4,1),center(j+4,2)];
    
    
    % quarto triangolo 
%     fig(j+7) = polyshape([foot_plan(j+7,3),foot_plan(j+7,5),foot_plan(j+7,7)],[foot_plan(j+7,4),foot_plan(j+7,6),foot_plan(j+7,8)]);
%     [x,y] = centroid(fig(j+6));
%     center(j+7,:) = [x,y]; % ultima prova
    center(j+7,:) = [center(j+6,1),center(j+6,2)];

end
