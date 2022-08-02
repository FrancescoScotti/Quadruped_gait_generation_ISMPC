clear; clc; close all; 

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

height = 0.56;   
disp_A = 100; % passo sulle x
disp_B = 0.15; % distanza sulle y dal asse centrale 
disp_C = 0.6;      % lunghezza corpo 
disp_i = disp_A/2;
disp_o = disp_A;
DELTA = 0.5*disp_A;
%disp_forw = 2*disp_A+DELTA;
disp_forw = 0.4;
disp_vertical = min (disp_i,min(disp_o,disp_B));


m = 50;        
g = 9.81;

N_gait = 20; % NF
foot_size = 0.05; % define the size of the foot, just for the plots

phi=pi/2;
x_passo = disp_A*cos(phi);
y_passo = disp_A*sin(phi);

if disp_A*sin(phi) > disp_vertical
    if phi>atan(disp_vertical/disp_forw) %oltre il vertice
        y_passo = disp_vertical
        x_passo = disp_vertical*cos(phi)/sin(phi)/2
    else
        x_passo = disp_forw;
        y_passo = disp_forw*sin(phi)/cos(phi);
    end
        
elseif disp_A*cos(phi) > disp_forw
    if phi>atan(disp_vertical/disp_forw) %oltre il vertice
        y_passo = disp_vertical;
        x_passo = disp_vertical*cos(phi)/sin(phi);
    else
        x_passo = disp_forw
        y_passo = disp_forw*sin(phi)/cos(phi)
    end
end

    
    

% initialize feet coordinates 
foot_des_back_left = [zeros(N_gait,1),disp_B*ones(N_gait,1)];
foot_des_back_right = [zeros(N_gait,1),-disp_B*ones(N_gait,1)];
foot_des_front_left = [disp_C*ones(N_gait,1),disp_B*ones(N_gait,1)];
foot_des_front_right = [disp_C*ones(N_gait,1),-disp_B*ones(N_gait,1)];

foot_des_back_left(2,1) = x_passo;
foot_des_front_right(2,1) = disp_C + x_passo;

foot_des_back_left(2,2) = disp_B + y_passo/2;
foot_des_front_right(2,2) = -disp_B + y_passo/2;

for j=3:N_gait
    if mod(j,2)==0
            foot_des_back_left(j,1) = foot_des_back_left(j-1,1) + x_passo;
            foot_des_front_right(j,1) = foot_des_front_right(j-1,1) + x_passo;
            %holding the rest
            foot_des_back_right(j,1) = foot_des_back_right(j-1,1);
            foot_des_front_left(j,1) = foot_des_front_left(j-1,1);
            foot_des_back_left(j,2) = foot_des_back_left(j-1,2) + y_passo;
            foot_des_front_right(j,2) = foot_des_front_right(j-1,2) + y_passo;
            %holding
            foot_des_back_right(j,2) = foot_des_back_right(j-1,2);
            foot_des_front_left(j,2) = foot_des_front_left(j-1,2);

    end
    if mod(j,2)==1
            foot_des_back_right(j,1) = foot_des_back_right(j-1,1) + x_passo;
            foot_des_front_left(j,1) = foot_des_front_left(j-1,1) + x_passo;
            %holding the rest
            foot_des_back_left(j,1) = foot_des_back_left(j-1,1);
            foot_des_front_right(j,1) = foot_des_front_right(j-1,1);
            
            foot_des_back_right(j,2) = foot_des_back_right(j-1,2) + y_passo;
            foot_des_front_left(j,2) = foot_des_front_left(j-1,2) + y_passo;
            
            foot_des_back_left(j,2) = foot_des_back_left(j-1,2);
            foot_des_front_right(j,2) = foot_des_front_right(j-1,2);

    end
end

% Possible rotation of the step direction 
theta = 0;
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

syms x y real;
for k=2:N_gait
    foots = foot_plan(k,:); % all four feet at gait "k"
    
    diagonal_1 = polyfit([foots(:,1),foots(:,5)],[foots(:,2),foots(:,6)],1);
    equation_1 = diagonal_1(1)*x + diagonal_1(2);
    
    diagonal_2 = polyfit([foots(:,3),foots(:,7)],[foots(:,4),foots(:,8)],1);
    equation_2 = diagonal_2(1)*x + diagonal_2(2);
    
    temp = solve([y-equation_1;y-equation_2],x,y);
    center(k,:) = [temp.x,temp.y];
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Rotate the foot steps for the diagonal walk
theta = 0;
R_z_theta = [cos(theta),-sin(theta); sin(theta),cos(theta)];
for i = 1 : N_gait
    foot_des_back_right(i,:) = (R_z_theta*(foot_des_back_right(i,:)'))';
    foot_des_back_left(i,:)  = (R_z_theta*(foot_des_back_left(i,:)'))';
    foot_des_front_right(i,:) = (R_z_theta*(foot_des_front_right(i,:)'))';
    foot_des_front_left(i,:)  = (R_z_theta*(foot_des_front_left(i,:)'))';
end
foot_plan=[foot_des_back_left ,foot_des_back_right,foot_des_front_right,foot_des_front_left]  %fs_plan

   
%% compute the center  of  the support poligon 
center=zeros(N_gait,2);
syms x y real;
for k=1:N_gait
    foots= foot_plan(k,:); %i 4 piendi di ogno quad-supporto 
    %piede 1 3 
    diagonal_1=polyfit([foots(:,1),foots(:,5)],[foots(:,2),foots(:,6)],1);
    equation_1=diagonal_1(1)*x+diagonal_1(2);
    %piede 2 4
    diagonal_2=polyfit([foots(:,3),foots(:,7)],[foots(:,4),foots(:,8)],1);
    equation_2=diagonal_2(1)*x+diagonal_2(2);
    temp = solve([y-equation_1;y-equation_2],x,y);
    center(k,:)=[temp.x,temp.y];
end 
pbaspect([7.5 1 1]);
 p=[zeros(2)]; %poligon 1, solco 1 


for k=1:N_gait
    h=[zeros(4)];
    z=[zeros(1)];
    d=[zeros(2)];
   
    
    foots= foot_plan(k,:);
    corners=zeros(16,2); 
    
    for n=1:4
    %Compute a square for the foot position
            center_foot = [foot_plan(k,2*n-1),foot_plan(k,2*n)];
            corner = center_foot + (R_z_theta*[foot_size/2;foot_size/2])';
            corners(n*4-3,:)=corner;
            square(1,:) = corner;
            
            corner = corner - (R_z_theta*[foot_size;0])';
            corners(n*4-2,:)=corner;
            square(2,:) = corner;
            
            corner = corner - (R_z_theta*[0;foot_size])';
            corners(n*4-1,:)=corner;
            square(3,:) = corner;
            
            corner = corner + (R_z_theta*[foot_size;0])';
            corners(n*4,:)=corner;
            square(4,:) = corner;
            
            corner = corner + (R_z_theta*[0; foot_size])';
            square(5,:) = corner;
            h(n)=plot( square(:,1),square(:,2),'Color','m','lineWidth', 2);
            ylim([-3,3]);
            xlim([-0.5,7]);
            hold on 
            grid
    end
    
    
    %disegno quadratino fittizio
    center_foot = [center(k,1),center(k,2)];
    corner = center_foot + (R_z_theta*[foot_size/2;foot_size/2])';
    square(1,:) = corner;
    corner = corner - (R_z_theta*[foot_size;0])';
    square(2,:) = corner;
    corner = corner - (R_z_theta*[0;foot_size])';
    square(3,:) = corner;
    corner = corner + (R_z_theta*[foot_size;0])';
    square(4,:) = corner;
    corner = corner + (R_z_theta*[0; foot_size])';
    square(5,:) = corner;
    z(1)=plot( square(:,1),square(:,2),'Color','b','lineWidth', 2);
    hold on
    grid
    
    
    d(1)=line([foots(:,1),foots(:,5)],[foots(:,2),foots(:,6)]);
    hold on 
    d(2)=line([foots(:,3),foots(:,7)],[foots(:,4),foots(:,8)]);
    if mod(k,2)==1
        polygon=[corners(2,:);corners(3,:);corners(7,:);corners(8,:);corners(11,:);corners(12,:);corners(16,:);corners(13,:);corners(14,:);corners(1,:);corners(2,:)];
        solco=[corners(14,:);corners(6,:);corners(7,:);corners(8,:);corners(16,:);corners(13,:);corners(14,:)];
    else
        
        polygon=[corners(2,:);corners(6,:);corners(7,:);corners(8,:);corners(11,:);corners(12,:);corners(9,:);corners(13,:);corners(14,:);corners(1,:);corners(2,:)];
        solco=[corners(2,:);corners(3,:);corners(11,:);corners(12,:);corners(9,:);corners(1,:);corners(2,:)];
    end
    
    p(1)=plot(polygon(:,1),polygon(:,2),'Color','r','lineWidth', 2);
    
    ylim([-3,3]);
    xlim([-0.5,7]);
    hold on 
    grid
    pause
    
    set(z(1),'Visible','off');
    

    %swing
    if mod(k,2)==0
        set(h(2),'Visible','off');
        set(h(4),'Visible','off');
        set(d(2),'Visible','off');
        set(p(1),'Visible','off');
        p(2)=plot(solco(:,1),solco(:,2),'Color','r','lineWidth', 2);
        
        center_foot = [(foots(:,1)+foots(:,5))/2,(foots(:,2)+foots(:,6))/2];
        corner = center_foot + (R_z_theta*[foot_size/2;foot_size/2])';
        square(1,:) = corner;
        corner = corner - (R_z_theta*[foot_size;0])';
        square(2,:) = corner;
        corner = corner - (R_z_theta*[0;foot_size])';
        square(3,:) = corner;
        corner = corner + (R_z_theta*[foot_size;0])';
        square(4,:) = corner;
        corner = corner + (R_z_theta*[0; foot_size])';
        square(5,:) = corner;
        z(1)=plot( square(:,1),square(:,2),'Color','b','lineWidth', 2);
        hold on
        grid


    else
        set(h(1),'Visible','off');
        set(h(3),'Visible','off');
        set(d(1),'Visible','off');
        set(p(1),'Visible','off');
        p(2)=plot(solco(:,1),solco(:,2),'Color','r','lineWidth', 2);
        
        center_foot = [(foots(:,3)+foots(:,7))/2,(foots(:,4)+foots(:,8))/2];
        corner = center_foot + (R_z_theta*[foot_size/2;foot_size/2])';
        square(1,:) = corner;
        corner = corner - (R_z_theta*[foot_size;0])';
        square(2,:) = corner;
        corner = corner - (R_z_theta*[0;foot_size])';
        square(3,:) = corner;
        corner = corner + (R_z_theta*[foot_size;0])';
        square(4,:) = corner;
        corner = corner + (R_z_theta*[0; foot_size])';
        square(5,:) = corner;
        z(1)=plot( square(:,1),square(:,2),'Color','b','lineWidth', 2);
        hold on
        grid
    end
         
  
    pause
    
    for l=1:4
        set(h(l),'Visible','off');
    end
    set(z(1),'Visible','off');
    set(d(1),'Visible','off');
    set(d(2),'Visible','off');
    set(p(2),'Visible','off');
    
    
end 
