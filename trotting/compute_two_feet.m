function [quattro_piedi,changed] = compute_two_feet(fsCounter,zmp,fixed)
     
     syms x y;
     diagonal_fix = polyfit([fixed(1,1),fixed(1,3)],[fixed(1,2),fixed(1,4)],1);
     equation_fix = diagonal_fix(1)*x + diagonal_fix(2);
    
     equation_2 = zmp(2)-diagonal_fix(1)*(x-zmp(1));
     
     temp = solve(equation_fix==equation_2,x);
     y = subs(equation_fix,x,temp);
     
     % dist_x = zmp(1)-double(temp.x);
     dist_y = zmp(2)-double(y)
     
     x1 = -(fixed(1,2) + dist_y -zmp(2) - diagonal_fix(1)*zmp(1))/diagonal_fix(1);
     x2 = -(fixed(1,4) + dist_y -zmp(2) - diagonal_fix(1)*zmp(1))/diagonal_fix(1);
         
        
    if mod(fsCounter,2) == 1
        quattro_piedi = double([x2,fixed(1,4)+dist_y, fixed(1,1),fixed(1,2), x1,fixed(1,2)+dist_y, fixed(1,3),fixed(1,4)]);
        if dist_y ~= 0
            changed = 1;
            
        else
            changed = 0;
        end
       
%        hold on
%        plot(x2, fixed(1,4)+dist_y, ['s' 'r'], 'Markersize',15, 'LineWidth',1.5) 
%        plot(fixed(1,1),fixed(1,2), ['s' 'r'], 'Markersize',15, 'LineWidth',1.5) 
%        plot(x1,fixed(1,2)+dist_y, ['s' 'r'], 'Markersize',15, 'LineWidth',1.5) 
%        plot(fixed(1,3),fixed(1,4), ['s' 'r'], 'Markersize',15, 'LineWidth',1.5) 
%        line([x2,x1],[fixed(1,4)+dist_y, fixed(1,2)+dist_y])
        
    elseif mod(fsCounter,2) == 0
        quattro_piedi = [fixed(1,1),fixed(1,2), x2,fixed(1,4)+dist_y, fixed(1,3),fixed(1,4), x1,fixed(1,2)+dist_y];
        if dist_y ~= 0
            changed = 1;
            
        else
            changed = 0;
        end
    end
end