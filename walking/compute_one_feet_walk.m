function [x_free,y_free,changed] = compute_one_feet(counter,zmp,fixed,free)
%     syms x y real
%     if counter == 6  
%         A = 1/2*(fixed(1)*fixed(4)-fixed(3)*fixed(2)...
%             +fixed(3)*y-x*fixed(4)...
%             +x*fixed(6)-fixed(5)*y...
%             +fixed(5)*fixed(2)-fixed(1)*fixed(6));
%         C_x = 1/6/A*((fixed(1)+fixed(3))*(fixed(1)*fixed(4)-fixed(3)*fixed(2))...
%             +(fixed(3)+x)*(fixed(3)*y-x*fixed(4))...
%             +(x+fixed(5))*(x*fixed(6)-fixed(5)*y)...
%             +(fixed(5)+fixed(1))*(fixed(5)*fixed(2)-fixed(1)*fixed(6)));
%         
%         C_y = 1/6/A*((fixed(2)+fixed(4))*(fixed(1)*fixed(4)-fixed(3)*fixed(2))...
%             +(fixed(4)+y)*(fixed(3)*y-x*fixed(4))...
%             +(y+fixed(6))*(x*fixed(6)-fixed(5)*y)...
%             +(fixed(6)+fixed(2))*(fixed(5)*fixed(2)-fixed(1)*fixed(6)));
%         
%     elseif counter == 8
%         A = 1/2*(x*fixed(2)-fixed(1)*y...
%             +fixed(1)*fixed(4)-fixed(3)*fixed(2)...
%             +fixed(3)*fixed(6)-fixed(5)*fixed(4)...
%             +fixed(5)*y-x*fixed(6));
%         C_x = 1/6/A*((fixed(1)+x)*(x*fixed(2)-fixed(1)*y)...
%             +(fixed(1)+fixed(3))*(fixed(1)*fixed(4)-fixed(3)*fixed(2))...
%             +(fixed(3)+fixed(5))*(fixed(3)*fixed(6)-fixed(5)*fixed(4))...
%             +(fixed(5)+x)*(fixed(5)*y-x*fixed(6)));
%         C_y = 1/6/A*((fixed(2)+y)*(x*fixed(2)-fixed(1)*y)...
%             +(fixed(2)+fixed(4))*(fixed(1)*fixed(4)-fixed(3)*fixed(2))...
%             +(fixed(4)+fixed(6))*(fixed(3)*fixed(6)-fixed(5)*fixed(4))...
%             +(fixed(6)+y)*(fixed(5)*y-x*fixed(6)));
%         
%     elseif counter == 2
%         A = 1/2*(fixed(1)*fixed(4)-fixed(2)*fixed(3)...
%             +fixed(3)*fixed(6)-fixed(5)*fixed(4)...
%             +fixed(5)*y-x*fixed(6)...
%             +fixed(2)*x-y*fixed(1));
%         
%         C_x = 1/6/A*((fixed(1)+fixed(3))*(fixed(1)*fixed(4)-fixed(2)*fixed(3))...
%             +(fixed(3)+fixed(5))*(fixed(3)*fixed(6)-fixed(5)*fixed(4))...
%             +(fixed(5)+x)*(fixed(5)*y-x*fixed(6))...
%             +(x+fixed(1))*(fixed(2)*x-y*fixed(1)));
%         
%         C_y = 1/6/A*((fixed(2)+fixed(4))*(fixed(1)*fixed(4)-fixed(2)*fixed(3))...
%             +(fixed(4)+fixed(6))*(fixed(3)*fixed(6)-fixed(5)*fixed(4))...
%             +(fixed(6)+y)*(fixed(5)*y-x*fixed(6))...
%             +(y+fixed(2))*(fixed(2)*x-y*fixed(1)));
%         
%     elseif counter == 4
%         A = 1/2*(fixed(1)*y-fixed(2)*x...
%             +x*fixed(4)-fixed(3)*y...
%             +fixed(3)*fixed(6)-fixed(5)*fixed(4)...
%             +fixed(5)*fixed(2)-fixed(1)*fixed(6));
% 
%         C_x =1/6/A* ((fixed(1)+x)*(fixed(1)*y-fixed(2)*x)...
%             +(x+fixed(3))*(x*fixed(4)-fixed(3)*y)...
%             +(fixed(3)+fixed(5))*(fixed(3)*fixed(6)-fixed(5)*fixed(4))...
%             +(fixed(5)+fixed(1))*(fixed(5)*fixed(2)-fixed(1)*fixed(6)));
%         
%         C_y = 1/6/A* ((fixed(2)+y)*(fixed(1)*y-fixed(2)*x)...
%             +(y+fixed(4))*(x*fixed(4)-fixed(3)*y)...
%             +(fixed(4)+fixed(6))*(fixed(3)*fixed(6)-fixed(5)*fixed(4))...
%             +(fixed(6)+fixed(2))*(fixed(5)*fixed(2)-fixed(1)*fixed(6)));
%     end
%     
%     sol = solve([C_x-zmp(1);C_y-zmp(2)],x,y);
%     min_piede = [sol.x(1),sol.y(1)];
%     minimum = 1000;
%     for i=1:1:length(sol.x)
%         dist = (free(1)-sol.x(i))^2+(free(2)-sol.y(i))^2;
%         if dist < minimum
%             minimum=dist;
%             min_piede = [sol.x(i),sol.y(i)];
%         end
%     end
%     x_free = min_piede(1);
%     y_free = min_piede(2);
%     
%     if min_piede(1)-free(1)<0.0001 || min_piede(2)-free(2)<0.0001
%         changed = 0;
%     else
%         changed = 1;
%     end

     syms x y;
     % retta passante per i due piedi fissi
     diagonal_fix = polyfit([fixed(1,1),fixed(1,3)],[fixed(1,2),fixed(1,4)],1);
     equation_fix = diagonal_fix(1)*x + diagonal_fix(2);
    
     % retta con coeff ang opposto passante per lo zmp
     equation_2 = zmp(2)-diagonal_fix(1)*(x-zmp(1));
     
     % intersezione
     temp = solve(equation_fix==equation_2,x);
     y = subs(equation_fix,x,temp);
     
     dist_x = zmp(1)-double(temp);
     dist_y = zmp(2)-double(y);
     
     x_free = free(1)+dist_x;
     y_free = free(2)+dist_y;
         
 
    % Check
    if counter == 2
        quattro_piedi = double([fixed(1,1),fixed(1,2),fixed(1,3),fixed(1,4),fixed(1,5),fixed(1,6),x_free,y_free]);
        if (dist_y ~= 0 || dist_x ~= 0) 
            changed = 1;
            
        else
            changed = 0;
        end
        
    elseif counter == 4
        quattro_piedi = double([fixed(1,1),fixed(1,2),x_free,y_free,fixed(1,3),fixed(1,4),fixed(1,5),fixed(1,6)]);
        if (dist_y ~= 0 || dist_x ~= 0)
            changed = 1;
           
        else
            changed = 0;
        end
        
   elseif counter == 6
        quattro_piedi = double([fixed(1,1),fixed(1,2),fixed(1,3),fixed(1,4),x_free,y_free,fixed(1,5),fixed(1,6)]);
        if (dist_y ~= 0 || dist_x ~= 0)
            changed = 1;

        else
            changed = 0;
        end
        
    elseif counter == 8
        quattro_piedi = double([x_free,y_free,fixed(1,1),fixed(1,2),fixed(1,3),fixed(1,4),fixed(1,5),fixed(1,6)]);
        if (dist_y ~= 0 || dist_x ~= 0)
            changed = 1;

        else
            changed = 0;
        end
        
    end
end