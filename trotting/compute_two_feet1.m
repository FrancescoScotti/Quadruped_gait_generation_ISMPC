function [quattro_piedi,changed] = compute_two_feet1(fsCounter,zmp,fixed,free,phi)
     
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
     
         
     % rette piedi fissi con angolazione data da phi (tan(phi))
     if phi==pi/2
         x1 = free(1,1);
         x2 = free(1,3);
         
         y1 = subs(equation_2,x,x1);
         y2 = subs(equation_2,x,x2);
     else
         eq1 = tan(phi)*(x-free(1,1))+free(1,2);
         eq2 = tan(phi)*(x-free(1,3))+free(1,4);
         sol1 = solve(eq1 == equation_2,x);
         sol2 = solve(eq2 == equation_2,x);
         
         y1 = subs(eq1,sol1);
         x1 = sol1;
         y2 = subs(eq2,sol2);
         x2 = sol2;
     end
 
    % Check
    if mod(fsCounter,2) == 1
        quattro_piedi = double([x1,y1, fixed(1,1),fixed(1,2), x2,y2, fixed(1,3),fixed(1,4)]);
        if abs(dist_y ~= 0 || dist_x ~= 0) 
            changed = 1;
            
        else
            changed = 0;
        end
        
    elseif mod(fsCounter,2) == 0
        quattro_piedi = [fixed(1,1),fixed(1,2), x1,y1, fixed(1,3),fixed(1,4), x2,y2];
        if (dist_y ~= 0 || dist_x ~= 0)
            changed = 1;
           
        else
            changed = 0;
        end
    end
end