%% Plot dei file .txt

fileID = fopen('ComTrajectory_walk_phi0.txt','r');
vec = fscanf(fileID,'%f %f %f\n');
vec1 = zeros(2000,3);
cont = 1;
for i=1:3:size(vec(:,1))
    vec1(cont,1) = vec(i);
    vec1(cont,2) = vec(i+1);
    vec1(cont,3) = vec(i+2);
    cont=cont+1;
end
plot([0:0.01:20-0.01],vec1(:,1));