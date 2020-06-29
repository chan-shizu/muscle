clc;
clear;
warning('on','all');
warning;
muscleNames = loadMuscleName();
muscle_name = [muscleNames{1}];
divisionMuscle = readmatrix("C:\Users\kou_0\OneDrive\ドキュメント\研究matlab\divisionMuscle.csv");
y=divisionMuscle(1);
t=divisionMuscle(2);
h=divisionMuscle(3);
% file_name_data0 = strcat("C:\Users\kou_0\OneDrive\ドキュメント\研究matlab\muscle\", muscle_name, "_min_final.csv");%"data0_rot_h.csv"
% data0=csvread(file_name_data0, 0, 0);%/1000;
data_x_name = ['output\', muscle_name,'_data_x.csv'];
data_y_name = ['output\', muscle_name,'_data_y.csv'];
data_z_name = ['output\', muscle_name,'_data_z.csv'];
datai_x=csvread(data_x_name, 0, 0);
datai_y=csvread(data_y_name, 0, 0);
datai_z=csvread(data_z_name, 0, 0);
dataSize = size(datai_x)

% tetraV1 = abs(dot((data0(26,:)-data0(1,:)),cross((data0(2,:)-data0(1,:)),(data0(6,:)-data0(1,:)))))/6;
% tetraV2 = abs(dot((data0(26,:)-data0(27,:)),cross((data0(6,:)-data0(27,:)),(data0(2,:)-data0(27,:)))))/6;
% tetraV3 = abs(dot((data0(26,:)-data0(27,:)),cross((data0(6,:)-data0(27,:)),(data0(31,:)-data0(27,:)))))/6;
% tetraV4 = abs(dot((data0(32,:)-data0(7,:)),cross((data0(6,:)-data0(7,:)),(data0(2,:)-data0(7,:)))))/6;
% tetraV5 = abs(dot((data0(32,:)-data0(2,:)),cross((data0(27,:)-data0(2,:)),(data0(6,:)-data0(2,:)))))/6;
% tetraV6 = abs(dot((data0(32,:)-data0(6,:)),cross((data0(27,:)-data0(6,:)),(data0(6,:)-data0(31,:)))))/6;

a = 6*4*4*5

for n=1:dataSize(1)
    n
        data0(:,1) = datai_x(n,:)';
        data0(:,2) = datai_y(n,:)';
        data0(:,3) = datai_z(n,:)';

    
    for i=1:h-1
        for j=1:t-1
            for k=1:y-1
                tetraV1(k+(j-1)*(y-1)+(i-1)*(y-1)*(t-1)) = abs(dot((data0(i*t*y+k+(j-1)*y,:)-data0((i-1)*t*y+k+(j-1)*y,:)),cross((data0((i-1)*t*y+1+k+(j-1)*y,:)-data0((i-1)*t*y+k+(j-1)*y,:)),(data0((i-1)*t*y+k+j*y,:)-data0((i-1)*t*y+k+(j-1)*y,:)))))/6;
                tetraV2(k+(j-1)*(y-1)+(i-1)*(y-1)*(t-1)) = abs(dot((data0(i*t*y+k+(j-1)*y,:)-data0(i*t*y+k+(j-1)*y+1,:)),cross((data0((i-1)*t*y+k+j*y,:)-data0(i*t*y+k+(j-1)*y+1,:)),(data0((i-1)*t*y+1+k+(j-1)*y,:)-data0(i*t*y+k+(j-1)*y+1,:)))))/6;
                tetraV3(k+(j-1)*(y-1)+(i-1)*(y-1)*(t-1)) = abs(dot((data0(i*t*y+k+(j-1)*y,:)-data0(i*t*y+k+(j-1)*y+1,:)),cross((data0((i-1)*t*y+k+j*y,:)-data0(i*t*y+k+(j-1)*y+1,:)),(data0(i*t*y+k+j*y,:)-data0(i*t*y+k+(j-1)*y+1,:)))))/6;
                tetraV4(k+(j-1)*(y-1)+(i-1)*(y-1)*(t-1)) = abs(dot((data0(i*t*y+k+j*y+1,:)-data0((i-1)*t*y+k+j*y+1,:)),cross((data0((i-1)*t*y+k+j*y,:)-data0((i-1)*t*y+k+j*y+1,:)),(data0((i-1)*t*y+1+k+(j-1)*y,:)-data0((i-1)*t*y+k+j*y+1,:)))))/6;
                tetraV5(k+(j-1)*(y-1)+(i-1)*(y-1)*(t-1)) = abs(dot((data0(i*t*y+k+j*y+1,:)-data0((i-1)*t*y+1+k+(j-1)*y,:)),cross((data0(i*t*y+k+(j-1)*y+1,:)-data0((i-1)*t*y+1+k+(j-1)*y,:)),(data0((i-1)*t*y+k+j*y,:)-data0((i-1)*t*y+1+k+(j-1)*y,:)))))/6;
                tetraV6(k+(j-1)*(y-1)+(i-1)*(y-1)*(t-1)) = abs(dot((data0(i*t*y+k+j*y+1,:)-data0((i-1)*t*y+k+j*y,:)),cross((data0(i*t*y+k+(j-1)*y+1,:)-data0((i-1)*t*y+k+j*y,:)),(data0((i-1)*t*y+k+j*y,:)-data0(i*t*y+k+(j-1)*y,:)))))/6;
            end
        end
    end
    sumTetraV1 = sum(tetraV1);
    sumTetraV2 = sum(tetraV2);
    sumTetraV3 = sum(tetraV3);
    sumTetraV4 = sum(tetraV4);
    sumTetraV5 = sum(tetraV5);
    sumTetraV6 = sum(tetraV6);
    volume(n) = sumTetraV1 + sumTetraV2 + sumTetraV3 + sumTetraV4 + sumTetraV5 + sumTetraV6;
    
end

plot(volume);


