clc;
clear;
warning('on','all');
warning;
muscleNames = loadMuscleName();
muscle_name = [muscleNames{1}];
divisionMuscle = readmatrix("divisionMuscle.csv");
y=divisionMuscle(1);
t=divisionMuscle(2);
h=divisionMuscle(3);

file_name_data0 = strcat("muscle\", muscle_name, "_min_final.csv");%"data0_rot_h.csv"
data0=csvread(file_name_data0, 0, 0);%/1000;
dataSize = size(data0)

% a = 6*4*4*5

% data0(:,1) = datai_x(n,:)';
% data0(:,2) = datai_y(n,:)';
% data0(:,3) = datai_z(n,:)';


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
volume = sumTetraV1 + sumTetraV2 + sumTetraV3 + sumTetraV4 + sumTetraV5 + sumTetraV6;

totalm = 1.1*volume*(10^6)
pointm = totalm/(y*t*h);

