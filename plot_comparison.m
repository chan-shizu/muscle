clc
clear
warning('on','all');
warning;
muscleNames = loadMuscleName();
muscle_name = [muscleNames{1}];
divisionMuscle = readmatrix("divisionMuscle.csv");
y=divisionMuscle(1);
t=divisionMuscle(2);
h=divisionMuscle(3);

MA0point1x = readmatrix("D:\Documents\matlabまとめ\研究matlab\output\gaikei_hikaku\biceps_ver3_MA=0.5_c=50_kamen\biceps_ver3_data_x.csv");
MA0point1y = readmatrix("D:\Documents\matlabまとめ\研究matlab\output\gaikei_hikaku\biceps_ver3_MA=0.5_c=50_kamen\biceps_ver3_data_y.csv");
MA0point1z = readmatrix("D:\Documents\matlabまとめ\研究matlab\output\gaikei_hikaku\biceps_ver3_MA=0.5_c=50_kamen\biceps_ver3_data_z.csv");
MA1x = readmatrix("D:\Documents\matlabまとめ\研究matlab\output\gaikei_hikaku\biceps_ver3_MA=0.5_c=50_zyoumen\biceps_ver3_data_x.csv");
MA1y = readmatrix("D:\Documents\matlabまとめ\研究matlab\output\gaikei_hikaku\biceps_ver3_MA=0.5_c=50_zyoumen\biceps_ver3_data_y.csv");
MA1z = readmatrix("D:\Documents\matlabまとめ\研究matlab\output\gaikei_hikaku\biceps_ver3_MA=0.5_c=50_zyoumen\biceps_ver3_data_z.csv")+0.02;
MA0x = readmatrix("D:\Documents\作業まとめ\2020_0916\活性度変形比較bicep_ver3_c=50\MA=0.1\biceps_ver3_data_x.csv");
MA0y = readmatrix("D:\Documents\作業まとめ\2020_0916\活性度変形比較bicep_ver3_c=50\MA=0.1\biceps_ver3_data_y.csv");
MA0z = readmatrix("D:\Documents\作業まとめ\2020_0916\活性度変形比較bicep_ver3_c=50\MA=0.1\biceps_ver3_data_z.csv");
% MA0point1x_arranged = zeros(1,y*t*h);
% MA0point1z_arranged = zeros(1,y*t*h);

for i=1:h
    MA0point1x_arranged(:,1+(i-1)*(2*y+2*(t-2)):y+(i-1)*(2*y+2*(t-2))) = MA0point1x(:,1+(i-1)*y*t:y+(i-1)*y*t);
    MA0point1x_arranged(:,1+y+(i-1)*(2*y+2*(t-2)):2:1+y+2*(y-3)+(i-1)*(2*y+2*(t-2))) = MA0point1x(:,1+y+(i-1)*y*t:y:1+y*(t-2)+(i-1)*y*t);
    MA0point1x_arranged(:,2+y+(i-1)*(2*y+2*(t-2)):2:2+y+2*(y-3)+(i-1)*(2*y+2*(t-2))) = MA0point1x(:,2*y+(i-1)*y*t:y:y*(t-1)+(i-1)*y*t);
    MA0point1x_arranged(:,y+2*(y-2)+1+(i-1)*(2*y+2*(t-2)):2*y+2*(y-2)+(i-1)*(2*y+2*(t-2))) = MA0point1x(:,y*(t-1)+1+(i-1)*y*t:y*t+(i-1)*y*t);    
    MA0point1y_arranged(:,1+(i-1)*(2*y+2*(t-2)):y+(i-1)*(2*y+2*(t-2))) = MA0point1y(:,1+(i-1)*y*t:y+(i-1)*y*t);
    MA0point1y_arranged(:,1+y+(i-1)*(2*y+2*(t-2)):2:1+y+2*(y-3)+(i-1)*(2*y+2*(t-2))) = MA0point1y(:,1+y+(i-1)*y*t:y:1+y*(t-2)+(i-1)*y*t);
    MA0point1y_arranged(:,2+y+(i-1)*(2*y+2*(t-2)):2:2+y+2*(y-3)+(i-1)*(2*y+2*(t-2))) = MA0point1y(:,2*y+(i-1)*y*t:y:y*(t-1)+(i-1)*y*t);
    MA0point1y_arranged(:,y+2*(y-2)+1+(i-1)*(2*y+2*(t-2)):2*y+2*(y-2)+(i-1)*(2*y+2*(t-2))) = MA0point1y(:,y*(t-1)+1+(i-1)*y*t:y*t+(i-1)*y*t);   
    MA0point1z_arranged(:,1+(i-1)*(2*y+2*(t-2)):y+(i-1)*(2*y+2*(t-2))) = MA0point1z(:,1+(i-1)*y*t:y+(i-1)*y*t);
    MA0point1z_arranged(:,1+y+(i-1)*(2*y+2*(t-2)):2:1+y+2*(y-3)+(i-1)*(2*y+2*(t-2))) = MA0point1z(:,1+y+(i-1)*y*t:y:1+y*(t-2)+(i-1)*y*t);
    MA0point1z_arranged(:,2+y+(i-1)*(2*y+2*(t-2)):2:2+y+2*(y-3)+(i-1)*(2*y+2*(t-2))) = MA0point1z(:,2*y+(i-1)*y*t:y:y*(t-1)+(i-1)*y*t);
    MA0point1z_arranged(:,y+2*(y-2)+1+(i-1)*(2*y+2*(t-2)):2*y+2*(y-2)+(i-1)*(2*y+2*(t-2))) = MA0point1z(:,y*(t-1)+1+(i-1)*y*t:y*t+(i-1)*y*t);
    
    MA1x_arranged(:,1+(i-1)*(2*y+2*(t-2)):y+(i-1)*(2*y+2*(t-2))) = MA1x(:,1+(i-1)*y*t:y+(i-1)*y*t);
    MA1x_arranged(:,1+y+(i-1)*(2*y+2*(t-2)):2:1+y+2*(y-3)+(i-1)*(2*y+2*(t-2))) = MA1x(:,1+y+(i-1)*y*t:y:1+y*(t-2)+(i-1)*y*t);
    MA1x_arranged(:,2+y+(i-1)*(2*y+2*(t-2)):2:2+y+2*(y-3)+(i-1)*(2*y+2*(t-2))) = MA1x(:,2*y+(i-1)*y*t:y:y*(t-1)+(i-1)*y*t);
    MA1x_arranged(:,y+2*(y-2)+1+(i-1)*(2*y+2*(t-2)):2*y+2*(y-2)+(i-1)*(2*y+2*(t-2))) = MA1x(:,y*(t-1)+1+(i-1)*y*t:y*t+(i-1)*y*t);    
    MA1y_arranged(:,1+(i-1)*(2*y+2*(t-2)):y+(i-1)*(2*y+2*(t-2))) = MA1y(:,1+(i-1)*y*t:y+(i-1)*y*t);
    MA1y_arranged(:,1+y+(i-1)*(2*y+2*(t-2)):2:1+y+2*(y-3)+(i-1)*(2*y+2*(t-2))) = MA1y(:,1+y+(i-1)*y*t:y:1+y*(t-2)+(i-1)*y*t);
    MA1y_arranged(:,2+y+(i-1)*(2*y+2*(t-2)):2:2+y+2*(y-3)+(i-1)*(2*y+2*(t-2))) = MA1y(:,2*y+(i-1)*y*t:y:y*(t-1)+(i-1)*y*t);
    MA1y_arranged(:,y+2*(y-2)+1+(i-1)*(2*y+2*(t-2)):2*y+2*(y-2)+(i-1)*(2*y+2*(t-2))) = MA1y(:,y*(t-1)+1+(i-1)*y*t:y*t+(i-1)*y*t);   
    MA1z_arranged(:,1+(i-1)*(2*y+2*(t-2)):y+(i-1)*(2*y+2*(t-2))) = MA1z(:,1+(i-1)*y*t:y+(i-1)*y*t);
    MA1z_arranged(:,1+y+(i-1)*(2*y+2*(t-2)):2:1+y+2*(y-3)+(i-1)*(2*y+2*(t-2))) = MA1z(:,1+y+(i-1)*y*t:y:1+y*(t-2)+(i-1)*y*t);
    MA1z_arranged(:,2+y+(i-1)*(2*y+2*(t-2)):2:2+y+2*(y-3)+(i-1)*(2*y+2*(t-2))) = MA1z(:,2*y+(i-1)*y*t:y:y*(t-1)+(i-1)*y*t);
    MA1z_arranged(:,y+2*(y-2)+1+(i-1)*(2*y+2*(t-2)):2*y+2*(y-2)+(i-1)*(2*y+2*(t-2))) = MA1z(:,y*(t-1)+1+(i-1)*y*t:y*t+(i-1)*y*t);
    
    MA0x_arranged(:,1+(i-1)*(2*y+2*(t-2)):y+(i-1)*(2*y+2*(t-2))) = MA0x(:,1+(i-1)*y*t:y+(i-1)*y*t);
    MA0x_arranged(:,1+y+(i-1)*(2*y+2*(t-2)):2:1+y+2*(y-3)+(i-1)*(2*y+2*(t-2))) = MA0x(:,1+y+(i-1)*y*t:y:1+y*(t-2)+(i-1)*y*t);
    MA0x_arranged(:,2+y+(i-1)*(2*y+2*(t-2)):2:2+y+2*(y-3)+(i-1)*(2*y+2*(t-2))) = MA0x(:,2*y+(i-1)*y*t:y:y*(t-1)+(i-1)*y*t);
    MA0x_arranged(:,y+2*(y-2)+1+(i-1)*(2*y+2*(t-2)):2*y+2*(y-2)+(i-1)*(2*y+2*(t-2))) = MA0x(:,y*(t-1)+1+(i-1)*y*t:y*t+(i-1)*y*t); 
    MA0y_arranged(:,1+(i-1)*(2*y+2*(t-2)):y+(i-1)*(2*y+2*(t-2))) = MA0y(:,1+(i-1)*y*t:y+(i-1)*y*t);
    MA0y_arranged(:,1+y+(i-1)*(2*y+2*(t-2)):2:1+y+2*(y-3)+(i-1)*(2*y+2*(t-2))) = MA0y(:,1+y+(i-1)*y*t:y:1+y*(t-2)+(i-1)*y*t);
    MA0y_arranged(:,2+y+(i-1)*(2*y+2*(t-2)):2:2+y+2*(y-3)+(i-1)*(2*y+2*(t-2))) = MA0y(:,2*y+(i-1)*y*t:y:y*(t-1)+(i-1)*y*t);
    MA0y_arranged(:,y+2*(y-2)+1+(i-1)*(2*y+2*(t-2)):2*y+2*(y-2)+(i-1)*(2*y+2*(t-2))) = MA0y(:,y*(t-1)+1+(i-1)*y*t:y*t+(i-1)*y*t);   
    MA0z_arranged(:,1+(i-1)*(2*y+2*(t-2)):y+(i-1)*(2*y+2*(t-2))) = MA0z(:,1+(i-1)*y*t:y+(i-1)*y*t);
    MA0z_arranged(:,1+y+(i-1)*(2*y+2*(t-2)):2:1+y+2*(y-3)+(i-1)*(2*y+2*(t-2))) = MA0z(:,1+y+(i-1)*y*t:y:1+y*(t-2)+(i-1)*y*t);
    MA0z_arranged(:,2+y+(i-1)*(2*y+2*(t-2)):2:2+y+2*(y-3)+(i-1)*(2*y+2*(t-2))) = MA0z(:,2*y+(i-1)*y*t:y:y*(t-1)+(i-1)*y*t);
    MA0z_arranged(:,y+2*(y-2)+1+(i-1)*(2*y+2*(t-2)):2*y+2*(y-2)+(i-1)*(2*y+2*(t-2))) = MA0z(:,y*(t-1)+1+(i-1)*y*t:y*t+(i-1)*y*t);
    end
% sizeData = size(MA0point1x);
plot(MA0point1y_arranged(end,:),MA0point1z_arranged(end,:),'r.');
hold on
plot(MA1y_arranged(end,:),MA1z_arranged(end,:),'g.');
plot(MA1y_arranged(1,:),MA1z_arranged(1,:),'b.');
axis equal