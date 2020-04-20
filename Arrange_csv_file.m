%function area_mive=2Darea():
clc;
clear;
warning('on','all');
warning;
muscleNames = loadMuscleName();
muscle_name = [muscleNames{1}];

y=5;
t=5;
h=6;

timeStep = 1000;

pathDataDisplacementPlane = strcat("C:\Users\kou_0\OneDrive\ドキュメント\研究matlab\muscle\data_", muscle_name, ".csv")%"data_rot_h.csv"
pathDataInitial = strcat("C:\Users\kou_0\OneDrive\ドキュメント\研究matlab\muscle\", muscle_name, "_min_final.csv")%"data0_rot_h.csv"
% file_name_data = strcat("C:\Users\bubbl\Documents\shizuya_M1\DefMuscle_for_TUS_2018\matlab\surface2grid_0608\muscle\min\urabe_arm_data.csv")%"data_rot_h.csv"
% file_name_data0 = strcat("C:\Users\bubbl\Documents\shizuya_M1\DefMuscle_for_TUS_2018\matlab\surface2grid_0608\muscle\min\urabe_arm_data0.csv")%"data0_rot_h.csv"
dataDisplacementPlane = csvread(pathDataDisplacementPlane, 0, 0);
dataInitial = csvread(pathDataInitial, 0, 0);

for i=1:h
    newDataInitial(1+(i-1)*y*t:i*y*t,1:3) = dataInitial(1+(h-i)*y*t:(h-i+1)*y*t,1:3);
end

for i=1:y*t
    newDataDisplacementPlane(1,1+(i-1)*3:i*3) = newDataInitial((h-1)*t*y+i,1:3);
end

for i=1:y*t
    xValue = newDataDisplacementPlane(1,1+3*(i-1));
    yValue = newDataDisplacementPlane(1,2+3*(i-1)); 
    newDataDisplacementPlane(1:timeStep,1+3*(i-1)) = xValue;
    newDataDisplacementPlane(1:timeStep,2+3*(i-1)) = yValue;
end

displacementStep = abs(newDataInitial(1,3) - newDataInitial(1+y*t,3))/(2*timeStep*3/10);

for i=2:timeStep*3/10
    zValue = newDataDisplacementPlane(i-1,3) + displacementStep;
    for j=1:y*t
        newDataDisplacementPlane(i,3*j) = zValue;
    end
end

for i=timeStep*3/10:2*timeStep*3/10
    zValue = newDataDisplacementPlane(i-1,3) - displacementStep;
    for j=1:y*t
        newDataDisplacementPlane(i,3*j) = zValue;
    end
end

for i=2*timeStep*3/10:timeStep
    zValue = newDataDisplacementPlane(i-1,3)
    for j=1:y*t
        newDataDisplacementPlane(i,3*j) = zValue;
    end
end

csvwrite(pathDataDisplacementPlane, newDataDisplacementPlane)
csvwrite(pathDataInitial, newDataInitial)
tic
%プログラムの高速化のために事前割り当て
% datai_x = zeros(timeNum(1),y*t*h);
% datai_y = datai_x;
% datai_z = datai_x;
% con_x = datai_x;
% con_y = datai_x;
% con_z = datai_x;
% Fm_xFm_y = datai_x;
% Fm_z = datai_x;
% Fn_x = datai_x;
% Fn_y = datai_x;
% Fn_z = datai_x;
% Fsl_x = datai_x;
% Fsl_y = datai_x;
% Fsl_z = datai_x;
% Fsr_x = datai_x;
% Fsr_y = datai_x;
% Fsr_z = datai_x;
% FsSum_x = datai_x;
% FsSum_y = datai_x;
% FsSum_z = datai_x;
% Fv_x = datai_x;
% Fv_y = datai_x;
% Fv_z = datai_x;
% vis_x = datai_x;
% vis_y = datai_x;
% vis_z = datai_x;
% vn_x = datai_x;
% vn_y = datai_x;
% vn_z = datai_x;
% Mr_x = datai_x;
% Mr_y = datai_x;
% Mr_z = datai_x;
% F_fib = zeros(timeNum(1),y*t*(h-1));
% f_Lce = F_fib;
% f_Vce = F_fib;
% HillActive = F_fib;
% HillPassive = F_fib;
% length_per = F_fib;
% Vce = F_fib;
% Fs_x = zeros(timeNum(1),y*t*h,6);
% Fs_y = Fs_x;
% Fs_z = Fs_x;
% length_s = zeros(timeNum(1),seNum(1));
% lengthen_s = length_s;
% springF = length_s;
