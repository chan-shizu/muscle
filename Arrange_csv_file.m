function Arrange_csv_file
warning('on','all');
warning;
muscleNames = loadMuscleName();
muscle_name = [muscleNames{1}];
divisionMuscle = readmatrix("divisionMuscle.csv");
y=divisionMuscle(1);
t=divisionMuscle(2);
h=divisionMuscle(3);

pathDataDisplacementPlane = strcat("muscle\data_", muscle_name, ".csv");%"data_rot_h.csv"
pathDataInitial = strcat("muscle\", muscle_name, "_min_final.csv");%"data0_rot_h.csv"
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
    newDataDisplacementPlane(1,1+3*(i-1)) = newDataDisplacementPlane(1,1+3*(i-1));
    newDataDisplacementPlane(1,2+3*(i-1)) = newDataDisplacementPlane(1,2+3*(i-1));
    newDataDisplacementPlane(1,3+3*(i-1)) = newDataDisplacementPlane(1,3+3*(i-1));
end

csvwrite(pathDataDisplacementPlane, newDataDisplacementPlane);
csvwrite(pathDataInitial, newDataInitial);
tic
end

