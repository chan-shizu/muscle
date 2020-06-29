clear
warning('on','all');
warning;
data = readmatrix("\xyz_coordinate.csv");
AD_starting_point = data(1,2:4);
PD_starting_point = data(2,2:4);
time = data(1:20,6);
stopping_point_x = data(1:20,7);
stopping_point_y = data(1:20,8);
stopping_point_z = data(1:20,9);
muscleNames = loadMuscleName();
muscleName = [muscleNames{1}];
divisionMuscle = readmatrix("divisionMuscle.csv");
y=divisionMuscle(1);
t=divisionMuscle(2);
h=divisionMuscle(3);

timeStep = 500;

%%嚥下中のADとPdの伸びとひずみを取得
for i=1:20
    
    AD_length(i) = sqrt((AD_starting_point(1)- stopping_point_x(i))^2+...
        (AD_starting_point(2) - stopping_point_y(i))^2+...
        (AD_starting_point(3) - stopping_point_z(i))^2);
    
    PD_length(i) = sqrt((PD_starting_point(1) - stopping_point_x(i))^2+...
        (PD_starting_point(2) - stopping_point_y(i))^2+...
        (PD_starting_point(3) - stopping_point_z(i))^2);
end

% distortion_AD = AD_length ./ AD_length(1);
% distortion_PD = PD_length ./ PD_length(1);

% distortion_AD = [1 0.99 0.98 0.97 0.96 0.95	0.96 0.97 0.98 0.99	1 0.99 0.98	0.97 0.96 0.95 0.96	0.97 0.985 1]
% distortion_AD = [1 0.995 0.99 0.985 0.980 0.975	0.970 0.965 0.960 0.955	0.950 0.955 0.960 0.965 0.970 0.975 0.980 0.985 0.990 1]
% distortion_AD = [1 0.98 0.96 0.94 0.92 0.90	0.88 0.86 0.84 0.82	0.80 0.82 0.84 0.86 0.88 0.90 0.920 0.95 0.97 1.00];
% distortion_AD = [1 0.995 0.99 0.985 0.980 0.975	0.970 0.965 0.960 0.955	0.950 0.945 0.940 0.935 0.930 0.925 0.920 0.915 0.910 0.905]
% distortion_AD = [1 0.98 0.96 0.94 0.92 0.90	0.92 0.94 0.96 0.98	1 0.98 0.96 0.94 0.92 0.90 0.92 0.94 0.970 1]
% distortion_AD = [1 1.02 1.04 1.06 1.08 1.10	1.12 1.14 1.16 1.18	1.20 1.18 1.16 1.14 1.12 1.10 1.08 1.06 1.03 1.00]
% distortion_AD = [1 0.99 0.98 0.97 0.96 0.95 0.94 0.93 0.92 0.91 0.90 0.91 0.92 0.93 0.94 0.95 0.96 0.97 0.985 1]
distortion_AD = [1 0.99 0.98 0.97 0.96 0.95 0.94 0.93 0.92 0.91 0.90 0.90 0.90 0.90 0.90 0.90 0.90 0.90 0.90 0.90]
% distortion_AD = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 ]

%%強制変位面の座標を作成
pathDataDisplacementPlane = ['muscle\data_', muscleName, '.csv'];%"data_rot_h.csv"
pathDataInitial = ['muscle\', muscleName, '_min_final.csv'];%data0_rot_h.csv"
dataDisplacementPlane = csvread(pathDataDisplacementPlane, 0, 0);
dataInitial = csvread(pathDataInitial, 0, 0);

for i=1:y*t
    newDataDisplacementPlane(1,1+(i-1)*3:i*3) = dataInitial((h-1)*t*y+i,1:3);
end

for i=1:y*t
    xValue = newDataDisplacementPlane(1,1+3*(i-1));
    yValue = newDataDisplacementPlane(1,2+3*(i-1));
    newDataDisplacementPlane(1:timeStep,1+3*(i-1)) = xValue;
    newDataDisplacementPlane(1:timeStep,2+3*(i-1)) = yValue;
end

%強制変位面の座標の時間ステップを多くするためにスプライン補間
distortedLength = (dataInitial(y*t*h,3)-dataInitial(1,3))*distortion_AD;
numberStepOriginal = 1:timeStep*0.8/20:timeStep*0.8;
numberStepInterp = 1:1:(timeStep*0.8-timeStep*0.8/20);
interpLength = interp1(numberStepOriginal,distortedLength,numberStepInterp,'spline');

for j=1:y*t
    newDataDisplacementPlane(1:timeStep*0.8-timeStep*0.8/20,3*j) = dataInitial(1,3) + interpLength.';
end

for j=1:y*t
    newDataDisplacementPlane(timeStep*0.8-timeStep*0.8/20+1:timeStep,3*j) = newDataDisplacementPlane(timeStep*0.8-timeStep*0.8/20,3);
end

pathDataDisplacementPlane = ['muscle\data_', muscleName, '.csv'];%"data_rot_h.csv"
csvwrite(pathDataDisplacementPlane, newDataDisplacementPlane);
figure1 = figure();
plot(1:20,distortion_AD);
pic_name1 = [muscleName,'Distortio_PD.jpg'];
% saveas(figure1,pic_name1);
figure2 = figure();
plot(1:timeStep,newDataDisplacementPlane(1:timeStep,3));
pic_name2 = [muscleName,'NewDataDisplacementPlane.jpg'];
% saveas(figure2,pic_name2);
