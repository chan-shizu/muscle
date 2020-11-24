clear
warning('on','all');
warning;
checkActivity = true;


L_Percent = readmatrix("muscle_length\Lpercent_rot_27m.csv")
L_Percent = L_Percent(1:10:end,21)
% L_Percent(1:20) = [];
sizeL = size(L_Percent);
timeStep = sizeL(1)*10;
ActivityLevel = readmatrix("muscle_length\Activity_27m.csv")
ActivityLevel = ActivityLevel(:,22)
ActivityLevel(end) = []
sizeA = size(ActivityLevel);
timeStepActivity = sizeA(1);

muscleNames = loadMuscleName();
muscleName = [muscleNames{1}];
divisionMuscle = readmatrix("divisionMuscle.csv");
y=divisionMuscle(1);
t=divisionMuscle(2);
h=divisionMuscle(3);


%%嚥下中のADとPdの伸びとひずみを取得
% for i=1:20
%     
%     AD_length(i) = sqrt((AD_starting_point(1)- stopping_point_x(i))^2+...
%         (AD_starting_point(2) - stopping_point_y(i))^2+...
%         (AD_starting_point(3) - stopping_point_z(i))^2);
%     
%     PD_length(i) = sqrt((PD_starting_point(1) - stopping_point_x(i))^2+...
%         (PD_starting_point(2) - stopping_point_y(i))^2+...
%         (PD_starting_point(3) - stopping_point_z(i))^2);
% end

% distortion_AD = AD_length ./ AD_length(1);
% distortion_PD = PD_length ./ PD_length(1);

% L_Percent = [1 0.99 0.98 0.97 0.96 0.95	0.96 0.97 0.98 0.99	1 0.99 0.98	0.97 0.96 0.95 0.96	0.97 0.985 1]
% L_Percent = [1 0.995 0.99 0.985 0.980 0.975	0.970 0.965 0.960 0.955	0.950 0.955 0.960 0.965 0.970 0.975 0.980 0.985 0.990 1]
% L_Percent = [1.00 0.98 0.96 0.94 0.92 0.90	0.88 0.86 0.84 0.82	0.80 0.82 0.84 0.86 0.88 0.90 0.920 0.95 0.97 1.00];
% L_Percent = [1 0.995 0.99 0.985 0.980 0.975	0.970 0.965 0.960 0.955	0.950 0.945 0.940 0.935 0.930 0.925 0.920 0.915 0.910 0.905]
% L_Percent = [1 0.98 0.96 0.94 0.92 0.90	0.92 0.94 0.96 0.98	1 0.98 0.96 0.94 0.92 0.90 0.92 0.94 0.970 1]
% L_Percent = [1 1.02 1.04 1.06 1.08 1.10	1.12 1.14 1.16 1.18	1.20 1.18 1.16 1.14 1.12 1.10 1.08 1.06 1.03 1.00]
% L_Percent = [1 0.99 0.98 0.97 0.96 0.95 0.94 0.93 0.92 0.91 0.90 0.91 0.92 0.93 0.94 0.95 0.96 0.97 0.985 1]
% L_Percent = [1 0.99 0.98 0.97 0.96 0.95 0.94 0.93 0.92 0.91 0.90 0.90 0.90 0.90 0.90 0.90 0.90 0.90 0.90 0.90 0.90 0.90 0.90 0.90 0.90 0.90 0.90 0.90 0.90 0.90]
% L_Percent = [1 0.95 0.90 0.85 0.80 0.75 0.70 0.65 0.60 0.55 0.50 0.50 0.50 0.50 0.50 0.50 0.50 0.50 0.50 0.50 0.50 0.50 0.50 0.50 0.50 0.50 0.50 0.50 0.50 0.50]
L_Percent = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 ]
sizeL = size(L_Percent);
sizeL(1) = sizeL(2)
timeStep = sizeL(1)*20;

%%強制変位面の座標を作成
pathDataDisplacementPlane = ['muscle\data_', muscleName, '.csv'];%"data_rot_h.csv"
pathDataInitial = ['muscle\', muscleName, '_min_final.csv'];%data0_rot_h.csv"
% dataDisplacementPlane = csvread(pathDataDisplacementPlane, 0, 0);
dataInitial = csvread(pathDataInitial, 0, 0);

if dataInitial(1,3)>dataInitial(y*t+1,3)
    Arrange_csv_file();%min_finalのデータを上面と下面入れ替え
%     dataDisplacementPlane = csvread(pathDataDisplacementPlane, 0, 0);
    dataInitial = csvread(pathDataInitial, 0, 0);
end

for i=1:y*t
    newDataDisplacementPlaneTop(1,1+(i-1)*3:i*3) = dataInitial((h-1)*t*y+i,1:3);
    newDataDisplacementPlaneBottom(1,1+(i-1)*3:i*3) = dataInitial(i,1:3);
end

% x座標とy座標は変わらないので初期座標をそのまま代入
for i=1:y*t
    xValue = newDataDisplacementPlaneTop(1,1+3*(i-1));
    yValue = newDataDisplacementPlaneTop(1,2+3*(i-1));
    newDataDisplacementPlaneTop(1:timeStep,1+3*(i-1)) = xValue;
    newDataDisplacementPlaneTop(1:timeStep,2+3*(i-1)) = yValue;
    
    xValue = newDataDisplacementPlaneBottom(1,1+3*(i-1));
    yValue = newDataDisplacementPlaneBottom(1,2+3*(i-1));
    newDataDisplacementPlaneBottom(1:timeStep,1+3*(i-1)) = xValue;
    newDataDisplacementPlaneBottom(1:timeStep,2+3*(i-1)) = yValue;
end

%強制変位面の座標の時間ステップを多くするためにスプライン補間
distortedLength = (dataInitial(y*t*h,3)-dataInitial(1,3))*(1-(1-L_Percent)/2);
numberStepOriginal = 1:timeStep/sizeL(1):timeStep;
numberStepInterp = 1:1:(timeStep-timeStep/sizeL(1));
interpLength = interp1(numberStepOriginal,distortedLength,numberStepInterp,'spline');
interpLength = movmean(interpLength,5);

for j=1:y*t
    newDataDisplacementPlaneTop(1:timeStep-timeStep/sizeL(1),3*j) = dataInitial(1,3) + interpLength.';
    newDataDisplacementPlaneBottom(1:timeStep-timeStep/sizeL(1),3*j) = dataInitial(h*t*y,3) - interpLength.';
end

for j=1:y*t
    newDataDisplacementPlaneTop(timeStep-timeStep/sizeL(1)+1:timeStep,3*j) = newDataDisplacementPlaneTop(timeStep-timeStep/sizeL(1),3);
    newDataDisplacementPlaneBottom(timeStep-timeStep/sizeL(1)+1:timeStep,3*j) = newDataDisplacementPlaneBottom(timeStep-timeStep/sizeL(1),3);
end

pathDataDisplacementPlaneTop = ['muscle\data_', muscleName, '_top.csv'];%"data_rot_h.csv"
pathDataDisplacementPlaneBottom = ['muscle\data_', muscleName, '_bottom.csv'];%"data_rot_h.csv"
csvwrite(pathDataDisplacementPlaneTop, newDataDisplacementPlaneTop);
csvwrite(pathDataDisplacementPlaneBottom, newDataDisplacementPlaneBottom);

%Activity levelもスプライン補間
if checkActivity == 1
%     numberStepOriginalActivity = 1:timeStepActivity/sizeA(1):timeStepActivity;
%     numberStepInterpActivity = 1:1:(timeStepActivity-timeStepActivity/sizeA(1));
%     interpActivity = interp1(numberStepOriginalActivity,ActivityLevel,numberStepInterp,'spline');
    ActivityLevel2 = ActivityLevel(1:10:end);
    fileNameActivityLevel = ['parameter\',muscleName, '_ActivityLevel.csv'];
    writematrix(ActivityLevel2,fileNameActivityLevel);
%     figure1 = figure();
%     plot(1:sizeA(1),ActivityLevel);
%     pic_name1 = [muscleName,'ActivityLevel.jpg'];
%     figure2 = figure();
%     plot(ActivityLevel2);
%     pic_name2 = [muscleName,'ActivityLevel2.jpg'];
end

figure1 = figure();
plot(1:sizeL(1),L_Percent);

figure1 = figure();
plot(1:timeStep,newDataDisplacementPlaneTop(1:timeStep,3));

figure1 = figure();
plot(1:timeStep,newDataDisplacementPlaneBottom(1:timeStep,3));

