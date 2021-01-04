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
prompt = '上腕二頭筋なら1, オトガイ舌骨筋なら2, 茎突舌骨筋なら3, 顎二腹筋前腹なら4, 顎二腹筋後腹なら5を押してね';
muscleNumber = inputdlg(prompt,...
    'choose muscle', [1 50])
muscleNumber = str2num(muscleNumber{1})
muscleName = [muscleNames{muscleNumber,1}];
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

% L_Percent = [1 0.995 0.99 0.985 0.980 0.975	0.970 0.965 0.960 0.955	0.950 0.945 0.940 0.935 0.930 0.925 0.920 0.915 0.910 0.905]
% L_Percent = [1 0.98 0.96 0.94 0.92 0.90	0.92 0.94 0.96 0.98	1 0.98 0.96 0.94 0.92 0.90 0.92 0.94 0.970 1]
% L_Percent = [1 1.02 1.04 1.06 1.08 1.10	1.12 1.14 1.16 1.18	1.20 1.18 1.16 1.14 1.12 1.10 1.08 1.06 1.03 1.00]
% L_Percent = [1 0.99 0.98 0.97 0.96 0.95 0.94 0.93 0.92 0.91 0.90 0.91 0.92 0.93 0.94 0.95 0.96 0.97 0.985 1]
% L_Percent = [1 0.99 0.98 0.97 0.96 0.95 0.94 0.93 0.92 0.91 0.90 0.90 0.90 0.90 0.90 0.90 0.90 0.90 0.90 0.90 0.90 0.90 0.90 0.90 0.90 0.90 0.90 0.90 0.90 0.90]
% L_Percent = [1 0.90 0.90 0.90 0.90 0.90 0.90 0.90 0.90 0.90 0.90 0.90 0.90 0.90 0.90 0.90 0.90 0.90 0.90 0.90]
% L_Percent = [1 0.99 0.98 0.97 0.96 0.95 0.94 0.93 0.92 0.91 0.90 0.90 0.90 0.90 0.90 0.90 0.90 0.90 0.90 0.90]
L_Percent = [0.90 0.90 0.90 0.90 0.90]
%  L_Percent = [1.1 1.1 1.1 1.1 1.1]
% L_Percent = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 ]

sizeL = size(L_Percent);
sizeL(1) = sizeL(2)
timeStep = sizeL(1)*300;

%%強制変位面の座標を作成
pathDataDisplacementPlane = ['muscle\data_', muscleName, '.csv'];%"data_rot_h.csv"
pathDataInitial = ['muscle\', muscleName, '_min_final.csv'];%data0_rot_h.csv"
% dataDisplacementPlane = csvread(pathDataDisplacementPlane, 0, 0);
dataInitial = csvread(pathDataInitial, 0, 0);

if dataInitial(1,3)<dataInitial(y*t+1,3)
    Arrange_csv_file();%min_finalのデータを上面と下面入れ替え
%     dataDisplacementPlane = csvread(pathDataDisplacementPlane, 0, 0);
    dataInitial = csvread(pathDataInitial, 0, 0);
end

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
distortedLength = (dataInitial(1,3)-dataInitial(y*t*h,3))*L_Percent
numberStepOriginal = 1:timeStep/sizeL(1):timeStep;
numberStepInterp = 1:1:(timeStep-timeStep/sizeL(1));
interpLength = interp1(numberStepOriginal,distortedLength,numberStepInterp,'spline');
interpLength = movmean(interpLength,5);

for j=1:y*t
    newDataDisplacementPlane(1:timeStep-timeStep/sizeL(1),3*j) = dataInitial(1,3) - interpLength.';
end

for j=1:y*t
    newDataDisplacementPlane(timeStep-timeStep/sizeL(1)+1:timeStep,3*j) = newDataDisplacementPlane(timeStep-timeStep/sizeL(1),3);
end

pathDataDisplacementPlane = ['muscle\data_', muscleName, '.csv'];%"data_rot_h.csv"
csvwrite(pathDataDisplacementPlane, newDataDisplacementPlane);

%Activity levelもスプライン補間
if checkActivity == 1
%     numberStepOriginalActivity = 1:timeStepActivity/sizeA(1):timeStepActivity;
%     numberStepInterpActivity = 1:1:(timeStepActivity-timeStepActivity/sizeA(1));
%     interpActivity = interp1(numberStepOriginalActivity,ActivityLevel,numberStepInterp,'spline');
    ActivityLevel2 = ActivityLevel(1:10:end);
    DataSize = size(newDataDisplacementPlane);
    ActivityLevel2(1:DataSize(1)) = 0; 
    fileNameActivityLevel = ['parameter\',muscleName, '_ActivityLevel.csv'];
    writematrix(ActivityLevel2,fileNameActivityLevel);
%     figure1 = figure();
%     plot(1:sizeA(1),ActivityLevel);
%     pic_name1 = [muscleName,'ActivityLevel.jpg'];
%     figure2 = figure();
%     plot(ActivityLevel2);
%     pic_name2 = [muscleName,'ActivityLevel2.jpg'];
end

figure3 = figure();
plot(1:sizeL(1),L_Percent);
pic_name1 = [muscleName,'Distortio_PD.jpg'];
% saveas(figure1,pic_name1);
figure4 = figure();
plot(1:timeStep,newDataDisplacementPlane(1:timeStep,3));
pic_name2 = [muscleName,'NewDataDisplacementPlane.jpg'];
% saveas(figure2,pic_name2);
