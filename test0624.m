clc
clear
dataX = readmatrix("output\test_column_data_x.csv");
dataY = readmatrix("output\test_column_data_y.csv");
dataZ = readmatrix("output\test_column_data_z.csv");
divisionMuscle = readmatrix("divisionMuscle.csv");
y=divisionMuscle(1);
t=divisionMuscle(2);
h=divisionMuscle(3);

muscleNames = loadMuscleName();
prompt = '上腕二頭筋なら1, オトガイ舌骨筋なら2, 茎突舌骨筋なら3, 顎二腹筋前腹なら4, 顎二腹筋後腹なら5を押してね';
muscleNumber = inputdlg(prompt,...
    'choose muscle', [1 50])
muscleNumber = str2num(muscleNumber{1})
muscle_name = [muscleNames{muscleNumber,1}];

%% x
for i=1:h
    if i==1
        newDataX(:,1:y*t) = dataX(:,y*t*(i-1)+1:y*t*(i-1)+y*t);
    elseif i==h
        newDataX(:,y*t+(2*y+2*(y-2))*(h-2)+1:y*t+(2*y+2*(y-2))*(h-2)+y*t) = dataX(:,y*t*(h-1)+1:y*t*h);
    else
        newDataX(:,y*t+(2*y+2*(t-2))*(i-2)+1:y*t+(2*y+2*(t-2))*(i-2)+y) = dataX(:,y*t*(i-1)+1:y*t*(i-1)+y);
        newDataX(:,y*t+(2*y+2*(t-2))*(i-2)+y+2*(t-2)+1:y*t+(2*y+2*(t-2))*(i-2)+y+2*(t-2)+y) = dataX(:,y*t*(i-1)+y*(t-1)+1:y*t*(i-1)+y*t);
    end
end

for i=2:h-1
    for j=2:t-1
        newDataX(:,y*t+(2*y+2*(t-2))*(i-2)+y+2*(j-2)+1) = dataX(:,y*t*(i-1)+y*(j-1)+1);
        newDataX(:,y*t+(2*y+2*(t-2))*(i-2)+y+2*(j-2)+2) = dataX(:,y*t*(i-1)+y*j);
    end
end

%% y
for i=1:h
    if i==1
        newDataY(:,1:y*t) = dataY(:,y*t*(i-1)+1:y*t*(i-1)+y*t);
    elseif i==h
        newDataY(:,y*t+(2*y+2*(y-2))*(h-2)+1:y*t+(2*y+2*(y-2))*(h-2)+y*t) = dataY(:,y*t*(h-1)+1:y*t*h);
    else
        newDataY(:,y*t+(2*y+2*(t-2))*(i-2)+1:y*t+(2*y+2*(t-2))*(i-2)+y) = dataY(:,y*t*(i-1)+1:y*t*(i-1)+y);
        newDataY(:,y*t+(2*y+2*(t-2))*(i-2)+y+2*(t-2)+1:y*t+(2*y+2*(t-2))*(i-2)+y+2*(t-2)+y) = dataY(:,y*t*(i-1)+y*(t-1)+1:y*t*(i-1)+y*t);
    end
end

for i=2:h-1
    for j=2:t-1
        newDataY(:,y*t+(2*y+2*(t-2))*(i-2)+y+2*(j-2)+1) = dataY(:,y*t*(i-1)+y*(j-1)+1);
        newDataY(:,y*t+(2*y+2*(t-2))*(i-2)+y+2*(j-2)+2) = dataY(:,y*t*(i-1)+y*j);
    end
end

%% z
for i=1:h
    if i==1
        newDataZ(:,1:y*t) = dataZ(:,y*t*(i-1)+1:y*t*(i-1)+y*t);
    elseif i==h
        newDataZ(:,y*t+(2*y+2*(y-2))*(h-2)+1:y*t+(2*y+2*(y-2))*(h-2)+y*t) = dataZ(:,y*t*(h-1)+1:y*t*h);
    else
        newDataZ(:,y*t+(2*y+2*(t-2))*(i-2)+1:y*t+(2*y+2*(t-2))*(i-2)+y) = dataZ(:,y*t*(i-1)+1:y*t*(i-1)+y);
        newDataZ(:,y*t+(2*y+2*(t-2))*(i-2)+y+2*(t-2)+1:y*t+(2*y+2*(t-2))*(i-2)+y+2*(t-2)+y) = dataZ(:,y*t*(i-1)+y*(t-1)+1:y*t*(i-1)+y*t);
    end
end

for i=2:h-1
    for j=2:t-1
        newDataZ(:,y*t+(2*y+2*(t-2))*(i-2)+y+2*(j-2)+1) = dataZ(:,y*t*(i-1)+y*(j-1)+1);
        newDataZ(:,y*t+(2*y+2*(t-2))*(i-2)+y+2*(j-2)+2) = dataZ(:,y*t*(i-1)+y*j);
    end
end

newDataPathX = ['output\3dsmax_', muscleName, '_x.csv'];%"data_rot_h.csv"
newDataPathY = ['output\3dsmax_', muscleName, '_y.csv'];%"data_rot_h.csv"
newDataPathZ = ['output\3dsmax_', muscleName, '_z.csv'];%"data_rot_h.csv"
csvwrite(newDataPathX, newDataX);
csvwrite(newDataPathY, newDataY);
csvwrite(newDataPathZ, newDataZ);