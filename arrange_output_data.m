%function area_mive=2Darea():
clc;
clear;
warning('on','all');
warning;
muscleNames = loadMuscleName();
prompt = '上腕二頭筋なら1, オトガイ舌骨筋なら2, 茎突舌骨筋なら3, 顎二腹筋前腹なら4, 顎二腹筋後腹なら5を押してね';
muscleNumber = inputdlg(prompt,...
             'choose muscle', [1 50])
muscleNumber = str2num(muscleNumber{1})
muscle_name = [muscleNames{muscleNumber,1}];
file_name_data_x = ['output\',muscle_name, '_data_x.csv'];
file_name_data_y = ['output\',muscle_name, '_data_y.csv'];
file_name_data_z = ['output\',muscle_name, '_data_z.csv'];

data_x = readmatrix(file_name_data_x);
data_y = readmatrix(file_name_data_y);
data_z = readmatrix(file_name_data_z);

divisionMuscle = readmatrix("divisionMuscle.csv");
y=divisionMuscle(1);
t=divisionMuscle(2);
h=divisionMuscle(3);

for i=1:h
    if i == 1
        x_arranged(:,1:y*t) = data_x(:,1:y*t);
        y_arranged(:,1:y*t) = data_y(:,1:y*t);
        z_arranged(:,1:y*t) = data_z(:,1:y*t);
    elseif i == h
        x_arranged(:,y*t+(i-2)*(2*y+2*(t-2))+1:2*y*t+(i-2)*(2*y+2*(t-2))) = data_x(:,y*t*(h-1)+1:y*t*h);
        y_arranged(:,y*t+(i-2)*(2*y+2*(t-2))+1:2*y*t+(i-2)*(2*y+2*(t-2))) = data_y(:,y*t*(h-1)+1:y*t*h);
        z_arranged(:,y*t+(i-2)*(2*y+2*(t-2))+1:2*y*t+(i-2)*(2*y+2*(t-2))) = data_z(:,y*t*(h-1)+1:y*t*h);
    else
        x_arranged(:,1+t*y+(i-2)*(2*y+2*(t-2)):y+t*y+(i-2)*(2*y+2*(t-2))) = data_x(:,1+(i-1)*y*t:y+(i-1)*y*t);
        x_arranged(:,1+t*y+y+(i-2)*(2*y+2*(t-2)):2:1+y+t*y+2*(y-3)+(i-2)*(2*y+2*(t-2))) = data_x(:,1+y+(i-1)*y*t:y:1+y*(t-2)+(i-1)*y*t);
        x_arranged(:,2+t*y+y+(i-2)*(2*y+2*(t-2)):2:2+y+t*y+2*(y-3)+(i-2)*(2*y+2*(t-2))) = data_x(:,2*y+(i-1)*y*t:y:y*(t-1)+(i-1)*y*t);
        x_arranged(:,y+t*y+2*(y-2)+1+(i-2)*(2*y+2*(t-2)):2*y+t*y+2*(y-2)+(i-2)*(2*y+2*(t-2))) = data_x(:,y*(t-1)+1+(i-1)*y*t:y*t+(i-1)*y*t);
        y_arranged(:,1+t*y+(i-2)*(2*y+2*(t-2)):y+t*y+(i-2)*(2*y+2*(t-2))) = data_y(:,1+(i-1)*y*t:y+(i-1)*y*t);
        y_arranged(:,1+t*y+y+(i-2)*(2*y+2*(t-2)):2:1+t*y+y+2*(y-3)+(i-2)*(2*y+2*(t-2))) = data_y(:,1+y+(i-1)*y*t:y:1+y*(t-2)+(i-1)*y*t);
        y_arranged(:,2+t*y+y+(i-2)*(2*y+2*(t-2)):2:2+t*y+y+2*(y-3)+(i-2)*(2*y+2*(t-2))) = data_y(:,2*y+(i-1)*y*t:y:y*(t-1)+(i-1)*y*t);
        y_arranged(:,y+t*y+2*(y-2)+1+(i-2)*(2*y+2*(t-2)):2*y+t*y+2*(y-2)+(i-2)*(2*y+2*(t-2))) = data_y(:,y*(t-1)+1+(i-1)*y*t:y*t+(i-1)*y*t);
        z_arranged(:,1+t*y+(i-2)*(2*y+2*(t-2)):y+t*y+(i-2)*(2*y+2*(t-2))) = data_z(:,1+(i-1)*y*t:y+(i-1)*y*t);
        z_arranged(:,1+t*y+y+(i-2)*(2*y+2*(t-2)):2:1+t*y+y+2*(y-3)+(i-2)*(2*y+2*(t-2))) = data_z(:,1+y+(i-1)*y*t:y:1+y*(t-2)+(i-1)*y*t);
        z_arranged(:,2+t*y+y+(i-2)*(2*y+2*(t-2)):2:2+y+t*y+2*(y-3)+(i-2)*(2*y+2*(t-2))) = data_z(:,2*y+(i-1)*y*t:y:y*(t-1)+(i-1)*y*t);
        z_arranged(:,y+t*y+2*(y-2)+1+(i-2)*(2*y+2*(t-2)):2*y+t*y+2*(y-2)+(i-2)*(2*y+2*(t-2))) = data_z(:,y*(t-1)+1+(i-1)*y*t:y*t+(i-1)*y*t);
    end
end

file_name_data_x = ['3ds_max\',muscle_name, '_data_x.csv'];
file_name_data_y = ['3ds_max\',muscle_name, '_data_y.csv'];
file_name_data_z = ['3ds_max\',muscle_name, '_data_z.csv'];

writematrix(x_arranged,file_name_data_x);
writematrix(y_arranged,file_name_data_y);
writematrix(z_arranged,file_name_data_z);

