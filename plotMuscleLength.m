warning('on','all');
warning;
clc;
clear;

muscleParam = readtable("muscle_length\MuscleParam_27m.csv",'ReadVariableNames',false);
Lpersent = readmatrix("muscle_length\Lpercent_rot_27m.csv");

muscleSize = size(muscleParam);
text = muscleParam.Var1

for i=1:muscleSize(1)/2
    figure;
    plot(Lpersent(:,2*(i-1)+1));
    title(text(2*(i-1)+1));
end