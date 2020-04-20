clc;
clear;
warning('on','all');
warning;
muscleNames = loadMuscleName();
muscle_name = [muscleNames{1}]
file_name_readcsv = ['C:\Users\kou_0\OneDrive\ドキュメント\研究matlab\aligned\', muscle_name, '_aligned.csv']
% file_name_readcsv = strcat('C:\Users\bubbl\Documents\shizuya_M1\DefMuscle_for_TUS_2018\matlab\surface2grid_0608\aligned\', muscle_name, '_aligned.csv')
file = csvread(file_name_readcsv,0,3);
file_size = size(file);
if file(1,file_size(2)) == 0
    file(:,file_size(2)) = []
end
file_size = size(file);

%最上層の平均の座標を取得
average_upper_x = mean(file(1:file_size(1),1));
average_upper_y = mean(file(1:file_size(1),2));
average_upper_z = mean(file(1:file_size(1),3));

%最下層の平均の座標を取得
average_lower_x = mean(file(1:file_size(1),file_size(2)-2));
average_lower_y = mean(file(1:file_size(1),file_size(2)-1));
average_lower_z = mean(file(1:file_size(1),file_size(2)));

%xy平面のxz平面の角度を計算
alpha = atan((average_upper_y - average_lower_y)/(average_upper_x - average_lower_x));
beta = atan((average_upper_z - average_lower_z)/(sqrt((average_upper_x - average_lower_x)^2+(average_upper_y - average_lower_y)^2)));
%theta = atan((average_upper_y - average_lower_y)/(average_upper_x - average_lower_x));

center_of_rotation = [average_lower_x-average_upper_x*2,...
    average_lower_y-average_upper_x*2*tan(alpha),...
    average_lower_z-average_upper_x*2*tan(beta)]

% center_of_rotation = [average_lower_x,average_lower_y,average_lower_z];

Coordinate_transformed_y = [cos(-beta) 0 sin(-beta);0 1 0;-sin(-beta) 0 cos(-beta)];%y軸まわりの回転
Coordinate_transformed_z = [cos(-alpha) sin(-alpha) 0 ; sin(-alpha) cos(-alpha) 0; 0 0 1];%z軸まわりの回転

for i=1:file_size(1)
    for j=1:file_size(2)/3
        x_y_z = file(i,1+3*(j-1):3*j) - center_of_rotation;
        Coordinate_transformed_file(i,1+3*(j-1):3*j) = Coordinate_transformed_y * Coordinate_transformed_z...
            * x_y_z.' + center_of_rotation.';
    end
end
file_size2 = size(file)
% file_name_readcsv = ['C:\Users\kou_0\OneDrive\ドキュメント\研究matlab\aligned\', muscle_name, '_aligned.csv']
% name = csvread(file_name_readcsv);
rotated_file(1:file_size2(1),1:file_size2(2)) = 0
rotated_file(:,4:file_size2(2)+3) = Coordinate_transformed_file
csvwrite(file_name_readcsv,rotated_file)

for i=1:file_size2(1)
    for j=1:file_size2(2)/3
        xx(file_size2(2)/3*(i-1)+j) = file(i,3*(j-1)+1);
    end
end
for i=1:file_size2(1)
    for j=1:file_size2(2)/3
        yy(file_size2(2)/3*(i-1)+j) = file(i,3*(j-1)+2);
    end
end
for i=1:file_size2(1)
    for j=1:file_size2(2)/3
        zz(file_size2(2)/3*(i-1)+j) = file(i,3*(j-1)+3);
    end
end

for i=1:file_size2(1)
    for j=1:file_size2(2)/3
        xx_2(file_size2(2)/3*(i-1)+j) = Coordinate_transformed_file(i,3*(j-1)+1);
    end
end
for i=1:file_size2(1)
    for j=1:file_size2(2)/3
        yy_2(file_size2(2)/3*(i-1)+j) = Coordinate_transformed_file(i,3*(j-1)+2);
    end
end
for i=1:file_size2(1)
    for j=1:file_size2(2)/3
        zz_2(file_size2(2)/3*(i-1)+j) = Coordinate_transformed_file(i,3*(j-1)+3);
    end
end
f1 = figure
plot3(xx(1,:), yy(1,:), zz(1,:), 'b.' );
xlabel('x')
ylabel('y')
zlabel('z')
axis equal
f2 = figure
plot3(xx_2(1,:), yy_2(1,:), zz_2(1,:), 'b.' );
xlabel('x')
ylabel('y')
zlabel('z')
axis equal;
