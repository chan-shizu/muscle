function rotation_traced_data()
clear
clc
warning('on','all');
warning;
muscleNames = loadMuscleName();
muscle_name = [muscleNames{1}];
file_name_readcsv = ['traced\', muscle_name, '.csv'];
% file_name_readcsv = strcat('C:\Users\bubbl\Documents\shizuya_M1\DefMuscle_for_TUS_2018\matlab\surface2grid_0608\aligned\', muscle_name, '_aligned.csv')
file = readmatrix(file_name_readcsv);

% 回転する前のデータを別名で保存, 文字データの取り扱いがよく分からないので，table使用
fileText = readtable(file_name_readcsv);
noRotationFileName = ['traced\', muscle_name, '_no_rotation.csv']
writetable(fileText, noRotationFileName, 'WriteVariableNames',false);

file(:, 1:3) = [];

file_size = size(file);
if file(1,file_size(2)) == 0
    file(:,file_size(2)) = [];
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

%% z軸まわりに回転
%xy平面のxz平面の角度を計算
alpha = atan((average_upper_y - average_lower_y)/(average_upper_x - average_lower_x));
% beta = atan((average_upper_z - average_lower_z)/(sqrt((average_upper_x - average_lower_x)^2+(average_upper_y - average_lower_y)^2)));
        
if average_upper_x - average_lower_x<= 0
            alpha = alpha+ pi;
end

center_of_rotation = [1.0*average_lower_x-0*average_upper_x,...
    1.0*average_lower_y-0*average_upper_y,...
    1.0*average_lower_z-0*average_upper_z]

% center_of_rotation = [(average_lower_x + average_upper_x)/2,...
%     (average_lower_y + average_upper_x*tan(alpha))/2,...
%     (average_lower_z-average_upper_x*2*tan(beta))/2];
% center_of_rotation = [average_lower_x,average_lower_y,average_lower_z];

% Coordinate_transformed_y = [cos(-beta) 0 sin(-beta);0 1 0;-sin(pi/2-beta) 0 cos(pi/2-beta)];%y軸まわりの回転
Coordinate_transformed_z = [cos(alpha) sin(alpha) 0 ; -sin(alpha) cos(alpha) 0; 0 0 1];%z軸まわりの回転
Coordinate_transformed_file = zeros(file_size(1),file_size(2));

for i=1:file_size(1)
    for j=1:file_size(2)/3
        x_y_z = file(i,1+3*(j-1):3*j) - center_of_rotation;
        Coordinate_transformed_file(i,1+3*(j-1):3*j) = Coordinate_transformed_z * x_y_z.' + center_of_rotation.';
    end
end

xx = zeros(1, file_size(1)*file_size(2)/3);
yy= zeros(1, file_size(1)*file_size(2)/3);
zz = zeros(1, file_size(1)*file_size(2)/3);
xx_2 = zeros(1, file_size(1)*file_size(2)/3);
yy_2 = zeros(1, file_size(1)*file_size(2)/3);
zz_2 = zeros(1, file_size(1)*file_size(2)/3);

for i=1:file_size(1)
    for j=1:file_size(2)/3
        xx(file_size(2)/3*(i-1)+j) = file(i,3*(j-1)+1);
    end
end
for i=1:file_size(1)
    for j=1:file_size(2)/3
        yy(file_size(2)/3*(i-1)+j) = file(i,3*(j-1)+2);
    end
end
for i=1:file_size(1)
    for j=1:file_size(2)/3
        zz(file_size(2)/3*(i-1)+j) = file(i,3*(j-1)+3);
    end
end

for i=1:file_size(1)
    for j=1:file_size(2)/3
        xx_2(file_size(2)/3*(i-1)+j) = Coordinate_transformed_file(i,3*(j-1)+1);
    end
end
for i=1:file_size(1)
    for j=1:file_size(2)/3
        yy_2(file_size(2)/3*(i-1)+j) = Coordinate_transformed_file(i,3*(j-1)+2);
    end
end
for i=1:file_size(1)
    for j=1:file_size(2)/3
        zz_2(file_size(2)/3*(i-1)+j) = Coordinate_transformed_file(i,3*(j-1)+3);
    end
end

f1 = figure;
plot3(xx(1,:), yy(1,:), zz(1,:), 'b.' );
xlabel('x')
ylabel('y')
zlabel('z')
axis equal
f2 = figure;
plot3(xx_2(1,:), yy_2(1,:), zz_2(1,:), 'b.' );
xlabel('x')
ylabel('y')
zlabel('z')
axis equal;

%% y軸まわりに回転
%最上層の平均の座標を取得
average_upper_x = mean(Coordinate_transformed_file(1:file_size(1),1));
average_upper_y = mean(Coordinate_transformed_file(1:file_size(1),2));
average_upper_z = mean(Coordinate_transformed_file(1:file_size(1),3));

%最下層の平均の座標を取得
average_lower_x = mean(Coordinate_transformed_file(1:file_size(1),file_size(2)-2));
average_lower_y = mean(Coordinate_transformed_file(1:file_size(1),file_size(2)-1));
average_lower_z = mean(Coordinate_transformed_file(1:file_size(1),file_size(2)));
beta = atan((average_upper_z - average_lower_z)/(sqrt((average_upper_x - average_lower_x)^2+(average_upper_y - average_lower_y)^2)));
        
% if average_upper_x - average_lower_x<= 0
%             beta = localAlpha(j) + pi;
% end

center_of_rotation = [1.0*average_lower_x-0*average_upper_x,...
    1.0*average_lower_y-0*average_upper_y,...
    1.0*average_lower_z-0*average_upper_z]

% center_of_rotation = [(average_lower_x + average_upper_x)/2,...
%     (average_lower_y + average_upper_x*tan(alpha))/2,...
%     (average_lower_z-average_upper_x*2*tan(beta))/2];
% center_of_rotation = [average_lower_x,average_lower_y,average_lower_z];

Coordinate_transformed_y = [cos(pi/2-beta) 0 -sin(pi/2-beta);0 1 0;sin(pi/2-beta) 0 cos(pi/2-beta)];%y軸まわりの回転
% Coordinate_transformed_z = [cos(-alpha) sin(-alpha) 0 ; sin(-alpha) cos(-alpha) 0; 0 0 1];%z軸まわりの回転

for i=1:file_size(1)
    for j=1:file_size(2)/3
        x_y_z = Coordinate_transformed_file(i,1+3*(j-1):3*j) - center_of_rotation;
        Coordinate_transformed_file(i,1+3*(j-1):3*j) = Coordinate_transformed_y * x_y_z.' + center_of_rotation.';
    end
end

for i=1:file_size(1)
    for j=1:file_size(2)/3
        xx3(file_size(2)/3*(i-1)+j) = Coordinate_transformed_file(i,3*(j-1)+1);
    end
end
for i=1:file_size(1)
    for j=1:file_size(2)/3
        yy3(file_size(2)/3*(i-1)+j) = Coordinate_transformed_file(i,3*(j-1)+2);
    end
end
for i=1:file_size(1)
    for j=1:file_size(2)/3
        zz3(file_size(2)/3*(i-1)+j) = Coordinate_transformed_file(i,3*(j-1)+3);
    end
end

f3 = figure;
plot3(xx3(1,:), yy3(1,:), zz3(1,:), 'b.' );
xlabel('x')
ylabel('y')
zlabel('z')
axis equal

CoordinateTransformedFileTable = array2table(Coordinate_transformed_file);
fileText = fileText(:, 1:3);
CoordinateTransformedFileTable = [fileText CoordinateTransformedFileTable];
writetable(CoordinateTransformedFileTable, file_name_readcsv, 'WriteVariableNames',false);
% rotated_file(1:file_size(1),1:file_size(2)) = 0;
% rotated_file(:,4:file_size(2)+3) = Coordinate_transformed_file;
% csvwrite(file_name_readcsv,rotated_file);

%% 層ごとに質点を回転させる．
% for i=1:file_size(2)/3
%     %center of mass points by each layers
%     localAverageX = mean(Coordinate_transformed_file(1:file_size(1),1+3*(i-1)));
%     localAverageY = mean(Coordinate_transformed_file(1:file_size(1),2+3*(i-1)));
%     localAverageZ = mean(Coordinate_transformed_file(1:file_size(1),3+3*(i-1)));
%     
%     localAlpha = zeros(file_size(1));
%     for j = 1:file_size(1)
%         r = sqrt((Coordinate_transformed_file(j,1+3*(i-1)) - localAverageX)^2 + (Coordinate_transformed_file(j,2+3*(i-1)) - localAverageY)^2 + (Coordinate_transformed_file(j,3+3*(i-1)) - localAverageZ)^2);%層の中心と各点からのそれぞれの距離
%         localAlpha(j) = atan((Coordinate_transformed_file(j,2+3*(i-1)) - localAverageY)/(Coordinate_transformed_file(j,1+3*(i-1)) - localAverageX));
%         %             localBeta(j) = atan((Coordinate_transformed_file(j,3+3*(i-1)) - localAverageY)/sqrt((Coordinate_transformed_file(j,1+3*(i-1)) - localAverageX)^2 + (Coordinate_transformed_file(j,2+3*(i-1)) - localAverageY)^2));
%         
%         if Coordinate_transformed_file(j,1+3*(i-1)) - localAverageX <= 0
%             localAlpha(j) = localAlpha(j) + pi;
%         end
%         %         localCoordinateTransformedY = [cos(localBeta(j)) 0 sin(localBeta(j));0 1 0;-sin(localBeta(j)) 0 cos(localBeta(j))];%y軸まわりの回転
%         %         localCoordinateTransformedZ = [cos(localAlpha(j)) sin(localAlpha(j)) 0 ; sin(localAlpha(j)) cos(localAlpha(j)) 0; 0 0 1];%z軸まわりの回転
%         
%         %         x_y_z = Coordinate_transformed_file(j,1+3*(i-1):3*i) - localCenter;
%         Coordinate_transformed_file(j,1+3*(i-1):3*i) = [localAverageX + r*cos(localAlpha(j)),...
%                                                         localAverageY + r*sin(localAlpha(j)),...
%                                                         localAverageZ];
%     end
% end
% 
% xx_3 = zeros(1, file_size(1)*file_size(2)/3);
% yy_3 = zeros(1, file_size(1)*file_size(2)/3);
% zz_3 = zeros(1, file_size(1)*file_size(2)/3);
% for i=1:file_size(1)
%     for j=1:file_size(2)/3
%         xx4(file_size(2)/3*(i-1)+j) = Coordinate_transformed_file(i,3*(j-1)+1);
%     end
% end
% for i=1:file_size(1)
%     for j=1:file_size(2)/3
%         yy4(file_size(2)/3*(i-1)+j) = Coordinate_transformed_file(i,3*(j-1)+2);
%     end
% end
% for i=1:file_size(1)
%     for j=1:file_size(2)/3
%         zz4(file_size(2)/3*(i-1)+j) = Coordinate_transformed_file(i,3*(j-1)+3);
%     end
% end
% 
% rotated_file(1:file_size(1),1:file_size(2)) = 0;
% rotated_file(:,4:file_size(2)+3) = Coordinate_transformed_file;
% csvwrite(file_name_readcsv,rotated_file);
% 
% f4 = figure;
% plot3(xx4(1,:), yy4(1,:), zz4(1,:), 'b.' );
% xlabel('x')
% ylabel('y')
% zlabel('z')
% axis equal
end