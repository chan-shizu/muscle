function no_meaning = rotate(file_name)
file = csvread(file_name_readcsv,0,3);
size = size(file);
size = size(file);

%最上層の平均の座標を取得
average_upper_x = mean(file(1:size(1),1));
average_upper_y = mean(file(1:size(1),2));
average_upper_z = mean(file(1:size(1),3));

%最下層の平均の座標を取得
average_lower_x = mean(file(1:size(1),size(2)-2));
average_lower_y = mean(file(1:size(1),size(2)-1));
average_lower_z = mean(file(1:size(1),size(2)));

%xy平面のxz平面の角度を計算
alpha =  atan((average_upper_y - average_lower_y)/(average_upper_x - average_lower_x));
beta = atan((average_upper_z - average_lower_z)/(sqrt((average_upper_x - average_lower_x)^2+(average_upper_y - average_lower_y)^2)));
%theta = atan((average_upper_y - average_lower_y)/(average_upper_x - average_lower_x));

center_of_rotation = [average_lower_x-average_upper_x*2,...
    average_lower_y-average_upper_x*2*tan(alpha),...
    average_lower_z-average_upper_x*2*tan(beta)]

% center_of_rotation = [average_lower_x,average_lower_y,average_lower_z];

Coordinate_transformed_y = [cos(-beta) 0 sin(-beta);0 1 0;-sin(-beta) 0 cos(-beta)];%y軸まわりの回転
Coordinate_transformed_z = [cos(-alpha) sin(-alpha) 0 ; sin(-alpha) cos(-alpha) 0; 0 0 1];%z軸まわりの回転

for i=1:size(1)
    for j=1:size(2)/3
        x_y_z = file(i,1+3*(j-1):3*j) - center_of_rotation;
        rotated_file(i,1+3*(j-1):3*j) = Coordinate_transformed_y * Coordinate_transformed_z...
            * x_y_z.' + center_of_rotation.';
    end
end

csvwrite(file_name, Coordinate_transformed_file,0,3)
no_meaning = 0;

for i=1:size(1)
    for j=1:size(2)/3
        xx(size(2)/3*(i-1)+j) = file(i,3*(j-1)+1);
    end
end
for i=1:size(1)
    for j=1:size(2)/3
        yy(size(2)/3*(i-1)+j) = file(i,3*(j-1)+2);
    end
end
for i=1:size(1)
    for j=1:size(2)/3
        zz(size(2)/3*(i-1)+j) = file(i,3*(j-1)+3);
    end
end

for i=1:size(1)
    for j=1:size(2)/3
        xx_2(size(2)/3*(i-1)+j) = Coordinate_transformed_file(i,3*(j-1)+1);
    end
end
for i=1:size(1)
    for j=1:size(2)/3
        yy_2(size(2)/3*(i-1)+j) = Coordinate_transformed_file(i,3*(j-1)+2);
    end
end
for i=1:size(1)
    for j=1:size(2)/3
        zz_2(size(2)/3*(i-1)+j) = Coordinate_transformed_file(i,3*(j-1)+3);
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
