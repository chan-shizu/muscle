clear
clc
checkRotation = true;%true or fablse
checkSurface = true;
checkMusclemake = true;
checkArrange = true;
checkParameter = true;
checkSwallowing = true;

%% rotation
if checkRotation == 1
    rotation_traced_data();
end

%% surface2grid_circle_n
if checkSurface == 1
    surface2grid_circle_n_notepc();
end

%% musclemake
if checkMusclemake == 1
    musclemake();
end

%% Arrange csv file
if checkArrange == 1
    Arrange_csv_file();
end

%% parameter
if checkParameter == 1
    parameter_test();
end

%% swallowing, 筋肉の伸縮、活性度を指定
if checkSwallowing == 1
    swallowing_data();
end