clc
clear
% muscleNames = loadMuscleName();
% prompt = '上腕二頭筋なら1, オトガイ舌骨筋なら2, 茎突舌骨筋なら3, 顎二腹筋前腹なら4, 顎二腹筋後腹なら5を押してね';
% muscleNumber = inputdlg(prompt,...
%     'choose muscle', [1 50])
% muscleNumber = str2num(muscleNumber{1})
% muscle_name = [muscleNames{muscleNumber,1}];
muscle_name = 'biceps_ver3'

min_n_path = strcat('muscle\',muscle_name, '_min_final.csv');
file_n_final_input = readmatrix(min_n_path);

% file_name_readcsv = strcat('aligned\','biceps_ver2', '_aligned.csv');
% file = readmatrix(file_name_readcsv);
% fileNum=size(file);

divisionMuscle = readmatrix("divisionMuscle.csv");
y=divisionMuscle(1);
t=divisionMuscle(2);
h=divisionMuscle(3);

file_data0=zeros((y*2+2*(t-2))*h,3);
file_n=ones(y*2+2*(t-2),3*h);

for layer =1:h
    file_data0((layer-1)*(2*y+2*(t-2))+1:(layer-1)*(2*y+2*(t-2))+6,1:3) = file_n_final_input((layer-1)*y*t+1:(layer-1)*y*t+6,1:3);
    file_data0((layer-1)*(2*y+2*(t-2))+7:(layer-1)*(2*y+2*(t-2))+8,1:3) = file_n_final_input((layer-1)*y*t+10:(layer-1)*y*t+11,1:3);
    file_data0((layer-1)*(2*y+2*(t-2))+9:(layer-1)*(2*y+2*(t-2))+10,1:3) = file_n_final_input((layer-1)*y*t+15:(layer-1)*y*t+16,1:3);
    file_data0((layer-1)*(2*y+2*(t-2))+11,1:3) = file_n_final_input((layer-1)*y*t+20,1:3);
    file_data0((layer-1)*(2*y+2*(t-2))+12:(layer-1)*(2*y+2*(t-2))+16,1:3) = file_n_final_input((layer-1)*y*t+21:(layer-1)*y*t+25,1:3);
end

file_data0_arranged = file_data0

% 最後から二番目の層の質点の順番を入れ替える
% file_data0_arranged((y*2+2*(t-2))*10+1,:) = file_data0((y*2+2*(t-2))*10+2,:);
% file_data0_arranged((y*2+2*(t-2))*10+2,:) = file_data0((y*2+2*(t-2))*10+3,:);
% file_data0_arranged((y*2+2*(t-2))*10+3,:) = file_data0((y*2+2*(t-2))*10+4,:);
% file_data0_arranged((y*2+2*(t-2))*10+4,:) = file_data0((y*2+2*(t-2))*10+5,:);
% file_data0_arranged((y*2+2*(t-2))*10+5,:) = file_data0((y*2+2*(t-2))*10+7,:);
% file_data0_arranged((y*2+2*(t-2))*10+6,:) = file_data0((y*2+2*(t-2))*10+1,:);
% file_data0_arranged((y*2+2*(t-2))*10+7,:) = file_data0((y*2+2*(t-2))*10+9,:);
% file_data0_arranged((y*2+2*(t-2))*10+8,:) = file_data0((y*2+2*(t-2))*10+6,:);
% file_data0_arranged((y*2+2*(t-2))*10+9,:) = file_data0((y*2+2*(t-2))*10+11,:);
% file_data0_arranged((y*2+2*(t-2))*10+10,:) = file_data0((y*2+2*(t-2))*10+8,:);
% file_data0_arranged((y*2+2*(t-2))*10+11,:) = file_data0((y*2+2*(t-2))*10+16,:);
% file_data0_arranged((y*2+2*(t-2))*10+12,:) = file_data0((y*2+2*(t-2))*10+10,:);
% file_data0_arranged((y*2+2*(t-2))*10+13,:) = file_data0((y*2+2*(t-2))*10+12,:);
% file_data0_arranged((y*2+2*(t-2))*10+14,:) = file_data0((y*2+2*(t-2))*10+13,:);
% file_data0_arranged((y*2+2*(t-2))*10+15,:) = file_data0((y*2+2*(t-2))*10+14,:);
% file_data0_arranged((y*2+2*(t-2))*10+16,:) = file_data0((y*2+2*(t-2))*10+15,:);

%最後の層の質点の順番を入れ替える
% file_data0_arranged((y*2+2*(t-2))*11+1,:) = file_data0((y*2+2*(t-2))*11+2,:);
% file_data0_arranged((y*2+2*(t-2))*11+2,:) = file_data0((y*2+2*(t-2))*11+3,:);
% file_data0_arranged((y*2+2*(t-2))*11+3,:) = file_data0((y*2+2*(t-2))*11+4,:);
% file_data0_arranged((y*2+2*(t-2))*11+4,:) = file_data0((y*2+2*(t-2))*11+5,:);
% file_data0_arranged((y*2+2*(t-2))*11+5,:) = file_data0((y*2+2*(t-2))*11+7,:);
% file_data0_arranged((y*2+2*(t-2))*11+6,:) = file_data0((y*2+2*(t-2))*11+1,:);
% file_data0_arranged((y*2+2*(t-2))*11+7,:) = file_data0((y*2+2*(t-2))*11+9,:);
% file_data0_arranged((y*2+2*(t-2))*11+8,:) = file_data0((y*2+2*(t-2))*11+6,:);
% file_data0_arranged((y*2+2*(t-2))*11+9,:) = file_data0((y*2+2*(t-2))*11+11,:);
% file_data0_arranged((y*2+2*(t-2))*11+10,:) = file_data0((y*2+2*(t-2))*11+8,:);
% file_data0_arranged((y*2+2*(t-2))*11+11,:) = file_data0((y*2+2*(t-2))*11+16,:);
% file_data0_arranged((y*2+2*(t-2))*11+12,:) = file_data0((y*2+2*(t-2))*11+10,:);
% file_data0_arranged((y*2+2*(t-2))*11+13,:) = file_data0((y*2+2*(t-2))*11+12,:);
% file_data0_arranged((y*2+2*(t-2))*11+14,:) = file_data0((y*2+2*(t-2))*11+13,:);
% file_data0_arranged((y*2+2*(t-2))*11+15,:) = file_data0((y*2+2*(t-2))*11+14,:);
% file_data0_arranged((y*2+2*(t-2))*11+16,:) = file_data0((y*2+2*(t-2))*11+15,:);

file_data0_arranged((y*2+2*(t-2))*11+1,:) = file_data0((y*2+2*(t-2))*11+6,:);
file_data0_arranged((y*2+2*(t-2))*11+2,:) = file_data0((y*2+2*(t-2))*11+1,:);
file_data0_arranged((y*2+2*(t-2))*11+3,:) = file_data0((y*2+2*(t-2))*11+2,:);
file_data0_arranged((y*2+2*(t-2))*11+4,:) = file_data0((y*2+2*(t-2))*11+3,:);
file_data0_arranged((y*2+2*(t-2))*11+5,:) = file_data0((y*2+2*(t-2))*11+4,:);
file_data0_arranged((y*2+2*(t-2))*11+6,:) = file_data0((y*2+2*(t-2))*11+8,:);
file_data0_arranged((y*2+2*(t-2))*11+7,:) = file_data0((y*2+2*(t-2))*11+5,:);
file_data0_arranged((y*2+2*(t-2))*11+8,:) = file_data0((y*2+2*(t-2))*11+10,:);
file_data0_arranged((y*2+2*(t-2))*11+9,:) = file_data0((y*2+2*(t-2))*11+7,:);
file_data0_arranged((y*2+2*(t-2))*11+10,:) = file_data0((y*2+2*(t-2))*11+12,:);
file_data0_arranged((y*2+2*(t-2))*11+11,:) = file_data0((y*2+2*(t-2))*11+9,:);
file_data0_arranged((y*2+2*(t-2))*11+12,:) = file_data0((y*2+2*(t-2))*11+13,:);
file_data0_arranged((y*2+2*(t-2))*11+13,:) = file_data0((y*2+2*(t-2))*11+14,:);
file_data0_arranged((y*2+2*(t-2))*11+14,:) = file_data0((y*2+2*(t-2))*11+15,:);
file_data0_arranged((y*2+2*(t-2))*11+15,:) = file_data0((y*2+2*(t-2))*11+16,:);
file_data0_arranged((y*2+2*(t-2))*11+16,:) = file_data0((y*2+2*(t-2))*11+11,:);


for layer2 =1:h
   file_n(1:2*y+2*(t-2),(layer2-1)*3+1:layer2*3) = file_data0_arranged((layer2-1)*(2*y+2*(t-2))+1:layer2*(2*y+2*(t-2)),1:3);
end

fileNum = size(file_n);


%同じ層のそれぞれのfiberの質点のz軸方向の高さを平均でそろえてる
% for i=1:h
%     file_n(1:y*4-4,3*i)=mean(file_n(1:y*4-4,3*i));
% end
% for i=1:h
%     file_n(1:(y*2+t*2-4),3*i)=mean(file_n(1:(y*2+t*2-4),3*i));
% end
% 
% %3列になるように並び替え
% for i=1:fileNum(2)/3
%     mix(1+(i-1)*(y*2+t*2-4):(y*2+t*2-4)*i,1:3)=file_n(1:y*2+t*2-4,1+3*(i-1):3*i);
% end
% 
% file_name_mix = ['muscle\',muscle_name, '_min_n.csv'];
% % file_name_mix = strcat("C:\Users\bubbl\Documents\shizuya_M1\DefMuscle_for_TUS_2018\matlab\surface2grid_0608\muscle\", muscle_name, "_min_n.csv")
% csvwrite(file_name_mix, mix);
%
%
%
%% 内部点を作成
%縦線，一次の多項式(線形)作成，xy平面
% for i=1:h
%     for j=1:y-2
%         line_t(j,:,i)=polyfit([file_n(j+1,1+3*(i-1));file_n(j+1+y+(t-2)*2,1+3*(i-1))],[file_n(j+1,2+3*(i-1));file_n(j+1+y+(t-2)*2,2+3*(i-1))],1);
%     end
% end
% 
% 
% %横線
% for i=1:h
%     for j=1:t-2
%         line_y(j,:,i)=polyfit([file_n(y+2*j-1,1+3*(i-1));file_n(y+1+2*j-1,1+3*(i-1))],[file_n(y+2*j-1,2+3*(i-1));file_n(y+1+2*j-1,2+3*(i-1))],1);
%     end
% end

for i=1:h
    for j=1:y-2
        line_t(j,:,i)=polyfit([file_n((y+2*j-1),1+3*(i-1)), file_n((y+2*j),1+3*(i-1))],[file_n((y+2*j-1),2+3*(i-1)), file_n((y+2*j),2+3*(i-1))],1);
    end
end

%横線
for i=1:h
    for j=1:t-2
        line_y(j,:,i)=polyfit([file_n((y+2*j-1),2+3*(i-1)), file_n((y+2*j),2+3*(i-1))],[file_n((y+2*j-1),1+3*(i-1)), file_n((y+2*j),1+3*(i-1))],1);
    end
end
%立方体のように角度が90度となる場合
if abs(line_y(1,1,1)) > 10^5
    file_mix(1:y, :) = file_n(1:y, :);
    file_mix(y*(t-1)+1:y*t, :) = file_n(y+2*(t-2)+1:2*y+2*(t-2), :);
    for i=1:h
        for j=1:t-2
            for k=1:y
                file_mix(y*j+k,1+3*(i-1)) = file_n(y+j*2-1,1+3*(i-1)) + (file_n(y+j*2,1+3*(i-1)) - file_n(y+j*2-1,1+3*(i-1)))/(y-1)*(k-1);
                file_mix(y*j+k,2+3*(i-1)) = file_n(y+j*2,2+3*(i-1));
            end
        end
%         file_mix(y*(t-1)+1:y*t,1+3*(i-1)) = file_n(y+2*(t-2)+1:2*y+2*(t-2),1+3*(i-1));
%         file_mix(y*(t-1)+1:y*t,1+3*(i-1)) = file_n(y+2*(t-2)+1:2*y+2*(t-2),2+3*(i-1));
        file_mix(1:y*t,3*i) = file_n(1,3*i);
    end
else
    %% 占部さんver 
%     file_n(y*(t-1):y*t,:)=file_n(3*y-4:4*y-4,:);
%     file_n(y*(t-2):y*(t-2)+1,:)=file(y+2*(t-3):y+2*(t-3)+1,:);%なぞ
%     file_n(y*(t-3):y*(t-3)+1,:)=file(y+2*(t-4):y+2*(t-4)+1,:);%なぞ
%     file_n(1:4,:)=file_n(1:4,:);
%     file_n(6,:)=file_n(5,:);
%     %交点=(切片の差2-1)/(傾きの差1-2)
%     file_p(1:y*t,:)=0;
%     for i=1:h
%         for j=1:y-2
%             for k=1:t-2
%                 file_p(k+y+1+y*(j-1),1+3*(i-1))=(line_t(k,2,i)-line_y(j,2,i))/(line_y(j,1,i)-line_t(k,1,i));
%                 file_p(k+y+1+y*(j-1),2+3*(i-1))=line_y(j,1,i)*file_p(k+y+1+y*(j-1),1+3*(i-1))+line_y(j,2,i);
%                 file_p(k+y+1+y*(j-1),3+3*(i-1))=file_n(k,3+3*(i-1));
%                 file_n(k+y+1+y*(j-1),1+3*(i-1):3+3*(i-1))=0;
%             end
%         end
%     end
%     
%     %交点=(切片の差2-1)/(傾きの差1-2)
%     fileSize = size(file_n);
%     if file_n(1,fileSize(2)) == 0
%         file_n(:,fileSize(2)) = [];
%     end
%     file_p(1:y*t,:)=0;
%     file_mix=file_n+file_p;

%% 靜谷ver
    file_mix(1:y, :) = file_n(1:y, :);
    file_mix(y*(t-1)+1:y*t, :) = file_n(y+2*(t-2)+1:2*y+2*(t-2), :);
    for i=1:h
        for j=1:t-2
            xarray = linspace(file_n(y+j*2-1,1+3*(i-1)), file_n(y+j*2,1+3*(i-1)),y);%左端と右端のx座標を等間隔で
            file_mix(y*j+1,1+3*(i-1)) = file_n(y+j*2-1,1+3*(i-1));
            file_mix(y*(j+1),1+3*(i-1)) = file_n(y+j*2,1+3*(i-1));
            file_mix(y*j+1,2+3*(i-1)) = file_n(y+j*2-1,2+3*(i-1));
            file_mix(y*(j+1),2+3*(i-1)) = file_n(y+j*2,2+3*(i-1));         
            for k=1:y-2
                file_mix(y*j+k+1,1+3*(i-1)) = xarray(k+1);
                file_mix(y*j+k+1,2+3*(i-1)) = line_t(j,1,i)*file_mix(y*j+k+1,1+3*(i-1)) + line_t(j,2,i);
            end
        end
        file_mix(1:y*t,3*i) = file_n(1,3*i);
    end
    %file_nのh*3+1列目に，なぜか値が0のセルが発生するから削除
    fileSize = size(file_n);
    if file_n(1,fileSize(2)) == 0
        file_n(:,fileSize(2)) = [];
    end

 %% 靜谷ver2 5x5を想定　断面が円の場合
%     file_mix(1:y, :) = file_n(1:y, :);
%     file_mix(y*(t-1)+1:y*t, :) = file_n(y+2*(t-2)+1:2*y+2*(t-2), :);%16行のfile_nを25行のfile_mixにする
%     for i=1:h
%         xarray1 = linspace(file_n(1,1+3*(i-1)), file_n(2*y+2*(t-2),1+3*(i-1)),y);%左端と右端のx座標を等間隔で
%         xarray2 = linspace(file_n(y+t-2,1+3*(i-1)), file_n(y+t-1,1+3*(i-1)),y);%左端と右端のx座標を等間隔で
%         xarray3 = linspace(file_n(y+2*(t-2)+1,1+3*(i-1)), file_n(y,1+3*(i-1)),y);%左端と右端のx座標を等間隔で
%         xarray4 = linspace(file_n(y+2*(t-2)+int8(y/2),1+3*(i-1)), file_n(int8(y/2),1+3*(i-1)),y);%左端と右端のx座標を等間隔で
%         yarray1 = linspace(file_n(1,2+3*(i-1)), file_n(2*y+2*(t-2),2+3*(i-1)),y);%左端と右端のx座標を等間隔で
%         yarray2 = linspace(file_n(y+t-2,2+3*(i-1)), file_n(y+t-1,2+3*(i-1)),y);%左端と右端のx座標を等間隔で
%         yarray3 = linspace(file_n(y+2*(t-2)+1,2+3*(i-1)), file_n(y,2+3*(i-1)),y);%左端と右端のx座標を等間隔で
%         yarray4 = linspace(file_n(y+2*(t-2)+int8(y/2),2+3*(i-1)), file_n(int8(y/2),2+3*(i-1)),y);%左端と右端のx座標を等間隔で
%         
%         for j=1:t-2
%             file_mix(y*j+1,1+3*(i-1)) = file_n(y+j*2-1,1+3*(i-1));%両端の座標はfile_nの値をそのまま使用できる
%             file_mix(y*(j+1),1+3*(i-1)) = file_n(y+j*2,1+3*(i-1));
%             file_mix(y*j+1,2+3*(i-1)) = file_n(y+j*2-1,2+3*(i-1));
%             file_mix(y*(j+1),2+3*(i-1)) = file_n(y+j*2,2+3*(i-1));
%         end
%         
%         
%         %2行目の内部補完
%         file_mix(y+2,1+3*(i-1)) = xarray1(2);
%         file_mix(y+2,2+3*(i-1)) = yarray1(2);
%         file_mix(y+3,1+3*(i-1)) = xarray4(4);
%         file_mix(y+3,2+3*(i-1)) = yarray4(4);
%         file_mix(y+4,1+3*(i-1)) = xarray3(4);
%         file_mix(y+4,2+3*(i-1)) = yarray3(4);
%         
%         %3行目の内部補完
%         file_mix(2*y+2,1+3*(i-1)) = xarray2(2);
%         file_mix(2*y+2,2+3*(i-1)) = yarray2(2);
%         file_mix(2*y+3,1+3*(i-1)) = xarray4(3);
%         file_mix(2*y+3,2+3*(i-1)) = yarray4(3);
%         file_mix(2*y+4,1+3*(i-1)) = xarray2(4);
%         file_mix(2*y+4,2+3*(i-1)) = yarray2(4);
%         
%         %4行目の内部補完
%         file_mix(3*y+2,1+3*(i-1)) = xarray3(2);
%         file_mix(3*y+2,2+3*(i-1)) = yarray3(2);
%         file_mix(3*y+3,1+3*(i-1)) = xarray4(2);
%         file_mix(3*y+3,2+3*(i-1)) = yarray4(2);
%         file_mix(3*y+4,1+3*(i-1)) = xarray1(4);
%         file_mix(3*y+4,2+3*(i-1)) = yarray1(4);
%         
%         file_mix(1:y*t,3*i) = file_n(1,3*i);
%     end
%     %file_nのh*3+1列目に，なぜか値が0のセルが発生するから削除
%     fileSize = size(file_n);
%     if file_n(1,fileSize(2)) == 0
%         file_n(:,fileSize(2)) = [];
%     end
% %% 靜谷ver3
%     file_mix(1:y, :) = file_n(1:y, :);
%     file_mix(y*(t-1)+1:y*t, :) = file_n(y+2*(t-2)+1:2*y+2*(t-2), :);
%     for i=1:h
%         for j=1:t-2
%             xarray = linspace(file_n(y+j*2-1,1+3*(i-1)), file_n(y+j*2,1+3*(i-1)),y);%左端と右端のx座標を等間隔で
%             file_mix(y*j+1,1+3*(i-1)) = file_n(y+j*2-1,1+3*(i-1));
%             file_mix(y*(j+1),1+3*(i-1)) = file_n(y+j*2,1+3*(i-1));
%             file_mix(y*j+1,2+3*(i-1)) = file_n(y+j*2-1,2+3*(i-1));
%             file_mix(y*(j+1),2+3*(i-1)) = file_n(y+j*2,2+3*(i-1));         
%             for k=1:y-2
%                 file_mix(y*j+k+1,1+3*(i-1)) = file_n(k+1,1+3*(i-1)) + (file_n(y+(t-2)*2+k+1,1+3*(i-1)) - file_n(k+1,1+3*(i-1)))*j/(t-1);
%                 file_mix(y*j+k+1,2+3*(i-1)) = file_n(k+1,2+3*(i-1)) + (file_n(y+(t-2)*2+k+1,2+3*(i-1)) - file_n(k+1,2+3*(i-1)))*j/(t-1);
%             end
%         end
%         file_mix(1:y*t,3*i) = file_n(1,3*i);
%     end
%     %file_nのh*3+1列目に，なぜか値が0のセルが発生するから削除
%     fileSize = size(file_n);
%     if file_n(1,fileSize(2)) == 0
%         file_n(:,fileSize(2)) = [];
%     end
end

for i=1:fileNum(2)/3
    muscle(1+(i-1)*y*t:y*t*i,1:3)=file_mix(1:y*t,1+3*(i-1):3*i);
end

file_name_muscle = strcat('muscle\',"biceps_ver4", '_min_final.csv');
% file_name_muscle = strcat("C:\Users\bubbl\Documents\shizuya_M1\DefMuscle_for_TUS_2018\matlab\surface2grid_0608\muscle\", muscle_name, "_min_final.csv")
csvwrite(file_name_muscle, muscle)

%% data移動
%muscle = csvread("C:\Users\bubbl\OneDrive\ドキュメント\shizuya_M1\DefMuscle_for_TUS_2018\matlab\surface2grid_0608\muscle\min\data0.csv",0,0);
for i=1:t*y
    data_move(1,1+3*(i-1):3*i)=muscle(t*y*h-t*y+i,1:3);
end

file_name_data_move = ['muscle\data_',muscle_name, '.csv'];
% file_name_data_move = strcat("C:\Users\bubbl\Documents\shizuya_M1\DefMuscle_for_TUS_2018\matlab\surface2grid_0608\muscle\data_", muscle_name, ".csv")
csvwrite(file_name_data_move, data_move)