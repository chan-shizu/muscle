% clc
% clear
muscleNames = loadMuscleName();
muscle_name = [muscleNames{1}];
%注意点!　zを筋線維方向にしておく
%点の配置はxmin yminからスタートして反時計回り
%muscleName = 'test_square'
file_name_readcsv = ['aligned\',muscle_name, '_aligned.csv'];
% file_name_readcsv = strcat('C:\Users\bubbl\Documents\shizuya_M1\DefMuscle_for_TUS_2018\matlab\surface2grid_0608\aligned\', muscle_name, '_aligned.csv')
file = csvread(file_name_readcsv,0,3);
fileNum=size(file);

divisionMuscle = readmatrix("divisionMuscle.csv");
y=divisionMuscle(1);
t=divisionMuscle(2);
h=divisionMuscle(3);

%かぶっている行を消してる?要質問
%file(y*4,:)=[];
%file(y*4-y+1,:)=[];
%file(2*y+1,:)=[];
%file(y+1,:)=[];

file(y*2+2*t,:)=[];
file(y*2+t*2-t+1,:)=[];
file(t+y+1,:)=[];
file(y+1,:)=[];


% for i=1:fileNum(2)/3
% mix(1+(i-1)*(y*4-4):(y*4-4)*i,1:3)=file(1:y*4-4,1+3*(i-1):3*i);
% end
%
% csvwrite('C:\Users\占部　麻里子\Desktop\DefMuscle_for_TUS_2018\matlab\surface2grid_0608\muscle\brachialis_min_muscle.csv',mix)
%
%
% %muscle.csvを確認して平面のnode番号順に割り振った後
% file = csvread('C:\Users\占部　麻里子\Desktop\DefMuscle_for_TUS_2018\matlab\surface2grid_0608\aligned\arm_k_aligned.csv',0,3);

% file(y*4,:)=[];
% file(y*4-y+1,:)=[];
% file(2*y+1,:)=[];
% file(y+1,:)=[];

% file_n(1,:)=file(1,:);
% file_n(2,:)=file(2,:);
% file_n(3,:)=file(3,:);
% file_n(4,:)=file(4,:);
% file_n(5,:)=file(5,:);
% file_n(6,:)=file(6,:);
% file_n(7,:)=file(7,:);
% file_n(8,:)=file(8,:);
% file_n(9,:)=file(9,:);
% file_n(10,:)=file(10,:);
% file_n(11,:)=file(11,:);
% file_n(12,:)=file(12,:);
% file_n(13,:)=file(13,:);
% file_n(14,:)=file(14,:);
% file_n(15,:)=file(15,:);
% file_n(16,:)=file(16,:);

file_n(1,:)=file(14,:);
file_n(2,:)=file(15,:);
file_n(3,:)=file(16,:);
file_n(4,:)=file(13,:);
file_n(5,:)=file(12,:);
file_n(6,:)=file(5,:);
file_n(7,:)=file(11,:);
file_n(8,:)=file(4,:);
file_n(9,:)=file(10,:);
file_n(10,:)=file(3,:);
file_n(11,:)=file(9,:);
file_n(12,:)=file(2,:);
file_n(13,:)=file(1,:);
file_n(14,:)=file(6,:);
file_n(15,:)=file(7,:);
file_n(16,:)=file(8,:);

% file_n(1,:)=file(7,:);
% file_n(2,:)=file(6,:);
% file_n(3,:)=file(5,:);
% file_n(4,:)=file(8,:);
% file_n(5,:)=file(4,:);
% file_n(6,:)=file(3,:);
% file_n(7,:)=file(2,:);
% file_n(8,:)=file(1,:);

%同じ層のそれぞれのfiberの質点のz軸方向の高さを平均でそろえてる
% for i=1:h
%     file_n(1:y*4-4,3*i)=mean(file_n(1:y*4-4,3*i));
% end
for i=1:h
    file_n(1:(y*2+t*2-4),3*i)=mean(file_n(1:(y*2+t*2-4),3*i));
end

%3列になるように並び替え
for i=1:fileNum(2)/3
    mix(1+(i-1)*(y*2+t*2-4):(y*2+t*2-4)*i,1:3)=file_n(1:y*2+t*2-4,1+3*(i-1):3*i);
end

file_name_mix = ['muscle\',muscle_name, '_min_n.csv'];
% file_name_mix = strcat("C:\Users\bubbl\Documents\shizuya_M1\DefMuscle_for_TUS_2018\matlab\surface2grid_0608\muscle\", muscle_name, "_min_n.csv")
csvwrite(file_name_mix, mix);
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
%                 file_mix(y*j+k+1,1+3*(i-1)) = xarray(k+1);
%                 file_mix(y*j+k+1,2+3*(i-1)) = line_t(j,1,i)*file_mix(y*j+k+1,1+3*(i-1)) + line_t(j,2,i);
%             end
%         end
%         file_mix(1:y*t,3*i) = file_n(1,3*i);
%     end
%     %file_nのh*3+1列目に，なぜか値が0のセルが発生するから削除
%     fileSize = size(file_n);
%     if file_n(1,fileSize(2)) == 0
%         file_n(:,fileSize(2)) = [];
%     end

 %% 靜谷ver2 5x5を想定　断面が円の場合
    file_mix(1:y, :) = file_n(1:y, :);
    file_mix(y*(t-1)+1:y*t, :) = file_n(y+2*(t-2)+1:2*y+2*(t-2), :);%16行のfile_nを25行のfile_mixにする
    for i=1:h
        xarray1 = linspace(file_n(1,1+3*(i-1)), file_n(2*y+2*(t-2),1+3*(i-1)),y);%左端と右端のx座標を等間隔で
        xarray2 = linspace(file_n(y+t-2,1+3*(i-1)), file_n(y+t-1,1+3*(i-1)),y);%左端と右端のx座標を等間隔で
        xarray3 = linspace(file_n(y+2*(t-2)+1,1+3*(i-1)), file_n(y,1+3*(i-1)),y);%左端と右端のx座標を等間隔で
        xarray4 = linspace(file_n(y+2*(t-2)+int8(y/2),1+3*(i-1)), file_n(int8(y/2),1+3*(i-1)),y);%左端と右端のx座標を等間隔で
        yarray1 = linspace(file_n(1,2+3*(i-1)), file_n(2*y+2*(t-2),2+3*(i-1)),y);%左端と右端のx座標を等間隔で
        yarray2 = linspace(file_n(y+t-2,2+3*(i-1)), file_n(y+t-1,2+3*(i-1)),y);%左端と右端のx座標を等間隔で
        yarray3 = linspace(file_n(y+2*(t-2)+1,2+3*(i-1)), file_n(y,2+3*(i-1)),y);%左端と右端のx座標を等間隔で
        yarray4 = linspace(file_n(y+2*(t-2)+int8(y/2),2+3*(i-1)), file_n(int8(y/2),2+3*(i-1)),y);%左端と右端のx座標を等間隔で
        
        for j=1:t-2
            file_mix(y*j+1,1+3*(i-1)) = file_n(y+j*2-1,1+3*(i-1));%両端の座標はfile_nの値をそのまま使用できる
            file_mix(y*(j+1),1+3*(i-1)) = file_n(y+j*2,1+3*(i-1));
            file_mix(y*j+1,2+3*(i-1)) = file_n(y+j*2-1,2+3*(i-1));
            file_mix(y*(j+1),2+3*(i-1)) = file_n(y+j*2,2+3*(i-1));
        end
        
        
        %2行目の内部補完
        file_mix(y+2,1+3*(i-1)) = xarray1(2);
        file_mix(y+2,2+3*(i-1)) = yarray1(2);
        file_mix(y+3,1+3*(i-1)) = xarray4(4);
        file_mix(y+3,2+3*(i-1)) = yarray4(4);
        file_mix(y+4,1+3*(i-1)) = xarray3(4);
        file_mix(y+4,2+3*(i-1)) = yarray3(4);
        
        %3行目の内部補完
        file_mix(2*y+2,1+3*(i-1)) = xarray2(2);
        file_mix(2*y+2,2+3*(i-1)) = yarray2(2);
        file_mix(2*y+3,1+3*(i-1)) = xarray4(3);
        file_mix(2*y+3,2+3*(i-1)) = yarray4(3);
        file_mix(2*y+4,1+3*(i-1)) = xarray2(4);
        file_mix(2*y+4,2+3*(i-1)) = yarray2(4);
        
        %4行目の内部補完
        file_mix(3*y+2,1+3*(i-1)) = xarray3(2);
        file_mix(3*y+2,2+3*(i-1)) = yarray3(2);
        file_mix(3*y+3,1+3*(i-1)) = xarray4(2);
        file_mix(3*y+3,2+3*(i-1)) = yarray4(2);
        file_mix(3*y+4,1+3*(i-1)) = xarray1(4);
        file_mix(3*y+4,2+3*(i-1)) = yarray1(4);
        
        file_mix(1:y*t,3*i) = file_n(1,3*i);
    end
    %file_nのh*3+1列目に，なぜか値が0のセルが発生するから削除
    fileSize = size(file_n);
    if file_n(1,fileSize(2)) == 0
        file_n(:,fileSize(2)) = [];
    end
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

file_name_muscle = ['muscle\',muscle_name, '_min_final.csv'];
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