% 両面が強制変位面、均等に縮む　筋張力は片方だけに
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
divisionMuscle = readmatrix("divisionMuscle.csv");
y=divisionMuscle(1);
t=divisionMuscle(2);
h=divisionMuscle(3);

powerFlag = 1;
if powerFlag == 2
    muscleFiberFInput = readmatrix("D:\Documents\matlabまとめ\研究matlab\muscleFiberPower.csv");
    %     muscleFiberFInput = zeros(300,2019);
end

mass=0.015;%5*10^(-3);
gravityG = -9.8*500*0;%重力加速度 [m/s^2]
bNum=(y-1)*(t-1)*(h-1)/8;   %blockNumber"C:\Users\bubbl\Documents\shizuya_M1\DefMuscle_for_TUS_2018\matlab\surface2grid_0608\muscle\min\data_test_square.csv"
%-9.98031度回転させる(筋肉ごとに代わる)
file_name_data_top = strcat("muscle\data_", muscle_name, "_top.csv");%"data_rot_h.csv"
file_name_data_bottom = strcat("muscle\data_", muscle_name, "_bottom.csv");%"data_rot_h.csv"
file_name_data0 = strcat("muscle\", muscle_name, "_min_final.csv");%"data0_rot_h.csv"
file_name_se = strcat("parameter\", muscle_name, "_se.csv");
file_name_tetra = strcat("parameter\", muscle_name, "_tetra.csv");
fileNameSprigk = ['parameter\',muscle_name, '_springk.csv'];
fileNameActivityLevel = ['parameter\',muscle_name, '_ActivityLevel.csv'];

data = csvread(file_name_data_top, 0, 0);%/1000;
dataBottom = csvread(file_name_data_bottom, 0, 0);%/1000;
data0=csvread(file_name_data0, 0, 0);%/1000;
se=csvread(file_name_se, 0, 0);
tetra=csvread(file_name_tetra, 0, 0);
springkVPk = readmatrix(fileNameSprigk);%各要素のばね定数と体積保存力

timeNum = size(data);%時間ステップ
ActivityLevel = readmatrix(fileNameActivityLevel);
actSize = size(ActivityLevel);
ActivityLevel(actSize(2)+1:timeNum) = ActivityLevel(end);
ActivityLevel(1:end) = 0.0;%0.2;

pointNum = size(data0);
seNum = size(se);
tetraNum = size(tetra);
dt=1*10^(-3); %10^(-3)
% dt = 2.0/timeNum(1);

% Hillモデル計算の準備
% d = 17; %各筋肉の番号，MuscleParam～っていうcsvファイルを参照
% Mscl = csvread('MuscleParam_27m.csv', 2, 1);
PCSA = str2num(muscleNames{muscleNumber,3});
Fmax = 0.7*10^6*PCSA/y/t/(h-1); % 筋線維1本当たりの最大筋力(筋肉ごとに算出)
% Fmax = 2940/y/t/(h-1);                 % 筋線維1本当たりの最大筋力(筋肉ごとに算出)
Vsh = 0.3;
Vshl = 0.23;
Vml = 1.3;
Ver = 0.5;
PEsh = 10;                  %筋肉の部位ごとの特性値
PExm = 0.4;                 %筋肉の部位ごとの特性値

tic
%プログラムの高速化のために事前割り当て
make_variable();

c = str2num(muscleNames{muscleNumber,2});
searchListK = [0.0];
searchListC = [c];
sizeSeachListK = size(searchListK);
sizeSeachListC = size(searchListC);

for searchNK=1:sizeSeachListK(2)
    ActivityLevel(1:end)=searchListK(searchNK);
    conc = searchListC(1);
    
    for i=1:timeNum(1)%iの値は時間ステップ
        i
        %プログラムの高速化のための事前割り当て
        fvA = zeros(3,tetraNum(2)+((tetraNum(1)/5-1)*tetraNum(2)));
        fvB = fvA;
        fvC = fvA;
        fvD = fvA;

        % 逆運動学の場合
        if (powerFlag == 0) || (powerFlag == 1)
            i
            % 固定面と強制変位面の座標を指定
            specify_coordinates_each_side();
            
            % 速度ベルレ法により座標の速度と位置を決定
            velocity_verlet_method();
            
        % 順運動学の場合
        else
            if i==1
            
            % シミュレーション開始時だけ座標をこちらで指定
            direct_kinematics_initial_condition();
            
            else
                
            % 速度ベルレ法により座標の速度と位置を決定。 強制変位させる面はなしで、固定面以外のすべての質点の位置が加わる力から求まる
            direct_kinematics();
            
            end
        end       
        
        %隣接する質点ごとの距離、角度を計算
        calculate_length_and_angle();

        %% 慣性力を計算
        inertial_force();
        
        %% 収縮力(筋力)
        %
        %     Fm_x(i,1:pointNum(1))=0;    Fm_y(i,1:pointNum(1))=0;    Fm_z(i,1:pointNum(1))=0;
        
        
        %% ばね力(横)
        for n=1:6%seの種類の数，縦，横，高さ，xy平面，yz平面，xz平面
            for j=1:seNum(1)
                if n==6
                    ul(1,n)=seNum(1);
                end
                if se(j,1)>n-1
                    ul(1,n)=j-1;  %upper limit
                    break
                end
            end
        end
        
        %筋張力が片方の質点にしか加わらないように変更
        for s=1:seNum(1)
            % 質点同士がつながったばねが発生する力を計算
            springF(i,s)=springkVPk(1,s)*lengthen_s(i,s);
            
            if se(s,1)==2 %高さ方向に連なった質点同士では筋線維の力を計算
                
                % Hill'sモデルより筋張力を計算
                muscle_fiber_force();
            end
        end
        
        %%質点ごとのspringforce
        %分布範囲
        for n=1:6
            % springFを各質点の各方向(x,y,z)に分配する
            distribute_spring_force();
        end
        
        %筋張力を追加
        % 逆運動学の場合
        if powerFlag == 0
            
            % 筋張力を各方向(x,y,z)に分配する
            distribute_muscle_fiber_force();
            
        elseif powerFlag == 1
            distribute_muscle_fiber_force_ver2
            
        % 順運動学の場合
        elseif powerFlag == 2
            for k=ul(2)+1:ul(3)
                
                % こちらが入力した筋張力を各方向(x,y,z)に分配する
                distribute_input_muscle_fiber_force();
            end
        end
        
        FsSum_x(i,:)=Fs_x(i,:,1)+Fs_x(i,:,3)+Fs_x(i,:,4)+Fs_x(i,:,5)+Fs_x(i,:,6)+muscleFiberF_x(i,:); %+Fs_x(i,:,2)
        FsSum_y(i,:)=Fs_y(i,:,1)+Fs_y(i,:,3)+Fs_y(i,:,4)+Fs_y(i,:,5)+Fs_y(i,:,6)+muscleFiberF_y(i,:); %+Fs_y(i,:,2
        FsSum_z(i,:)=Fs_z(i,:,1)+Fs_z(i,:,3)+Fs_z(i,:,4)+Fs_z(i,:,5)+Fs_z(i,:,6)+muscleFiberF_z(i,:); %+Fs_z(i,:,2)
        
        
        %% 拘束条件
        % 粘性力を計算
        viscous_force();
        
        %% 体積保存力
        
        for j=1:tetraNum(1)/5
            
            for k=1:tetraNum(2)
                if i==1
                    V0(k+(j-1)*tetraNum(2),1)=1/6*abs(dot(cross((data0(tetra(3+5*(j-1),k),:).'-data0(tetra(2+5*(j-1),k),:).'),(data0(tetra(4+5*(j-1),k),:).'-data0(tetra(2+5*(j-1),k),:).')),((data0(tetra(5+5*(j-1),k),:).'-data0(tetra(2+5*(j-1),k),:).'))));
                end
                
                % 四面体を構成する4つのベクトルを計算
                calculate_vec();
                       
                % 四面体の体積を計算
                V(k+(j-1)*tetraNum(2),1)=1/6*abs(dot(cross((-1)*Vec2(:,3+4*(k-1)+(j-1)*tetraNum(2)*4,i),Vec1(:,2+4*(k-1)+(j-1)*tetraNum(2)*4,i)),Vec2(:,2+4*(k-1)+(j-1)*tetraNum(2)*4,i)));
                
                %体積保存力を計算
                calculate_volune_preserving_force();
            end
        end
        
        % fvA,fvB,fvC,fvDをそれぞれの質点に方向ごと(x,y,z)に配分
        distribute_fvA_to_mass_points();
        distribute_fvB_to_mass_points();
        distribute_fvC_to_mass_points();
        distribute_fvD_to_mass_points();
        
        if i==1
            Fv_x(i,1:pointNum(1))=0;
            Fv_y(i,1:pointNum(1))=0;
            Fv_z(i,1:pointNum(1))=0;
        else
            
            Fv_x(i,1:pointNum(1)) = FvA(1:pointNum(1),1).'+FvB(1:pointNum(1),1).'+FvC(1:pointNum(1),1).'+FvD(1:pointNum(1),1).';
            Fv_y(i,1:pointNum(1)) = FvA(1:pointNum(1),2).'+FvB(1:pointNum(1),2).'+FvC(1:pointNum(1),2).'+FvD(1:pointNum(1),2).';
            Fv_z(i,1:pointNum(1)) = FvA(1:pointNum(1),3).'+FvB(1:pointNum(1),3).'+FvC(1:pointNum(1),3).'+FvD(1:pointNum(1),3).';
            
            %z軸方向の力を0にして，x,y軸方向に分散(xy平面に写像)
            if (powerFlag == 0) || (powerFlag == 1)
                FvTotal(i,1:pointNum(1)) = sqrt(Fv_x(i,1:pointNum(1)).^2 + Fv_y(i,1:pointNum(1)).^2 + Fv_z(i,1:pointNum(1)).^2);
                alpha = atan(Fv_y(i,1:pointNum(1)) ./ Fv_x(i,1:pointNum(1)));
                alpha = alpha + pi.*(Fv_x(i,1:pointNum(1)) < 0);
                Fv_x(i,1:pointNum(1)) =  FvTotal(i,1:pointNum(1)).*cos(alpha);
                Fv_y(i,1:pointNum(1)) =  FvTotal(i,1:pointNum(1)).*sin(alpha);
                Fv_x(i,1:pointNum(1)) = fillmissing(Fv_x(i,1:pointNum(1)),'constant',0);
                Fv_y(i,1:pointNum(1)) = fillmissing(Fv_y(i,1:pointNum(1)),'constant',0);
            end
            
        end
        
        %0：横　1：縦　2：高さ　3：xy平面斜め　4：xz平面斜め　5：yz平面斜め　6：空間斜め       
        
        %% 質点に加わる力と速度の算出
        %ここまで計算してきたばね, 筋張力, 体積保存力などの力を足す。その値を使用して速度を更新する。
        calculate_force_And_velocity();
        
        plot_volume();
        if volumeRatio(i) > 2.0
            break
        end
        
    end
    toc
    volumeRatio(end-1:end)=[]
    figure1 = figure()
    plot([1:(i-2)],volumeRatio);
    %         graphName = "springk=" + searchListK(searchNK) + "_c=" + searchListC(searchNC);
    graphName = "Number" + searchNK;
    saveas(figure1,graphName);
    VolumeHistroy(searchNK) = volumeRatio(end);
    volumeRatio = [];
    close
end

F_fib=springF(:,ul(1,2)+1:ul(1,2)+y*t*(h-1));

datai_x_name=strcat('output\', muscle_name, '_data_x.csv');
datai_y_name=strcat('output\', muscle_name, '_data_y.csv');
datai_z_name=strcat('output\', muscle_name, '_data_z.csv');
csvwrite(datai_x_name,datai_x);
csvwrite(datai_y_name,datai_y);
csvwrite(datai_z_name,datai_z);
csvwrite("VolumeHistory.csv",VolumeHistroy);

if (powerFlag == 0) || (powerFlag == 1)
    csvwrite("muscleFiberPower.csv",muscleFiberF);
end

% v = VideoWriter('model_3d_7.avi');
%
% open(v);
% writeVideo(v,Frame);
% close(v);
%


