% ���ʂ������ψʖʁA�ϓ��ɏk�ށ@�ؒ��͕͂Е�������
%function area_mive=2Darea():
clc;
clear;
warning('on','all');
warning;
muscleNames = loadMuscleName();
prompt = '��r�񓪋؂Ȃ�1, �I�g�K�C�㍜�؂Ȃ�2, �s�ː㍜�؂Ȃ�3, �{�񕠋ؑO���Ȃ�4, �{�񕠋،㕠�Ȃ�5�������Ă�';
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
    muscleFiberFInput = readmatrix("D:\Documents\matlab�܂Ƃ�\����matlab\muscleFiberPower.csv");
    %     muscleFiberFInput = zeros(300,2019);
end

mass=0.015;%5*10^(-3);
gravityG = -9.8*500*0;%�d�͉����x [m/s^2]
bNum=(y-1)*(t-1)*(h-1)/8;   %blockNumber"C:\Users\bubbl\Documents\shizuya_M1\DefMuscle_for_TUS_2018\matlab\surface2grid_0608\muscle\min\data_test_square.csv"
%-9.98031�x��]������(�ؓ����Ƃɑ���)
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
springkVPk = readmatrix(fileNameSprigk);%�e�v�f�̂΂˒萔�Ƒ̐ϕۑ���

timeNum = size(data);%���ԃX�e�b�v
ActivityLevel = readmatrix(fileNameActivityLevel);
actSize = size(ActivityLevel);
ActivityLevel(actSize(2)+1:timeNum) = ActivityLevel(end);
ActivityLevel(1:end) = 0.0;%0.2;

pointNum = size(data0);
seNum = size(se);
tetraNum = size(tetra);
dt=1*10^(-3); %10^(-3)
% dt = 2.0/timeNum(1);

% Hill���f���v�Z�̏���
% d = 17; %�e�ؓ��̔ԍ��CMuscleParam�`���Ă���csv�t�@�C�����Q��
% Mscl = csvread('MuscleParam_27m.csv', 2, 1);
PCSA = str2num(muscleNames{muscleNumber,3});
Fmax = 0.7*10^6*PCSA/y/t/(h-1); % �ؐ���1�{������̍ő�ؗ�(�ؓ����ƂɎZ�o)
% Fmax = 2940/y/t/(h-1);                 % �ؐ���1�{������̍ő�ؗ�(�ؓ����ƂɎZ�o)
Vsh = 0.3;
Vshl = 0.23;
Vml = 1.3;
Ver = 0.5;
PEsh = 10;                  %�ؓ��̕��ʂ��Ƃ̓����l
PExm = 0.4;                 %�ؓ��̕��ʂ��Ƃ̓����l

tic
%�v���O�����̍������̂��߂Ɏ��O���蓖��
make_variable();

c = str2num(muscleNames{muscleNumber,2});
searchListK = [0.0];
searchListC = [c];
sizeSeachListK = size(searchListK);
sizeSeachListC = size(searchListC);

for searchNK=1:sizeSeachListK(2)
    ActivityLevel(1:end)=searchListK(searchNK);
    conc = searchListC(1);
    
    for i=1:timeNum(1)%i�̒l�͎��ԃX�e�b�v
        i
        %�v���O�����̍������̂��߂̎��O���蓖��
        fvA = zeros(3,tetraNum(2)+((tetraNum(1)/5-1)*tetraNum(2)));
        fvB = fvA;
        fvC = fvA;
        fvD = fvA;

        % �t�^���w�̏ꍇ
        if (powerFlag == 0) || (powerFlag == 1)
            i
            % �Œ�ʂƋ����ψʖʂ̍��W���w��
            specify_coordinates_each_side();
            
            % ���x�x�����@�ɂ����W�̑��x�ƈʒu������
            velocity_verlet_method();
            
        % ���^���w�̏ꍇ
        else
            if i==1
            
            % �V�~�����[�V�����J�n���������W��������Ŏw��
            direct_kinematics_initial_condition();
            
            else
                
            % ���x�x�����@�ɂ����W�̑��x�ƈʒu������B �����ψʂ�����ʂ͂Ȃ��ŁA�Œ�ʈȊO�̂��ׂĂ̎��_�̈ʒu�������͂��狁�܂�
            direct_kinematics();
            
            end
        end       
        
        %�אڂ��鎿�_���Ƃ̋����A�p�x���v�Z
        calculate_length_and_angle();

        %% �����͂��v�Z
        inertial_force();
        
        %% ���k��(�ؗ�)
        %
        %     Fm_x(i,1:pointNum(1))=0;    Fm_y(i,1:pointNum(1))=0;    Fm_z(i,1:pointNum(1))=0;
        
        
        %% �΂˗�(��)
        for n=1:6%se�̎�ނ̐��C�c�C���C�����Cxy���ʁCyz���ʁCxz����
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
        
        %�ؒ��͂��Е��̎��_�ɂ��������Ȃ��悤�ɕύX
        for s=1:seNum(1)
            % ���_���m���Ȃ������΂˂���������͂��v�Z
            springF(i,s)=springkVPk(1,s)*lengthen_s(i,s);
            
            if se(s,1)==2 %���������ɘA�Ȃ������_���m�ł͋ؐ��ۂ̗͂��v�Z
                
                % Hill's���f�����ؒ��͂��v�Z
                muscle_fiber_force();
            end
        end
        
        %%���_���Ƃ�springforce
        %���z�͈�
        for n=1:6
            % springF���e���_�̊e����(x,y,z)�ɕ��z����
            distribute_spring_force();
        end
        
        %�ؒ��͂�ǉ�
        % �t�^���w�̏ꍇ
        if powerFlag == 0
            
            % �ؒ��͂��e����(x,y,z)�ɕ��z����
            distribute_muscle_fiber_force();
            
        elseif powerFlag == 1
            distribute_muscle_fiber_force_ver2
            
        % ���^���w�̏ꍇ
        elseif powerFlag == 2
            for k=ul(2)+1:ul(3)
                
                % �����炪���͂����ؒ��͂��e����(x,y,z)�ɕ��z����
                distribute_input_muscle_fiber_force();
            end
        end
        
        FsSum_x(i,:)=Fs_x(i,:,1)+Fs_x(i,:,3)+Fs_x(i,:,4)+Fs_x(i,:,5)+Fs_x(i,:,6)+muscleFiberF_x(i,:); %+Fs_x(i,:,2)
        FsSum_y(i,:)=Fs_y(i,:,1)+Fs_y(i,:,3)+Fs_y(i,:,4)+Fs_y(i,:,5)+Fs_y(i,:,6)+muscleFiberF_y(i,:); %+Fs_y(i,:,2
        FsSum_z(i,:)=Fs_z(i,:,1)+Fs_z(i,:,3)+Fs_z(i,:,4)+Fs_z(i,:,5)+Fs_z(i,:,6)+muscleFiberF_z(i,:); %+Fs_z(i,:,2)
        
        
        %% �S������
        % �S���͂��v�Z
        viscous_force();
        
        %% �̐ϕۑ���
        
        for j=1:tetraNum(1)/5
            
            for k=1:tetraNum(2)
                if i==1
                    V0(k+(j-1)*tetraNum(2),1)=1/6*abs(dot(cross((data0(tetra(3+5*(j-1),k),:).'-data0(tetra(2+5*(j-1),k),:).'),(data0(tetra(4+5*(j-1),k),:).'-data0(tetra(2+5*(j-1),k),:).')),((data0(tetra(5+5*(j-1),k),:).'-data0(tetra(2+5*(j-1),k),:).'))));
                end
                
                % �l�ʑ̂��\������4�̃x�N�g�����v�Z
                calculate_vec();
                       
                % �l�ʑ̂̑̐ς��v�Z
                V(k+(j-1)*tetraNum(2),1)=1/6*abs(dot(cross((-1)*Vec2(:,3+4*(k-1)+(j-1)*tetraNum(2)*4,i),Vec1(:,2+4*(k-1)+(j-1)*tetraNum(2)*4,i)),Vec2(:,2+4*(k-1)+(j-1)*tetraNum(2)*4,i)));
                
                %�̐ϕۑ��͂��v�Z
                calculate_volune_preserving_force();
            end
        end
        
        % fvA,fvB,fvC,fvD�����ꂼ��̎��_�ɕ�������(x,y,z)�ɔz��
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
            
            %z�������̗͂�0�ɂ��āCx,y�������ɕ��U(xy���ʂɎʑ�)
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
        
        %0�F���@1�F�c�@2�F�����@3�Fxy���ʎ΂߁@4�Fxz���ʎ΂߁@5�Fyz���ʎ΂߁@6�F��Ԏ΂�       
        
        %% ���_�ɉ����͂Ƒ��x�̎Z�o
        %�����܂Ōv�Z���Ă����΂�, �ؒ���, �̐ϕۑ��͂Ȃǂ̗͂𑫂��B���̒l���g�p���đ��x���X�V����B
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


