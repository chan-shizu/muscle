% clc
% clear
muscleNames = loadMuscleName();
muscle_name = [muscleNames{1}];
%���ӓ_!�@z���ؐ��ە����ɂ��Ă���
%�_�̔z�u��xmin ymin����X�^�[�g���Ĕ����v���
%muscleName = 'test_square'
file_name_readcsv = ['aligned\',muscle_name, '_aligned.csv'];
% file_name_readcsv = strcat('C:\Users\bubbl\Documents\shizuya_M1\DefMuscle_for_TUS_2018\matlab\surface2grid_0608\aligned\', muscle_name, '_aligned.csv')
file = csvread(file_name_readcsv,0,3);
fileNum=size(file);

divisionMuscle = readmatrix("divisionMuscle.csv");
y=divisionMuscle(1);
t=divisionMuscle(2);
h=divisionMuscle(3);

%���Ԃ��Ă���s�������Ă�?�v����
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
% csvwrite('C:\Users\�蕔�@�����q\Desktop\DefMuscle_for_TUS_2018\matlab\surface2grid_0608\muscle\brachialis_min_muscle.csv',mix)
%
%
% %muscle.csv���m�F���ĕ��ʂ�node�ԍ����Ɋ���U������
% file = csvread('C:\Users\�蕔�@�����q\Desktop\DefMuscle_for_TUS_2018\matlab\surface2grid_0608\aligned\arm_k_aligned.csv',0,3);

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

%�����w�̂��ꂼ���fiber�̎��_��z�������̍����𕽋ςł��낦�Ă�
% for i=1:h
%     file_n(1:y*4-4,3*i)=mean(file_n(1:y*4-4,3*i));
% end
for i=1:h
    file_n(1:(y*2+t*2-4),3*i)=mean(file_n(1:(y*2+t*2-4),3*i));
end

%3��ɂȂ�悤�ɕ��ёւ�
for i=1:fileNum(2)/3
    mix(1+(i-1)*(y*2+t*2-4):(y*2+t*2-4)*i,1:3)=file_n(1:y*2+t*2-4,1+3*(i-1):3*i);
end

file_name_mix = ['muscle\',muscle_name, '_min_n.csv'];
% file_name_mix = strcat("C:\Users\bubbl\Documents\shizuya_M1\DefMuscle_for_TUS_2018\matlab\surface2grid_0608\muscle\", muscle_name, "_min_n.csv")
csvwrite(file_name_mix, mix);
%
%
%
%% �����_���쐬
%�c���C�ꎟ�̑�����(���`)�쐬�Cxy����
% for i=1:h
%     for j=1:y-2
%         line_t(j,:,i)=polyfit([file_n(j+1,1+3*(i-1));file_n(j+1+y+(t-2)*2,1+3*(i-1))],[file_n(j+1,2+3*(i-1));file_n(j+1+y+(t-2)*2,2+3*(i-1))],1);
%     end
% end
% 
% 
% %����
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

%����
for i=1:h
    for j=1:t-2
        line_y(j,:,i)=polyfit([file_n((y+2*j-1),2+3*(i-1)), file_n((y+2*j),2+3*(i-1))],[file_n((y+2*j-1),1+3*(i-1)), file_n((y+2*j),1+3*(i-1))],1);
    end
end
%�����̂̂悤�Ɋp�x��90�x�ƂȂ�ꍇ
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
    %% �蕔����ver 
%     file_n(y*(t-1):y*t,:)=file_n(3*y-4:4*y-4,:);
%     file_n(y*(t-2):y*(t-2)+1,:)=file(y+2*(t-3):y+2*(t-3)+1,:);%�Ȃ�
%     file_n(y*(t-3):y*(t-3)+1,:)=file(y+2*(t-4):y+2*(t-4)+1,:);%�Ȃ�
%     file_n(1:4,:)=file_n(1:4,:);
%     file_n(6,:)=file_n(5,:);
%     %��_=(�ؕЂ̍�2-1)/(�X���̍�1-2)
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
%     %��_=(�ؕЂ̍�2-1)/(�X���̍�1-2)
%     fileSize = size(file_n);
%     if file_n(1,fileSize(2)) == 0
%         file_n(:,fileSize(2)) = [];
%     end
%     file_p(1:y*t,:)=0;
%     file_mix=file_n+file_p;

%% �ΒJver
%     file_mix(1:y, :) = file_n(1:y, :);
%     file_mix(y*(t-1)+1:y*t, :) = file_n(y+2*(t-2)+1:2*y+2*(t-2), :);
%     for i=1:h
%         for j=1:t-2
%             xarray = linspace(file_n(y+j*2-1,1+3*(i-1)), file_n(y+j*2,1+3*(i-1)),y);%���[�ƉE�[��x���W�𓙊Ԋu��
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
%     %file_n��h*3+1��ڂɁC�Ȃ����l��0�̃Z�����������邩��폜
%     fileSize = size(file_n);
%     if file_n(1,fileSize(2)) == 0
%         file_n(:,fileSize(2)) = [];
%     end

 %% �ΒJver2 5x5��z��@�f�ʂ��~�̏ꍇ
    file_mix(1:y, :) = file_n(1:y, :);
    file_mix(y*(t-1)+1:y*t, :) = file_n(y+2*(t-2)+1:2*y+2*(t-2), :);%16�s��file_n��25�s��file_mix�ɂ���
    for i=1:h
        xarray1 = linspace(file_n(1,1+3*(i-1)), file_n(2*y+2*(t-2),1+3*(i-1)),y);%���[�ƉE�[��x���W�𓙊Ԋu��
        xarray2 = linspace(file_n(y+t-2,1+3*(i-1)), file_n(y+t-1,1+3*(i-1)),y);%���[�ƉE�[��x���W�𓙊Ԋu��
        xarray3 = linspace(file_n(y+2*(t-2)+1,1+3*(i-1)), file_n(y,1+3*(i-1)),y);%���[�ƉE�[��x���W�𓙊Ԋu��
        xarray4 = linspace(file_n(y+2*(t-2)+int8(y/2),1+3*(i-1)), file_n(int8(y/2),1+3*(i-1)),y);%���[�ƉE�[��x���W�𓙊Ԋu��
        yarray1 = linspace(file_n(1,2+3*(i-1)), file_n(2*y+2*(t-2),2+3*(i-1)),y);%���[�ƉE�[��x���W�𓙊Ԋu��
        yarray2 = linspace(file_n(y+t-2,2+3*(i-1)), file_n(y+t-1,2+3*(i-1)),y);%���[�ƉE�[��x���W�𓙊Ԋu��
        yarray3 = linspace(file_n(y+2*(t-2)+1,2+3*(i-1)), file_n(y,2+3*(i-1)),y);%���[�ƉE�[��x���W�𓙊Ԋu��
        yarray4 = linspace(file_n(y+2*(t-2)+int8(y/2),2+3*(i-1)), file_n(int8(y/2),2+3*(i-1)),y);%���[�ƉE�[��x���W�𓙊Ԋu��
        
        for j=1:t-2
            file_mix(y*j+1,1+3*(i-1)) = file_n(y+j*2-1,1+3*(i-1));%���[�̍��W��file_n�̒l�����̂܂܎g�p�ł���
            file_mix(y*(j+1),1+3*(i-1)) = file_n(y+j*2,1+3*(i-1));
            file_mix(y*j+1,2+3*(i-1)) = file_n(y+j*2-1,2+3*(i-1));
            file_mix(y*(j+1),2+3*(i-1)) = file_n(y+j*2,2+3*(i-1));
        end
        
        
        %2�s�ڂ̓����⊮
        file_mix(y+2,1+3*(i-1)) = xarray1(2);
        file_mix(y+2,2+3*(i-1)) = yarray1(2);
        file_mix(y+3,1+3*(i-1)) = xarray4(4);
        file_mix(y+3,2+3*(i-1)) = yarray4(4);
        file_mix(y+4,1+3*(i-1)) = xarray3(4);
        file_mix(y+4,2+3*(i-1)) = yarray3(4);
        
        %3�s�ڂ̓����⊮
        file_mix(2*y+2,1+3*(i-1)) = xarray2(2);
        file_mix(2*y+2,2+3*(i-1)) = yarray2(2);
        file_mix(2*y+3,1+3*(i-1)) = xarray4(3);
        file_mix(2*y+3,2+3*(i-1)) = yarray4(3);
        file_mix(2*y+4,1+3*(i-1)) = xarray2(4);
        file_mix(2*y+4,2+3*(i-1)) = yarray2(4);
        
        %4�s�ڂ̓����⊮
        file_mix(3*y+2,1+3*(i-1)) = xarray3(2);
        file_mix(3*y+2,2+3*(i-1)) = yarray3(2);
        file_mix(3*y+3,1+3*(i-1)) = xarray4(2);
        file_mix(3*y+3,2+3*(i-1)) = yarray4(2);
        file_mix(3*y+4,1+3*(i-1)) = xarray1(4);
        file_mix(3*y+4,2+3*(i-1)) = yarray1(4);
        
        file_mix(1:y*t,3*i) = file_n(1,3*i);
    end
    %file_n��h*3+1��ڂɁC�Ȃ����l��0�̃Z�����������邩��폜
    fileSize = size(file_n);
    if file_n(1,fileSize(2)) == 0
        file_n(:,fileSize(2)) = [];
    end
% %% �ΒJver3
%     file_mix(1:y, :) = file_n(1:y, :);
%     file_mix(y*(t-1)+1:y*t, :) = file_n(y+2*(t-2)+1:2*y+2*(t-2), :);
%     for i=1:h
%         for j=1:t-2
%             xarray = linspace(file_n(y+j*2-1,1+3*(i-1)), file_n(y+j*2,1+3*(i-1)),y);%���[�ƉE�[��x���W�𓙊Ԋu��
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
%     %file_n��h*3+1��ڂɁC�Ȃ����l��0�̃Z�����������邩��폜
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

%% data�ړ�
%muscle = csvread("C:\Users\bubbl\OneDrive\�h�L�������g\shizuya_M1\DefMuscle_for_TUS_2018\matlab\surface2grid_0608\muscle\min\data0.csv",0,0);
for i=1:t*y
    data_move(1,1+3*(i-1):3*i)=muscle(t*y*h-t*y+i,1:3);
end

file_name_data_move = ['muscle\data_',muscle_name, '.csv'];
% file_name_data_move = strcat("C:\Users\bubbl\Documents\shizuya_M1\DefMuscle_for_TUS_2018\matlab\surface2grid_0608\muscle\data_", muscle_name, ".csv")
csvwrite(file_name_data_move, data_move)