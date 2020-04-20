clear
muscleNames = loadMuscleName();
muscle_name = [muscleNames{1}]
%���ӓ_!�@z���ؐ��ە����ɂ��Ă���
%�_�̔z�u��xmin ymin����X�^�[�g���Ĕ����v���
%muscleName = 'test_square'
file_name_readcsv = ['C:\Users\kou_0\OneDrive\�h�L�������g\����matlab\aligned\',muscle_name, '_aligned.csv']
% file_name_readcsv = strcat('C:\Users\bubbl\Documents\shizuya_M1\DefMuscle_for_TUS_2018\matlab\surface2grid_0608\aligned\', muscle_name, '_aligned.csv')
file = csvread(file_name_readcsv,0,3);
fileNum=size(file)

y = 5;
t = 5;
h = 6;

%���Ԃ��Ă���s�������Ă�?�v����
file(y*4,:)=[];
file(y*4-y+1,:)=[];
file(2*y+1,:)=[];
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

file_n(1,:)=file(5,:);
file_n(2,:)=file(4,:);
file_n(3,:)=file(3,:);
file_n(4,:)=file(2,:);
file_n(5,:)=file(1,:);
file_n(6,:)=file(14,:);
file_n(7,:)=file(6,:);
file_n(8,:)=file(15,:);
file_n(9,:)=file(7,:);
file_n(10,:)=file(16,:);
file_n(11,:)=file(8,:);
file_n(12,:)=file(13,:);
file_n(13,:)=file(12,:);
file_n(14,:)=file(11,:);
file_n(15,:)=file(10,:);
file_n(16,:)=file(9,:);

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

% file_n(1,:)=file(7,:);
% file_n(2,:)=file(6,:);
% file_n(3,:)=file(5,:);
% file_n(4,:)=file(8,:);
% file_n(5,:)=file(4,:);
% file_n(6,:)=file(3,:);
% file_n(7,:)=file(2,:);
% file_n(8,:)=file(1,:);

%�����w�̂��ꂼ���fiber�̎��_��z�������̍����𕽋ςł��낦�Ă�
for i=1:h
    file_n(1:y*4-4,3*i)=mean(file_n(1:y*4-4,3*i));
end


%3��ɂȂ�悤�ɕ��ёւ�
for i=1:fileNum(2)/3
    mix(1+(i-1)*(y*4-4):(y*4-4)*i,1:3)=file_n(1:y*4-4,1+3*(i-1):3*i);
end

file_name_mix = ['C:\Users\kou_0\OneDrive\�h�L�������g\����matlab\muscle\',muscle_name, '_min_n.csv']
% file_name_mix = strcat("C:\Users\bubbl\Documents\shizuya_M1\DefMuscle_for_TUS_2018\matlab\surface2grid_0608\muscle\", muscle_name, "_min_n.csv")
csvwrite(file_name_mix, mix)
%
%
%
%% �����_���쐬
%�c���C�ꎟ�̑�����(���`)�쐬�Cxy����
for i=1:h
    for j=1:y-2
        line_t(j,:,i)=polyfit([file_n(j+1,1+3*(i-1));file_n(j+1+y+(t-2)*2,1+3*(i-1))],[file_n(j+1,2+3*(i-1));file_n(j+1+y+(t-2)*2,2+3*(i-1))],1);
    end
end


%����
for i=1:h
    for j=1:t-2
        line_y(j,:,i)=polyfit([file_n(y+2*j-1,1+3*(i-1));file_n(y+1+2*j-1,1+3*(i-1))],[file_n(y+2*j-1,2+3*(i-1));file_n(y+1+2*j-1,2+3*(i-1))],1);
    end
end

file_n(y*(t-1):y*t,:)=file_n(3*y-4:4*y-4,:);
% file_n(y*(t-2):y*(t-2)+1,:)=file(y+2*(t-3):y+2*(t-3)+1,:);
% file_n(y*(t-3):y*(t-3)+1,:)=file(y+2*(t-4):y+2*(t-4)+1,:);
file_n(1:4,:)=file_n(1:4,:);
file_n(6,:)=file_n(5,:);


%��_=(�ؕЂ̍�2-1)/(�X���̍�1-2)
file_p(1:y*t,:)=0;
for i=1:h
    for j=1:y-2
        for k=1:t-2
            file_p(k+y+1+y*(j-1),1+3*(i-1))=(line_t(k,2,i)-line_y(j,2,i))/(line_y(j,1,i)-line_t(k,1,i));
            file_p(k+y+1+y*(j-1),2+3*(i-1))=line_y(j,1,i)*file_p(k+y+1+y*(j-1),1+3*(i-1))+line_y(j,2,i);
            file_p(k+y+1+y*(j-1),3+3*(i-1))=file_n(k,3+3*(i-1));
            file_n(k+y+1+y*(j-1),1+3*(i-1):3+3*(i-1))=0;
        end
    end
end

%file_n��h*3+1��ڂɁC�Ȃ����l��0�̃Z�����������邩��폜
file_size = size(file_n)
if file_n(1,size(2)) == 0
    file_n(:,size(2)) = []
end

file_mix=file_n+file_p;

for i=1:fileNum(2)/3
    muscle(1+(i-1)*y*t:y*t*i,1:3)=file_mix(1:y*t,1+3*(i-1):3*i);
end

file_name_muscle = ['C:\Users\kou_0\OneDrive\�h�L�������g\����matlab\muscle\',muscle_name, '_min_final.csv']
% file_name_muscle = strcat("C:\Users\bubbl\Documents\shizuya_M1\DefMuscle_for_TUS_2018\matlab\surface2grid_0608\muscle\", muscle_name, "_min_final.csv")
csvwrite(file_name_muscle, muscle)

%% data�ړ�
%muscle = csvread("C:\Users\bubbl\OneDrive\�h�L�������g\shizuya_M1\DefMuscle_for_TUS_2018\matlab\surface2grid_0608\muscle\min\data0.csv",0,0);
for i=1:t*y
    data_move(1,1+3*(i-1):3*i)=muscle(t*y*h-t*y+i,1:3);
end

file_name_data_move = ['C:\Users\kou_0\OneDrive\�h�L�������g\����matlab\muscle\data_',muscle_name, '.csv']
% file_name_data_move = strcat("C:\Users\bubbl\Documents\shizuya_M1\DefMuscle_for_TUS_2018\matlab\surface2grid_0608\muscle\data_", muscle_name, ".csv")
csvwrite(file_name_data_move, data_move)