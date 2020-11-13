%function area_mive=2Darea():
clc;
clear;
warning('on','all');
warning;
muscleNames = loadMuscleName();
muscle_name = [muscleNames{1}];
divisionMuscle = readmatrix("divisionMuscle.csv");
y=divisionMuscle(1);
t=divisionMuscle(2);
h=divisionMuscle(3);

%����
% a=csvread('a_r.csv', 0, 0);
% springk=a(1,1)/10000;    %1900
% springk = 8983400/10000;
% springk=1000;  %100
% kv = 0.10731*50000;
% kv=a(2,1)*50000;   %85

kv=0.25;    %0.25
scale = 1;
mass=0.015;%5*10^(-3);

gravityG = -9.8*500*0;%�d�͉����x [m/s^2]
conk=60;%60
conc=0;%60
bNum=(y-1)*(t-1)*(h-1)/8;   %blockNumber"C:\Users\bubbl\Documents\shizuya_M1\DefMuscle_for_TUS_2018\matlab\surface2grid_0608\muscle\min\data_test_square.csv"
%-9.98031�x��]������(�ؓ����Ƃɑ���)
file_name_data = strcat("muscle\data_", muscle_name, ".csv");%"data_rot_h.csv"
file_name_data0 = strcat("muscle\", muscle_name, "_min_final.csv");%"data0_rot_h.csv"
% file_name_data = strcat("C:\Users\bubbl\Documents\shizuya_M1\DefMuscle_for_TUS_2018\matlab\surface2grid_0608\muscle\min\urabe_arm_data.csv")%"data_rot_h.csv"
% file_name_data0 = strcat("C:\Users\bubbl\Documents\shizuya_M1\DefMuscle_for_TUS_2018\matlab\surface2grid_0608\muscle\min\urabe_arm_data0.csv")%"data0_rot_h.csv"
file_name_se = strcat("parameter\", muscle_name, "_se.csv");
file_name_tetra = strcat("parameter\", muscle_name, "_tetra.csv");
fileNameSprigk = ['parameter\',muscle_name, '_springk.csv'];
fileNameActivityLevel = ['parameter\',muscle_name, '_ActivityLevel.csv'];

data= csvread(file_name_data, 0, 0)*scale;%/1000;
data0=csvread(file_name_data0, 0, 0)*scale;%/1000;
se=csvread(file_name_se, 0, 0);
tetra=csvread(file_name_tetra, 0, 0);
springkVPk = readmatrix(fileNameSprigk);%�e�v�f�̂΂˒萔�Ƒ̐ϕۑ���
% springkVPk(1,:)=6.3343*(10^7)
% springkVPk(2,:)=1.1168*(10^-5)

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
Fmax = 2940/y/t/(h-1);                 % �ؐ���1�{������̍ő�ؗ�(�ؓ����ƂɎZ�o)
Vsh = 0.3;
Vshl = 0.23;
Vml = 1.3;
Ver = 0.5;
alpha = 0.1;         %�؂̍ő劈�������x����\���ϐ�
PEsh = 10;                  %�ؓ��̕��ʂ��Ƃ̓����l
PExm = 0.4;                 %�ؓ��̕��ʂ��Ƃ̓����l

tic
%�v���O�����̍������̂��߂Ɏ��O���蓖��
datai_x = zeros(timeNum(1),y*t*h);
datai_y = datai_x;
datai_z = datai_x;
con_x = datai_x;
con_y = datai_x;
con_z = datai_x;
Fm_xFm_y = datai_x;
Fm_z = datai_x;
Fn_x = datai_x;
Fn_y = datai_x;
Fn_z = datai_x;
Fsl_x = datai_x;
Fsl_y = datai_x;
Fsl_z = datai_x;
Fsr_x = datai_x;
Fsr_y = datai_x;
Fsr_z = datai_x;
FsSum_x = datai_x;
FsSum_y = datai_x;
FsSum_z = datai_x;
Fv_x = datai_x;
Fv_y = datai_x;
Fv_z = datai_x;
vis_x = datai_x;
vis_y = datai_x;
vis_z = datai_x;
vn_x = datai_x;
vn_y = datai_x;
vn_z = datai_x;
Mr_x = datai_x;
Mr_y = datai_x;
Mr_z = datai_x;
FCorrectx = datai_x;
FCorrecty = datai_x;
F_fib = zeros(timeNum(1),y*t*(h-1));
f_Lce = F_fib;
% f_Vce = F_fib;
HillActive = F_fib;
HillPassive = F_fib;
length_per = F_fib;
Vce = F_fib;
Fs_x = zeros(timeNum(1),y*t*h,6);
Fs_y = Fs_x;
Fs_z = Fs_x;
length_s = zeros(timeNum(1),seNum(1));
lengthen_s = length_s;
springF = length_s;

searchListK = [1]
searchListC = [0.001]
sizeSeachListK = size(searchListK);
sizeSeachListC = size(searchListC);

for searchNK=1:sizeSeachListK(2)
%     springkVPk(1,:) = searchListK(searchNK)
    for searchNC=1:sizeSeachListC(2)
        ActivityLevel(1:end)=searchListC(searchNC);
        conc = 50;
        conk = 0;
        
        for i=1:timeNum(1)%i�̒l�͎��ԃX�e�b�v
            i
            %�v���O�����̍������̂��߂̎��O���蓖��
            fvA = zeros(3,tetraNum(2)+((tetraNum(1)/5-1)*tetraNum(2)));
            fvB = fvA;
            fvC = fvA;
            fvD = fvA;
            %�������W
            if i==1
                datai_x(i,1:pointNum)=data0(:,1).';  %i���ɂ�����e�_�̍��W
                datai_y(i,1:pointNum)=data0(:,2).';
                datai_z(i,1:pointNum)=data0(:,3).';
                
            else
                %�Œ��
                datai_x(i,1:y*t)=data0(1:y*t,1);
                datai_y(i,1:y*t)=data0(1:y*t,2);
                datai_z(i,1:y*t)=data0(1:y*t,3);
                
                %�����ψʂ�^�����
                for k=1:y*t
                    datai_x(i,(pointNum(1)-k+1))=data(i,3*(y*t-k)+1);  %x���W
                    datai_y(i,(pointNum(1)-k+1))=data(i,3*(y*t-k)+2);
                    datai_z(i,(pointNum(1)-k+1))=data(i,3*(y*t-k)+3);
                    %             vn_x(i,(pointNum(1)-k+1))=vn_x(i-1,(pointNum(1)-k+1))+Fn_x(i-1,(pointNum(1)-k+1))*dt/2;
                    %             vn_y(i,(pointNum(1)-k+1))=vn_y(i-1,(pointNum(1)-k+1))+Fn_y(i-1,(pointNum(1)-k+1))*dt/2;
                    %             vn_z(i,(pointNum(1)-k+1))=vn_z(i-1,(pointNum(1)-k+1))+Fn_z(i-1,(pointNum(1)-k+1))*dt/2;
                    %             datai_x(i,(pointNum(1)-k+1))=datai_x(i-1,(pointNum(1)-k+1))+dt*vn_x(i,(pointNum(1)-k+1));
                    %             datai_y(i,(pointNum(1)-k+1))=datai_y(i-1,(pointNum(1)-k+1))+dt*vn_y(i,(pointNum(1)-k+1));
                    %             datai_z(i,(pointNum(1)-k+1))=datai_z(i-1,(pointNum(1)-k+1))+dt*vn_z(i,(pointNum(1)-k+1));
                end
            end
            
            if i==1
                vn_x(i,1:pointNum(1))=0;   %1���:node4 2���:node5
                vn_y(i,1:pointNum(1))=0;
                vn_z(i,1:pointNum(1))=0;
                
            else
                %�Œ�ʂƋ����ψʖʈȊO�𕽓��ɕψʂ�����(���x�x�����@) %���������U�̌���
                vn_x(i,(y*t+1):(pointNum(1)-y*t))=vn_x(i-1,(y*t+1):(pointNum(1)-y*t))+Fn_x(i-1,(y*t+1):(pointNum(1)-y*t))*dt/2;
                vn_y(i,(y*t+1):(pointNum(1)-y*t))=vn_y(i-1,(y*t+1):(pointNum(1)-y*t))+Fn_y(i-1,(y*t+1):(pointNum(1)-y*t))*dt/2;
                vn_z(i,(y*t+1):(pointNum(1)-y*t))=vn_z(i-1,(y*t+1):(pointNum(1)-y*t))+Fn_z(i-1,(y*t+1):(pointNum(1)-y*t))*dt/2;
                
                datai_x(i,(y*t+1):(pointNum(1)-y*t))=datai_x(i-1,(y*t+1):(pointNum(1)-y*t))+dt*vn_x(i,(y*t+1):(pointNum(1)-y*t));
                datai_y(i,(y*t+1):(pointNum(1)-y*t))=datai_y(i-1,(y*t+1):(pointNum(1)-y*t))+dt*vn_y(i,(y*t+1):(pointNum(1)-y*t));
                %         datai_z(i,(y*t+1):(pointNum(1)-y*t))=datai_z(i-1,(y*t+1):(pointNum(1)-y*t))+dt*vn_z(i,(y*t+1):(pointNum(1)-y*t));
                for n=1:(h-1)
                    if data0(1,3)<data0(y*t+1,3)
                        datai_z(i,n*y*t+1:(n+1)*y*t) = data0(1,3) + n/(h-1)*(data(i,3)-data0(1,3));
                    else
                        datai_z(i,n*y*t+1:(n+1)*y*t) = data0(1,3) - n/(h-1)*(data0(1,3)-data(i,3));
                    end
                end
                
                
                % ��(�������Ccsv��defmuscle��x������y�����𓝈ꂵ�Ă��Ȃ�����C�ԈႦ���N���肻��)
                %         for j=1:h-1
                %             for k=1:y
                %                 if datai_y(i,y*(t-1)+k+j*y*t)>data0(y*(t-1)+k+j*y*t,2)
                %                     datai_y(i,y*(t-1)+k+j*y*t)=data0(y*(t-1)+k+j*y*t,2);
                %                 end
                %             end
                %         end
                %                 %���ɕt���Ɖ���
                %�ȉ��ΒJ�ǉ�
                %         for j=1:h-1
                %             for k=1:t
                %                 datai_x(i,t*k+j*y*t)=data0(t*k+j*y*t,1);
                %             end
                %         end
                %
                %         for j=1:h-1
                %             for k=1:t
                %                 datai_y(i,20+k+j*y*t)=data0(20+k+j*y*t,2);
                %             end
                %         end
                %
                %         for j=1:h-1
                %             datai_x(i,6+j*y*t)=data0(6+j*y*t,1);
                %             datai_x(i,11+j*y*t)=data0(11+j*y*t,1);
                %             datai_x(i,16+j*y*t)=data0(16+j*y*t,1);
                %             datai_x(i,15+j*y*t)=data0(15+j*y*t,1);
                %             datai_x(i,20+j*y*t)=data0(20+j*y*t,1);
                %             datai_x(i,25+j*y*t)=data0(25+j*y*t,1);
                %         end
                %4�w�ڂ܂ō�
                % for j=1:h-3
                %         for k=1:y
                %             if datai_y(i,y*(t-1)+k+j*y*t)>data0(y*(t-1)+k+j*y*t,2)
                %         datai_y(i,y*(t-1)+k+j*y*t)=data0(y*(t-1)+k+j*y*t,2);
                %             end
                %         end
                % end
                
                % �㉺���E�̎��_�Ƃ̕��ϒl�����W�Ƃ���
                %         fivePointAverageDatax = datai_x(i,:);
                %         fivePointAverageDatay = datai_y(i,:);
                %         %         fivePointAverageDataz = datai_z(i,:);
                %         for j=1:h-2
                %             for  n =1:y-2
                %                 %                 fivePointAverageDatax(1+n+j*y*t) = (distancex(1+n+j*y*t,2)*datai_x(i,1+n+(j-1)*y*t) +  distancex(1+n+j*y*t,1)*datai_x(i,1+n+(j+1)*y*t))/distancex(1+n+j*y*t,3);
                %                 %                 fivePointAverageDatay(1+n+j*y*t) = (distancey(1+n+j*y*t,2)*datai_y(i,1+n+(j-1)*y*t) +  distancey(1+n+j*y*t,1)*datai_y(i,1+n+(j+1)*y*t))/distancey(1+n+j*y*t,3);
                %                 if (j==1 | j==h-2)
                %                     fivePointAverageDatax(1+n+j*y*t) = (datai_x(i,1+n+(j-1)*y*t) +  datai_x(i,1+n+(j+1)*y*t))/2;
                %                     fivePointAverageDatay(1+n+j*y*t) = (datai_y(i,1+n+(j-1)*y*t) +  datai_y(i,1+n+(j+1)*y*t))/2;
                %                 else
                %                     fivePointAverageDatax(1+n+j*y*t) = (datai_x(i,1+n+(j-2)*y*t) + datai_x(i,1+n+(j-1)*y*t) +  datai_x(i,1+n+(j+1)*y*t) + datai_x(i,1+n+(j+2)*y*t))/4;
                %                     fivePointAverageDatay(1+n+j*y*t) = (datai_y(i,1+n+(j-2)*y*t) + datai_y(i,1+n+(j-1)*y*t) +  datai_y(i,1+n+(j+1)*y*t) + datai_y(i,1+n+(j+2)*y*t))/4;
                %                 end
                %             end
                %         end
                %
                %         for j=1:h-2
                %             for  n =1:y-2
                %                 %                 fivePointAverageDatax(1+n+(t-1)*y+j*y*t) = (distancex(1+n+(t-1)*y+j*y*t,2)*datai_x(i,1+n+(t-1)*y+(j-1)*y*t) + distancex(1+n+(t-1)*y+j*y*t,1)*datai_x(i,1+n+(t-1)*y+(j+1)*y*t))/distancex(1+n+(t-1)*y+(j-1)*y*t,3);
                %                 %                 fivePointAverageDatay(1+n+(t-1)*y+j*y*t) = (distancey(1+n+(t-1)*y+j*y*t,2)*datai_y(i,1+n+(t-1)*y+(j-1)*y*t) + distancey(1+n+(t-1)*y+j*y*t,1)*datai_y(i,1+n+(t-1)*y+(j+1)*y*t))/distancey(1+n+(t-1)*y+(j-1)*y*t,3);
                %                 if  (j==1 | j==h-2)
                %                     fivePointAverageDatax(1+n+(t-1)*y+j*y*t) = (datai_x(i,1+n+(t-1)*y+(j-1)*y*t) + datai_x(i,1+n+(t-1)*y+(j+1)*y*t))/2;
                %                     fivePointAverageDatay(1+n+(t-1)*y+j*y*t) = (datai_y(i,1+n+(t-1)*y+(j-1)*y*t) + datai_y(i,1+n+(t-1)*y+(j+1)*y*t))/2;
                %                 else
                %                     fivePointAverageDatax(1+n+(t-1)*y+j*y*t) = (datai_x(i,1+n+(t-1)*y+(j-2)*y*t) + datai_x(i,1+n+(t-1)*y+(j-1)*y*t) + datai_x(i,1+n+(t-1)*y+(j+1)*y*t) + datai_x(i,1+n+(t-1)*y+(j+2)*y*t))/4;
                %                     fivePointAverageDatay(1+n+(t-1)*y+j*y*t) = (datai_y(i,1+n+(t-1)*y+(j-2)*y*t) + datai_y(i,1+n+(t-1)*y+(j-1)*y*t) + datai_y(i,1+n+(t-1)*y+(j+1)*y*t) + datai_y(i,1+n+(t-1)*y+(j+2)*y*t))/4;
                %                 end
                %             end
                %         end
                %
                %         for n=1:h-2
                %             for j =1:t-2
                %                 %                 fivePointAverageDatax(1+j*y+n*y*t) = (distancex(1+j*y+n*y*t,2)*datai_x(i,1+j*y+(n-1)*y*t) +  distancex(1+j*y+n*y*t,1)*datai_x(i,1+j*y+(n+1)*y*t))/distancex(1+j*y+n*y*t,3);
                %                 %                 fivePointAverageDatay(1+j*y+n*y*t) = (distancey(1+j*y+n*y*t,2)*datai_y(i,1+j*y+(n-1)*y*t) +  distancey(1+j*y+n*y*t,1)*datai_y(i,1+j*y+(n+1)*y*t))/distancey(1+j*y+n*y*t,3);
                %                 if (n==1 | n==h-2)
                %                     fivePointAverageDatax(1+j*y+n*y*t) = (datai_x(i,1+j*y+(n-1)*y*t) + datai_x(i,1+j*y+(n+1)*y*t))/2;
                %                     fivePointAverageDatay(1+j*y+n*y*t) = (datai_y(i,1+j*y+(n-1)*y*t) + datai_y(i,1+j*y+(n+1)*y*t))/2;
                %                 else
                %                     fivePointAverageDatax(1+j*y+n*y*t) = (datai_x(i,1+j*y+(n-2)*y*t) + datai_x(i,1+j*y+(n-1)*y*t) + datai_x(i,1+j*y+(n+1)*y*t) + datai_x(i,1+j*y+(n+2)*y*t))/4;
                %                     fivePointAverageDatay(1+j*y+n*y*t) = (datai_y(i,1+j*y+(n-2)*y*t) + datai_y(i,1+j*y+(n-1)*y*t) + datai_y(i,1+j*y+(n+1)*y*t) + datai_y(i,1+j*y+(n+2)*y*t))/4;
                %                 end
                %             end
                %         end
                %
                %         for n=1:h-2
                %             for j =1:t-2
                %                 %                 fivePointAverageDatax(y+j*y+n*y*t) = (distancex(y+j*y+n*y*t,2)*datai_x(i,y+j*y+(n-1)*y*t) + distancex(y+j*y+n*y*t,1)*datai_x(i,y+j*y+(n+1)*y*t))/distancex(y+j*y+n*y*t,3);
                %                 %                 fivePointAverageDatay(y+j*y+n*y*t) = (distancey(y+j*y+n*y*t,2)*datai_y(i,y+j*y+(n-1)*y*t) + distancey(y+j*y+n*y*t,1)*datai_y(i,y+j*y+(n+1)*y*t))/distancey(y+j*y+n*y*t,3);
                %                 if(n==1 | n==h-2)
                %                     fivePointAverageDatax(y+j*y+n*y*t) = (datai_x(i,y+j*y+(n-1)*y*t) + datai_x(i,y+j*y+(n+1)*y*t))/2;
                %                     fivePointAverageDatay(y+j*y+n*y*t) = (datai_y(i,y+j*y+(n-1)*y*t) + datai_y(i,y+j*y+(n+1)*y*t))/2;
                %                 else
                %                     fivePointAverageDatax(y+j*y+n*y*t) = (datai_x(i,y+j*y+(n-2)*y*t) + datai_x(i,y+j*y+(n-1)*y*t) + datai_x(i,y+j*y+(n+1)*y*t) + datai_x(i,y+j*y+(n+2)*y*t))/4;
                %                     fivePointAverageDatay(y+j*y+n*y*t) = (datai_y(i,y+j*y+(n-2)*y*t) + datai_y(i,y+j*y+(n-1)*y*t) + datai_y(i,y+j*y+(n+1)*y*t) + datai_y(i,y+j*y+(n+2)*y*t))/4;
                %                 end
                %             end
                %         end
                %         datai_x(i,:) = fivePointAverageDatax;
                %         datai_y(i,:) = fivePointAverageDatay;
                %         %         datai_z(i,:) = fivePointAverageDataz;
            end
            
            %�΂˗v�f����,se�őΉ��t����ꂽ���ׂĂ̎��_�̃L�������v�Z
            for s=1:seNum(1)
                length_s(i,s)=sqrt((datai_x(i,se(s,3))-datai_x(i,se(s,2)))^2+(datai_y(i,se(s,3))-datai_y(i,se(s,2)))^2+(datai_z(i,se(s,3))-datai_z(i,se(s,2)))^2);
                lengthen_s(i,s)=length_s(i,s)-length_s(1,s);%���R������̒���
                
                %�p�x�v�Z
                if datai_y(i,se(s,3))-datai_y(i,se(s,2))==0&&datai_x(i,se(s,3))-datai_x(i,se(s,2))>=0
                    se(s,5)=0;
                elseif datai_y(i,se(s,3))-datai_y(i,se(s,2))==0&&datai_x(i,se(s,3))-datai_x(i,se(s,2))<=0
                    se(s,5)=pi;
                elseif datai_x(i,se(s,3))-datai_x(i,se(s,2))==0&&datai_y(i,se(s,3))-datai_y(i,se(s,2))>=0
                    se(s,5)=pi/2;
                elseif datai_x(i,se(s,3))-datai_x(i,se(s,2))==0&&datai_y(i,se(s,3))-datai_y(i,se(s,2))<=0
                    se(s,5)=-pi/2;
                elseif datai_x(i,se(s,3))-datai_x(i,se(s,2))<=0&&datai_y(i,se(s,3))-datai_y(i,se(s,2))>=0
                    se(s,5)=atan((datai_y(i,se(s,3))-datai_y(i,se(s,2)))/(datai_x(i,se(s,3))-datai_x(i,se(s,2))))+pi;   %�p�x��(xy����)
                elseif datai_x(i,se(s,3))-datai_x(i,se(s,2))<=0&&datai_y(i,se(s,3))-datai_y(i,se(s,2))<=0
                    se(s,5)=atan((datai_y(i,se(s,3))-datai_y(i,se(s,2)))/(datai_x(i,se(s,3))-datai_x(i,se(s,2))))+pi;
                else
                    se(s,5)=atan((datai_y(i,se(s,3))-datai_y(i,se(s,2)))/(datai_x(i,se(s,3))-datai_x(i,se(s,2))));
                end
                
                se(s,4)=atan((datai_z(i,se(s,3))-datai_z(i,se(s,2)))/sqrt((datai_x(i,se(s,3))-datai_x(i,se(s,2)))^2+(datai_y(i,se(s,3))-datai_y(i,se(s,2)))^2));   %�p�x��(xz����
            end
            
            
            %% ������
            if i==1
                Mr_x(i,1:pointNum(1))=0;
                Mr_y(i,1:pointNum(1))=0;
                Mr_z(i,1:pointNum(1))=0;
            else
                Mr_x(i,1:pointNum(1))=mass*(datai_x(i,1:pointNum(1))-datai_x(i-1,1:pointNum(1)))/dt/dt;
                Mr_y(i,1:pointNum(1))=mass*(datai_y(i,1:pointNum(1))-datai_y(i-1,1:pointNum(1)))/dt/dt;
                Mr_z(i,1:pointNum(1))=mass*(datai_z(i,1:pointNum(1))-datai_z(i-1,1:pointNum(1)))/dt/dt;
            end
            
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
            
            for s=1:seNum(1)
                springF(i,s)=springkVPk(1,s)*lengthen_s(i,s);%�΂˂��Ƃɂ΂˒萔��ݒ�
                
                if se(s,1)==2%2  %�ؐ��ۂ���Hill�̃��f���v�Z
                    
                    Lcesh(1,s-ul(1,2)) = 0.5 * length_s(1,s);
                    f_Lce(i,s-ul(1,2)) = exp(-( lengthen_s(i,s)/Lcesh(1,s-ul(1,2))).^2);%�蕔����̏C�_�̎�(3.31)
                    if i == 1
                        Vce(i,s-ul(1,2)) = 0;
                    else
                        %                 Vce(i,s-ul(1,2)) = -(length_s(i,s) - length_s(i-1,s))/0.1;%0.1�ł͂Ȃ�dt�ł�?
                        Vce(i,s-ul(1,2)) = -(length_s(i,s) - length_s(i-1,s))/dt;
                    end
                    Vvm = 6 * length_s(1,s);%V_vm�͍ő哙�ڐ����k���̍ő呬�x��V_vm=6?l_ce
                    Vmax(1,s-ul(1,2)) = Vvm .* ( 1 - Ver .* (1 - ActivityLevel(i) * f_Lce(i,s-ul(1,2))));
                    
                    %�蕔����̏C�_��(3-32)
                    if Vce(i,s-ul(1,2)) <= -Vmax(1,s-ul(1,2))
                        f_Vce(i,s-ul(1,2)) = 0;
                    elseif (Vce(i,s-ul(1,2)) > -Vmax(1,s-ul(1,2)))  &&  (Vce(i,s-ul(1,2)) < 0)
                        f_Vce(i,s-ul(1,2)) = (Vsh * Vmax(1,s-ul(1,2)) + Vsh * Vce(i,s-ul(1,2))) / (Vsh *Vmax(1,s-ul(1,2)) - Vce(i,s-ul(1,2)));
                    elseif Vce(i,s-ul(1,2)) >= 0
                        f_Vce(i,s-ul(1,2)) = (Vsh * Vshl * Vmax(1,s-ul(1,2)) + Vml * Vce(i,s-ul(1,2))) / (Vsh * Vshl * Vmax(1,s-ul(1,2)) + Vce(i,s-ul(1,2)));%vmax�����ĂȂ�?
                        %             else
                        %                 f_Vce(i,s-ul(1,2)) = 0
                    end
                    
                    length_per(i,s-ul(1,2)) =lengthen_s(i,s)/ length_s(1,s);
                    if lengthen_s(i,s)<=0
                        HillPassive(i,s-ul(1,2)) =0;
                    else
                        HillPassive(i,s-ul(1,2)) = Fmax/exp(PEsh)* (exp(lengthen_s(i,s)/ length_s(1,s)* PEsh/PExm) -1); %
                    end
                    HillActive(i,s-ul(1,2)) = ActivityLevel(i) .* f_Lce(i,s-ul(1,2)) .* f_Vce(i,s-ul(1,2)) .* Fmax'; %�蕔����̏C�_�C��(3-29)
                    if i==1
                        springF(i,s)=0;
                    else
%                         springF(i,s)=HillActive(i,s-ul(1,2))+HillPassive(i,s-ul(1,2));%*scale;%�󓮓I�ȗ͂Ɣ\���I�ȗ̘͂a
                        springF(i,s)=-(HillActive(i,s-ul(1,2))+HillPassive(i,s-ul(1,2)));%*scale;%�󓮓I�ȗ͂Ɣ\���I�ȗ̘͂a
                    end
                end
            end
            
            %%���_���Ƃ�springforce
            %���z�͈�
            for n=1:6
                Fs_x(i,1:pointNum(1),n)=0;
                Fs_y(i,1:pointNum(1),n)=0;
                Fs_z(i,1:pointNum(1),n)=0;
                senum=(1:seNum);
                
                if n==1
                    Fsl_x(i,1:pointNum(1))=0;
                    Fsl_y(i,1:pointNum(1))=0;
                    Fsl_z(i,1:pointNum(1))=0;
                    Fsr_x(i,1:pointNum(1))=0;
                    Fsr_y(i,1:pointNum(1))=0;
                    Fsr_z(i,1:pointNum(1))=0;
                    
                    %����������
                    for k=1:ul(1,n)
                        for j=1:ul(1,n)
                            %             if n ==1
                            %                 ul_number = 1;
                            %             else
                            %                 ul_number = n-1;
                            %             end
                            %             for k=ul_number:ul(1,n)
                            %                 for j=ul_number:ul(1,n)
                            match_l(j,1)=se(j,2)==se(k,2);%(false(0)��true(1)�̔z��)
                        end
                        match_lNum=senum(match_l);   %�������_��L����O�pNo�𒊏o
                        sN=size(match_lNum);
                        if sN(2)>1
                            matchse_l=[];
                            for s=1:sN(2)   %matchnum�ɓ����Ă���O�pNo��xy�����͂�matchf_A�ɕ\��
                                matchse_l(1,s)=springF(i,match_lNum(1,s))*cos(se(match_lNum(1,s),4))*cos(se(match_lNum(1,s),5));  %x
                                matchse_l(2,s)=springF(i,match_lNum(1,s))*cos(se(match_lNum(1,s),4))*sin(se(match_lNum(1,s),5));   %y
                                matchse_l(3,s)=springF(i,match_lNum(1,s))*sin(se(match_lNum(1,s),4));   %z
                            end
                            Fsl_x(i,se(k,2))=sum(matchse_l(1,:));
                            Fsl_y(i,se(k,2))=sum(matchse_l(2,:));
                            Fsl_z(i,se(k,2))=sum(matchse_l(3,:));
                            
                        else
                            Fsl_x(i,se(k,2))=springF(i,k)*cos(se(k,4))*cos(se(k,5));
                            Fsl_y(i,se(k,2))=springF(i,k)*cos(se(k,4))*sin(se(k,5));
                            Fsl_z(i,se(k,2))=springF(i,k)*sin(se(k,4));
                        end
                        match_l(:,1)=0;
                    end
                    
                    %�E��������
                    for k=1:ul(1,n)
                        for j=1:ul(1,n)
                            %             if n ==1
                            %                 ul_number = 1;
                            %             else
                            %                 ul_number = n-1;
                            %             end
                            %             for k=ul_number:ul(1,n)
                            %                 for j=ul_number:ul(1,n)
                            match_r(j,1)=se(j,3)==se(k,3);
                        end
                        match_rNum=senum(match_r);   %�������_��L����O�pNo�𒊏o
                        sN=size(match_rNum);
                        if sN(2)>1
                            matchse_r=[];
                            for s=1:sN(2)   %matchnum�ɓ����Ă���O�pNo��xy�����͂�matchf_A�ɕ\��
                                matchse_r(1,s)=-springF(i,match_rNum(1,s))*cos(se(match_rNum(1,s),4))*cos(se(match_rNum(1,s),5));  %x
                                matchse_r(2,s)=-springF(i,match_rNum(1,s))*cos(se(match_rNum(1,s),4))*sin(se(match_rNum(1,s),5));   %y
                                matchse_r(3,s)=-springF(i,match_rNum(1,s))*sin(se(match_rNum(1,s),4));   %z
                            end
                            Fsr_x(i,se(k,3))=sum(matchse_r(1,:));
                            Fsr_y(i,se(k,3))=sum(matchse_r(2,:));
                            Fsr_z(i,se(k,3))=sum(matchse_r(3,:));
                            
                        else
                            Fsr_x(i,se(k,3))=-springF(i,k)*cos(se(k,4))*cos(se(k,5));
                            Fsr_y(i,se(k,3))=-springF(i,k)*cos(se(k,4))*sin(se(k,5));
                            Fsr_z(i,se(k,3))=-springF(i,k)*sin(se(k,4));
                        end
                        match_r(:,1)=0;
                    end
                    Fs_x(i,:,1)=Fsl_x(i,:)+Fsr_x(i,:);
                    Fs_y(i,:,1)=Fsl_y(i,:)+Fsr_y(i,:);
                    Fs_z(i,:,1)=Fsl_z(i,:)+Fsr_z(i,:);
                    
                else
                    Fsl_x(i,1:pointNum(1))=0;
                    Fsl_y(i,1:pointNum(1))=0;
                    Fsl_z(i,1:pointNum(1))=0;
                    Fsr_x(i,1:pointNum(1))=0;
                    Fsr_y(i,1:pointNum(1))=0;
                    Fsr_z(i,1:pointNum(1))=0;
                    
                    
                    %����������
                    for k=ul(1,n-1)+1:ul(1,n)
                        for j=ul(1,n-1)+1:ul(1,n)
                            %             if n ==1
                            %                 ul_number = 1;
                            %             else
                            %                 ul_number = n-1;
                            %             end
                            %             for k=ul_number:ul(1,n)
                            %                 for j=ul_number:ul(1,n)
                            match_l(j,1)=se(j,2)==se(k,2);
                        end
                        match_lNum=senum(match_l);   %�������_��L����O�pNo�𒊏o
                        sN=size(match_lNum);
                        if sN(2)>1
                            matchse_l=[];
                            for s=1:sN(2)  %matchnum�ɓ����Ă���O�pNo��xy�����͂�matchf_A�ɕ\��
                                matchse_l(1,s)=springF(i,match_lNum(1,s))*cos(se(match_lNum(1,s),4))*cos(se(match_lNum(1,s),5));  %x
                                matchse_l(2,s)=springF(i,match_lNum(1,s))*cos(se(match_lNum(1,s),4))*sin(se(match_lNum(1,s),5));   %y
                                matchse_l(3,s)=springF(i,match_lNum(1,s))*sin(se(match_lNum(1,s),4));   %z
                            end
                            Fsl_x(i,se(k,2))=sum(matchse_l(1,:));
                            Fsl_y(i,se(k,2))=sum(matchse_l(2,:));
                            Fsl_z(i,se(k,2))=sum(matchse_l(3,:));
                            
                        else
                            Fsl_x(i,se(k,2))=springF(i,k)*cos(se(k,4))*cos(se(k,5));
                            Fsl_y(i,se(k,2))=springF(i,k)*cos(se(k,4))*sin(se(k,5));
                            Fsl_z(i,se(k,2))=springF(i,k)*sin(se(k,4));
                        end
                        match_l(:,1)=0;
                    end
                    
                    %�E��������
                    for k=ul(1,n-1)+1:ul(1,n)
                        for j=ul(1,n-1)+1:ul(1,n)
                            %             if n ==1
                            %                 ul_number = 1;
                            %             else
                            %                 ul_number = n-1;
                            %             end
                            %             for k=ul_number:ul(1,n)
                            %                 for j=ul(1,n-1):ul(1,n)
                            match_r(j,1)=se(j,3)==se(k,3);
                        end
                        match_rNum=senum(match_r);   %�������_��L����O�pNo�𒊏o
                        sN=size(match_rNum);
                        if sN(2)>1
                            matchse_r=[];
                            for s=1:sN(2)   %matchnum�ɓ����Ă���O�pNo��xy�����͂�matchf_A�ɕ\��
                                matchse_r(1,s)=-springF(i,match_rNum(1,s))*cos(se(match_rNum(1,s),4))*cos(se(match_rNum(1,s),5));  %x
                                matchse_r(2,s)=-springF(i,match_rNum(1,s))*cos(se(match_rNum(1,s),4))*sin(se(match_rNum(1,s),5));   %y
                                matchse_r(3,s)=-springF(i,match_rNum(1,s))*sin(se(match_rNum(1,s),4));   %z
                            end
                            Fsr_x(i,se(k,3))=sum(matchse_r(1,:));
                            Fsr_y(i,se(k,3))=sum(matchse_r(2,:));
                            Fsr_z(i,se(k,3))=sum(matchse_r(3,:));
                            
                        else
                            Fsr_x(i,se(k,3))=-springF(i,k)*cos(se(k,4))*cos(se(k,5));
                            Fsr_y(i,se(k,3))=-springF(i,k)*cos(se(k,4))*sin(se(k,5));
                            Fsr_z(i,se(k,3))=-springF(i,k)*sin(se(k,4));
                        end
                        match_r(:,1)=0;
                    end
                    Fs_x(i,:,n)=Fsl_x(i,:)+Fsr_x(i,:);
                    Fs_y(i,:,n)=Fsl_y(i,:)+Fsr_y(i,:);
                    Fs_z(i,:,n)=Fsl_z(i,:)+Fsr_z(i,:);
                end
            end
            
            FsSum_x(i,:)=Fs_x(i,:,1)+Fs_x(i,:,2)+Fs_x(i,:,3)+Fs_x(i,:,4)+Fs_x(i,:,5)+Fs_x(i,:,6);
            FsSum_y(i,:)=Fs_y(i,:,1)+Fs_y(i,:,2)+Fs_y(i,:,3)+Fs_y(i,:,4)+Fs_y(i,:,5)+Fs_y(i,:,6);
            FsSum_z(i,:)=Fs_z(i,:,1)+Fs_z(i,:,2)+Fs_z(i,:,3)+Fs_z(i,:,4)+Fs_z(i,:,5)+Fs_z(i,:,6);
            
            
            %% �S������
            
%             con_x(i,1:pointNum(1))=-1*conk*(datai_x(i,1:pointNum(1))-datai_x(1,1:pointNum(1)));
%             con_y(i,1:pointNum(1))=-1*conk*(datai_y(i,1:pointNum(1))-datai_y(1,1:pointNum(1)));
%             con_z(i,1:pointNum(1))=-1*conk*(datai_z(i,1:pointNum(1))-datai_z(1,1:pointNum(1)));
            
            if i==1
                vis_x(i,1:pointNum(1))=0;
                vis_y(i,1:pointNum(1))=0;
                vis_z(i,1:pointNum(1))=0;
            else
                vis_x(i,1:pointNum(1))=-1*conc*(datai_x(i,1:pointNum(1))-datai_x(i-1,1:pointNum(1)))/dt;
                vis_y(i,1:pointNum(1))=-1*conc*(datai_y(i,1:pointNum(1))-datai_y(i-1,1:pointNum(1)))/dt;
                vis_z(i,1:pointNum(1))=-1*conc*(datai_z(i,1:pointNum(1))-datai_z(i-1,1:pointNum(1)))/dt;
            end
            
            %% �̐ϕۑ���
            
            for j=1:tetraNum(1)/5
                for k=1:tetraNum(2)
                    if i==1
                        V0(k+(j-1)*tetraNum(2),1)=1/6*abs(dot(cross((data0(tetra(3+5*(j-1),k),:).'-data0(tetra(2+5*(j-1),k),:).'),(data0(tetra(4+5*(j-1),k),:).'-data0(tetra(2+5*(j-1),k),:).')),((data0(tetra(5+5*(j-1),k),:).'-data0(tetra(2+5*(j-1),k),:).'))));
                    end
                    
                    %                 if i == 1
                    %fvA�����߂�̂ɕK�v���޸��
                    Vec1(:,1+4*(k-1)+(j-1)*tetraNum(2)*4,i)=[datai_x(i,tetra(3+5*(j-1),k));datai_y(i,tetra(3+5*(j-1),k));datai_z(i,tetra(3+5*(j-1),k))]-[datai_x(i,tetra(5+5*(j-1),k));datai_y(i,tetra(5+5*(j-1),k));datai_z(i,tetra(5+5*(j-1),k))];   %B-D
                    Vec2(:,1+4*(k-1)+(j-1)*tetraNum(2)*4,i)=[datai_x(i,tetra(4+5*(j-1),k));datai_y(i,tetra(4+5*(j-1),k));datai_z(i,tetra(4+5*(j-1),k))]-[datai_x(i,tetra(5+5*(j-1),k));datai_y(i,tetra(5+5*(j-1),k));datai_z(i,tetra(5+5*(j-1),k))];   %C-D
                    
                    %fvB
                    Vec1(:,2+4*(k-1)+(j-1)*tetraNum(2)*4,i)=[datai_x(i,tetra(4+5*(j-1),k));datai_y(i,tetra(4+5*(j-1),k));datai_z(i,tetra(4+5*(j-1),k))]-[datai_x(i,tetra(2+5*(j-1),k));datai_y(i,tetra(2+5*(j-1),k));datai_z(i,tetra(2+5*(j-1),k))];   %C-A
                    Vec2(:,2+4*(k-1)+(j-1)*tetraNum(2)*4,i)=[datai_x(i,tetra(5+5*(j-1),k));datai_y(i,tetra(5+5*(j-1),k));datai_z(i,tetra(5+5*(j-1),k))]-[datai_x(i,tetra(2+5*(j-1),k));datai_y(i,tetra(2+5*(j-1),k));datai_z(i,tetra(2+5*(j-1),k))];   %D-A
                    
                    %fvC
                    Vec1(:,3+4*(k-1)+(j-1)*tetraNum(2)*4,i)=[datai_x(i,tetra(5+5*(j-1),k));datai_y(i,tetra(5+5*(j-1),k));datai_z(i,tetra(5+5*(j-1),k))]-[datai_x(i,tetra(3+5*(j-1),k));datai_y(i,tetra(3+5*(j-1),k));datai_z(i,tetra(3+5*(j-1),k))];   %D-B
                    Vec2(:,3+4*(k-1)+(j-1)*tetraNum(2)*4,i)=[datai_x(i,tetra(2+5*(j-1),k));datai_y(i,tetra(2+5*(j-1),k));datai_z(i,tetra(2+5*(j-1),k))]-[datai_x(i,tetra(3+5*(j-1),k));datai_y(i,tetra(3+5*(j-1),k));datai_z(i,tetra(3+5*(j-1),k))];   %A-B
                    
                    %fvD
                    Vec1(:,4+4*(k-1)+(j-1)*tetraNum(2)*4,i)=[datai_x(i,tetra(2+5*(j-1),k));datai_y(i,tetra(2+5*(j-1),k));datai_z(i,tetra(2+5*(j-1),k))]-[datai_x(i,tetra(4+5*(j-1),k));datai_y(i,tetra(4+5*(j-1),k));datai_z(i,tetra(4+5*(j-1),k))];   %A-C
                    Vec2(:,4+4*(k-1)+(j-1)*tetraNum(2)*4,i)=[datai_x(i,tetra(3+5*(j-1),k));datai_y(i,tetra(3+5*(j-1),k));datai_z(i,tetra(3+5*(j-1),k))]-[datai_x(i,tetra(4+5*(j-1),k));datai_y(i,tetra(4+5*(j-1),k));datai_z(i,tetra(4+5*(j-1),k))];   %B-C
                    
                    %                 else
                    %                     Vec1(:,1+4*(k-1)+(j-1)*tetraNum(2)*4,i)=Vec1(:,1+4*(k-1)+(j-1)*tetraNum(2)*4,i-1);
                    %                     Vec2(:,1+4*(k-1)+(j-1)*tetraNum(2)*4,i)=Vec2(:,1+4*(k-1)+(j-1)*tetraNum(2)*4,i-1);
                    %
                    %                     Vec1(:,2+4*(k-1)+(j-1)*tetraNum(2)*4,i)=Vec1(:,2+4*(k-1)+(j-1)*tetraNum(2)*4,i-1);
                    %                     Vec2(:,2+4*(k-1)+(j-1)*tetraNum(2)*4,i)=Vec2(:,2+4*(k-1)+(j-1)*tetraNum(2)*4,i-1);
                    %
                    %                     Vec1(:,3+4*(k-1)+(j-1)*tetraNum(2)*4,i)=Vec1(:,3+4*(k-1)+(j-1)*tetraNum(2)*4,i-1);
                    %                     Vec2(:,3+4*(k-1)+(j-1)*tetraNum(2)*4,i)=Vec2(:,3+4*(k-1)+(j-1)*tetraNum(2)*4,i-1);
                    %
                    %                     Vec1(:,4+4*(k-1)+(j-1)*tetraNum(2)*4,i)=Vec1(:,4+4*(k-1)+(j-1)*tetraNum(2)*4,i-1);
                    %                     Vec2(:,4+4*(k-1)+(j-1)*tetraNum(2)*4,i)=Vec2(:,4+4*(k-1)+(j-1)*tetraNum(2)*4,i-1);
                    %                 end
                    
                    
                    V(k+(j-1)*tetraNum(2),1)=1/6*abs(dot(cross((-1)*Vec2(:,3+4*(k-1)+(j-1)*tetraNum(2)*4,i),Vec1(:,2+4*(k-1)+(j-1)*tetraNum(2)*4,i)),Vec2(:,2+4*(k-1)+(j-1)*tetraNum(2)*4,i)));
                    
                    %�̐ώ��n��
                    if data(i,3)-data(1,3)>0    %���������W��,���ꂼ��̎l�ʑ̂��Ƃɑ̐ϕۑ��W����ݒ�
                        fvA(:,k+(j-1)*tetraNum(2))=1/6*springkVPk(2,k+(j-1)*tetraNum(2))*(V(k+(j-1)*tetraNum(2),1)-V0(k+(j-1)*tetraNum(2),1))/V0(k+(j-1)*tetraNum(2),1)^2*cross(Vec2(:,1+4*(k-1)+(j-1)*tetraNum(2)*4,i),Vec1(:,1+4*(k-1)+(j-1)*tetraNum(2)*4,i));
                        fvB(:,k+(j-1)*tetraNum(2))=1/6*springkVPk(2,k+(j-1)*tetraNum(2))*(V(k+(j-1)*tetraNum(2),1)-V0(k+(j-1)*tetraNum(2),1))/V0(k+(j-1)*tetraNum(2),1)^2*cross(Vec1(:,2+4*(k-1)+(j-1)*tetraNum(2)*4,i),Vec2(:,2+4*(k-1)+(j-1)*tetraNum(2)*4,i));
                        fvC(:,k+(j-1)*tetraNum(2))=1/6*springkVPk(2,k+(j-1)*tetraNum(2))*(V(k+(j-1)*tetraNum(2),1)-V0(k+(j-1)*tetraNum(2),1))/V0(k+(j-1)*tetraNum(2),1)^2*cross(Vec2(:,3+4*(k-1)+(j-1)*tetraNum(2)*4,i),Vec1(:,3+4*(k-1)+(j-1)*tetraNum(2)*4,i));
                        fvD(:,k+(j-1)*tetraNum(2))=1/6*springkVPk(2,k+(j-1)*tetraNum(2))*(V(k+(j-1)*tetraNum(2),1)-V0(k+(j-1)*tetraNum(2),1))/V0(k+(j-1)*tetraNum(2),1)^2*cross(Vec1(:,4+4*(k-1)+(j-1)*tetraNum(2)*4,i),Vec2(:,4+4*(k-1)+(j-1)*tetraNum(2)*4,i));
                    else
                        fvA(:,k+(j-1)*tetraNum(2))=1/6*springkVPk(2,k+(j-1)*tetraNum(2))*(V(k+(j-1)*tetraNum(2),1)-V0(k+(j-1)*tetraNum(2),1))/V0(k+(j-1)*tetraNum(2),1)^2*cross(Vec1(:,1+4*(k-1)+(j-1)*tetraNum(2)*4,i),Vec2(:,1+4*(k-1)+(j-1)*tetraNum(2)*4,i));
                        fvB(:,k+(j-1)*tetraNum(2))=1/6*springkVPk(2,k+(j-1)*tetraNum(2))*(V(k+(j-1)*tetraNum(2),1)-V0(k+(j-1)*tetraNum(2),1))/V0(k+(j-1)*tetraNum(2),1)^2*cross(Vec2(:,2+4*(k-1)+(j-1)*tetraNum(2)*4,i),Vec1(:,2+4*(k-1)+(j-1)*tetraNum(2)*4,i));
                        fvC(:,k+(j-1)*tetraNum(2))=1/6*springkVPk(2,k+(j-1)*tetraNum(2))*(V(k+(j-1)*tetraNum(2),1)-V0(k+(j-1)*tetraNum(2),1))/V0(k+(j-1)*tetraNum(2),1)^2*cross(Vec1(:,3+4*(k-1)+(j-1)*tetraNum(2)*4,i),Vec2(:,3+4*(k-1)+(j-1)*tetraNum(2)*4,i));
                        fvD(:,k+(j-1)*tetraNum(2))=1/6*springkVPk(2,k+(j-1)*tetraNum(2))*(V(k+(j-1)*tetraNum(2),1)-V0(k+(j-1)*tetraNum(2),1))/V0(k+(j-1)*tetraNum(2),1)^2*cross(Vec2(:,4+4*(k-1)+(j-1)*tetraNum(2)*4,i),Vec1(:,4+4*(k-1)+(j-1)*tetraNum(2)*4,i));
                    end
                end
            end
            
            FvA(1:pointNum(1),1:3)=0;
            tetranum=(1:tetraNum(1)*tetraNum(2)/5);
            %a�sb��ڂ�node�ԍ���L����tetraNum��1
            for m=1:tetraNum(1)/5
                for k=1:tetraNum(2)
                    for j=1:tetraNum(1)/5
                        matchA(1+(j-1)*tetraNum(2):tetraNum(2)*j,1)=tetra(2+5*(j-1),:).'==tetra(2+5*(m-1),k);    %�������_�̂���O�pNo���@����ꍇ�F1 �Ȃ��F0
                    end
                    matchAnum=tetranum(matchA);   %�������_��L����O�pNo�𒊏o
                    mN=size(matchAnum);      %1�ȏ�Ȃ瓯�����_��L����O�pNo���ق��ɂ���Ƃ�������
                    if mN(2)>1
                        matchf_A=[];
                        %�v���O�����������̂��߂̎��O���蓖��
                        matcn_A = zeros(3, mN(2));
                        for n=1:mN(2)   %matchnum�ɓ����Ă���O�pNo��xy�����͂�matchf_A�ɕ\��
                            matchf_A(1,n)=fvA(1,matchAnum(1,n));   %x
                            matchf_A(2,n)=fvA(2,matchAnum(1,n));   %y
                            matchf_A(3,n)=fvA(3,matchAnum(1,n));   %z
                        end
                        
                        FvA(tetra(2+5*(m-1),k),1)=sum(matchf_A(1,:));
                        FvA(tetra(2+5*(m-1),k),2)=sum(matchf_A(2,:));
                        FvA(tetra(2+5*(m-1),k),3)=sum(matchf_A(3,:));
                        
                    else
                        FvA(tetra(2+5*(m-1),k),1)=fvA(1,matchAnum);
                        FvA(tetra(2+5*(m-1),k),2)=fvA(2,matchAnum);
                        FvA(tetra(2+5*(m-1),k),3)=fvA(3,matchAnum);
                    end
                end
            end
            
            FvB(1:pointNum(1),1:3)=0;
            %a�sb��ڂ�node�ԍ���L����tetraNum��1
            for m=1:tetraNum(1)/5
                for k=1:tetraNum(2)
                    for j=1:tetraNum(1)/5
                        matchB(1+(j-1)*tetraNum(2):tetraNum(2)*j,1)=tetra(3+5*(j-1),:).'==tetra(3+5*(m-1),k);    %�������_�̂���O�pNo���@����ꍇ�F1 �Ȃ��F0
                    end
                    matchBnum=tetranum(matchB);   %�������_��L����O�pNo�𒊏o
                    mN=size(matchBnum);      %1�ȏ�Ȃ瓯�����_��L����O�pNo���ق��ɂ���Ƃ�������
                    if mN(2)>1
                        matchf_B=[];
                        %�v���O�����������̂��߂̎��O���蓖��
                        matcn_B = zeros(3, mN(2));
                        for n=1:mN(2)  %matchnum�ɓ����Ă���O�pNo��xy�����͂�matchf_B�ɕ\��
                            matchf_B(1,n)=fvB(1,matchBnum(1,n));   %x
                            matchf_B(2,n)=fvB(2,matchBnum(1,n));   %y
                            matchf_B(3,n)=fvB(3,matchBnum(1,n));   %z
                        end
                        
                        FvB(tetra(3+5*(m-1),k),1)=sum(matchf_B(1,:));
                        FvB(tetra(3+5*(m-1),k),2)=sum(matchf_B(2,:));
                        FvB(tetra(3+5*(m-1),k),3)=sum(matchf_B(3,:));
                        
                    else
                        FvB(tetra(3+5*(m-1),k),1)=fvB(1,matchBnum);
                        FvB(tetra(3+5*(m-1),k),2)=fvB(2,matchBnum);
                        FvB(tetra(3+5*(m-1),k),3)=fvB(3,matchBnum);
                    end
                end
            end
            
            FvC(1:pointNum(1),1:3)=0;
            %a�sb��ڂ�node�ԍ���L����tetraNum��1
            for m=1:tetraNum(1)/5
                for k=1:tetraNum(2)
                    for j=1:tetraNum(1)/5
                        matchC(1+(j-1)*tetraNum(2):tetraNum(2)*j,1)=tetra(4+5*(j-1),:).'==tetra(4+5*(m-1),k);    %�������_�̂���O�pNo���@����ꍇ�F1 �Ȃ��F0
                    end
                    matchCnum=tetranum(matchC);   %�������_��L����O�pNo�𒊏o
                    mN=size(matchCnum);      %1�ȏ�Ȃ瓯�����_��L����O�pNo���ق��ɂ���Ƃ�������
                    if mN(2)>1
                        matchf_C=[];
                        %�v���O�����������̂��߂̎��O���蓖��
                        matcn_C = zeros(3, mN(2));
                        for n=1:mN(2)   %matchnum�ɓ����Ă���O�pNo��xy�����͂�matchf_C�ɕ\��
                            matchf_C(1,n)=fvC(1,matchCnum(1,n));   %x
                            matchf_C(2,n)=fvC(2,matchCnum(1,n));   %y
                            matchf_C(3,n)=fvC(3,matchCnum(1,n));   %z
                        end
                        
                        FvC(tetra(4+5*(m-1),k),1)=sum(matchf_C(1,:));
                        FvC(tetra(4+5*(m-1),k),2)=sum(matchf_C(2,:));
                        FvC(tetra(4+5*(m-1),k),3)=sum(matchf_C(3,:));
                        
                    else
                        FvC(tetra(4+5*(m-1),k),1)=fvC(1,matchCnum);
                        FvC(tetra(4+5*(m-1),k),2)=fvC(2,matchCnum);
                        FvC(tetra(4+5*(m-1),k),3)=fvC(3,matchCnum);
                    end
                end
            end
            
            
            FvD(1:pointNum(1),1:3)=0;
            %a�sb��ڂ�node�ԍ���L����tetraNum��1
            for m=1:tetraNum(1)/5
                for k=1:tetraNum(2)
                    for j=1:tetraNum(1)/5
                        matchD(1+(j-1)*tetraNum(2):tetraNum(2)*j,1)=tetra(5+5*(j-1),:).'==tetra(5+5*(m-1),k);    %�������_�̂���O�pNo���@����ꍇ�F1 �Ȃ��F0
                    end
                    matchDnum=tetranum(matchD);   %�������_��L����O�pNo�𒊏o
                    mN=size(matchDnum);      %1�ȏ�Ȃ瓯�����_��L����O�pNo���ق��ɂ���Ƃ�������
                    if mN(2)>1
                        matchf_D=[];
                        %�������̂��߂̎��O���蓖��
                        matcn_D = zeros(3, mN(2));
                        for n=1:mN(2)   %matchnum�ɓ����Ă���O�pNo��xy�����͂�matchf_D�ɕ\��
                            matchf_D(1,n)=fvD(1,matchDnum(1,n));   %x
                            matchf_D(2,n)=fvD(2,matchDnum(1,n));   %y
                            matchf_D(3,n)=fvD(3,matchDnum(1,n));   %z
                        end
                        
                        FvD(tetra(5+5*(m-1),k),1)=sum(matchf_D(1,:));
                        FvD(tetra(5+5*(m-1),k),2)=sum(matchf_D(2,:));
                        FvD(tetra(5+5*(m-1),k),3)=sum(matchf_D(3,:));
                        
                    else
                        FvD(tetra(5+5*(m-1),k),1)=fvD(1,matchDnum);
                        FvD(tetra(5+5*(m-1),k),2)=fvD(2,matchDnum);
                        FvD(tetra(5+5*(m-1),k),3)=fvD(3,matchDnum);
                    end
                end
            end
            
            
            if i==1
                Fv_x(i,1:pointNum(1))=0;
                Fv_y(i,1:pointNum(1))=0;
                Fv_z(i,1:pointNum(1))=0;
            else
                
                Fv_x(i,1:pointNum(1)) = FvA(1:pointNum(1),1).'+FvB(1:pointNum(1),1).'+FvC(1:pointNum(1),1).'+FvD(1:pointNum(1),1).';
                Fv_y(i,1:pointNum(1)) = FvA(1:pointNum(1),2).'+FvB(1:pointNum(1),2).'+FvC(1:pointNum(1),2).'+FvD(1:pointNum(1),2).';
                Fv_z(i,1:pointNum(1)) = FvA(1:pointNum(1),3).'+FvB(1:pointNum(1),3).'+FvC(1:pointNum(1),3).'+FvD(1:pointNum(1),3).';
                
                %z�������̗͂�0�ɂ��āCx,y�������ɕ��U(xy���ʂɎʑ�)
                FvTotal(i,1:pointNum(1)) = sqrt(Fv_x(i,1:pointNum(1)).^2 + Fv_y(i,1:pointNum(1)).^2 + Fv_z(i,1:pointNum(1)).^2);
                alpha = atan(Fv_y(i,1:pointNum(1)) ./ Fv_x(i,1:pointNum(1)));
                alpha = alpha + pi.*(Fv_x(i,1:pointNum(1)) < 0);
                Fv_x(i,1:pointNum(1)) =  FvTotal(i,1:pointNum(1)).*cos(alpha);
                Fv_y(i,1:pointNum(1)) =  FvTotal(i,1:pointNum(1)).*sin(alpha);
                Fv_x(i,1:pointNum(1)) = fillmissing(Fv_x(i,1:pointNum(1)),'constant',0);
                Fv_y(i,1:pointNum(1)) = fillmissing(Fv_y(i,1:pointNum(1)),'constant',0);
                
            end
            
            %0�F���@1�F�c�@2�F�����@3�Fxy���ʎ΂߁@4�Fxz���ʎ΂߁@5�Fyz���ʎ΂߁@6�F��Ԏ΂�
            
            %% �d��
            %     for j=1:h
            %         Fg_x(i,y*t*(j-1)+1:y*t*j) = 0;
            %         Fg_y(i,y*t*(j-1)+1:y*t*j) = 0;
            %         Fg_z(i,y*t*(j-1)+1:y*t*j) = (h-j)*mass*gravityG;
            %     end
            
            
            %% ���_�ɉ����͂Ƒ��x�̎Z�o
            Fn_x(i,:)=FsSum_x(i,:)+Fv_x(i,:)+Mr_x(i,:)+con_x(i,:)+vis_x(i,:);%Fg_x(i,:)+Fv_x(i,:)+FsSum_x(i,:)+
            Fn_y(i,:)=FsSum_y(i,:)+Fv_y(i,:)+Mr_y(i,:)+con_y(i,:)+vis_y(i,:);%Fg_y(i,:)+Fv_y(i,:)+FsSum_y(i,:)
            Fn_z(i,:)=FsSum_z(i,:)+Fv_z(i,:)+Mr_z(i,:)+con_z(i,:)+vis_z(i,:);%+Fg_z(i,:)+Fv_z(i,:)+
            
            vn_x(i,1:pointNum(1))=vn_x(i,1:pointNum(1))+Fn_x(i,1:pointNum(1))/2*dt;
            vn_y(i,1:pointNum(1))=vn_y(i,1:pointNum(1))+Fn_y(i,1:pointNum(1))/2*dt;
            vn_z(i,1:pointNum(1))=vn_z(i,1:pointNum(1))+Fn_z(i,1:pointNum(1))/2*dt;
            
            volumeData(:,1) = datai_x(i,:)';
            volumeData(:,2) = datai_y(i,:)';
            volumeData(:,3) = datai_z(i,:)';
            
            for n=1:h-1
                for j=1:t-1
                    for k=1:y-1
                        tetraV1(k+(j-1)*(y-1)+(n-1)*(y-1)*(t-1)) = abs(dot((volumeData(n*t*y+k+(j-1)*y,:)-volumeData((n-1)*t*y+k+(j-1)*y,:)),cross((volumeData((n-1)*t*y+1+k+(j-1)*y,:)-volumeData((n-1)*t*y+k+(j-1)*y,:)),(volumeData((n-1)*t*y+k+j*y,:)-volumeData((n-1)*t*y+k+(j-1)*y,:)))))/6;
                        tetraV2(k+(j-1)*(y-1)+(n-1)*(y-1)*(t-1)) = abs(dot((volumeData(n*t*y+k+(j-1)*y,:)-volumeData(n*t*y+k+(j-1)*y+1,:)),cross((volumeData((n-1)*t*y+k+j*y,:)-volumeData(n*t*y+k+(j-1)*y+1,:)),(volumeData((n-1)*t*y+1+k+(j-1)*y,:)-volumeData(n*t*y+k+(j-1)*y+1,:)))))/6;
                        tetraV3(k+(j-1)*(y-1)+(n-1)*(y-1)*(t-1)) = abs(dot((volumeData(n*t*y+k+(j-1)*y,:)-volumeData(n*t*y+k+(j-1)*y+1,:)),cross((volumeData((n-1)*t*y+k+j*y,:)-volumeData(n*t*y+k+(j-1)*y+1,:)),(volumeData(n*t*y+k+j*y,:)-volumeData(n*t*y+k+(j-1)*y+1,:)))))/6;
                        tetraV4(k+(j-1)*(y-1)+(n-1)*(y-1)*(t-1)) = abs(dot((volumeData(n*t*y+k+j*y+1,:)-volumeData((n-1)*t*y+k+j*y+1,:)),cross((volumeData((n-1)*t*y+k+j*y,:)-volumeData((n-1)*t*y+k+j*y+1,:)),(volumeData((n-1)*t*y+1+k+(j-1)*y,:)-volumeData((n-1)*t*y+k+j*y+1,:)))))/6;
                        tetraV5(k+(j-1)*(y-1)+(n-1)*(y-1)*(t-1)) = abs(dot((volumeData(n*t*y+k+j*y+1,:)-volumeData((n-1)*t*y+1+k+(j-1)*y,:)),cross((volumeData(n*t*y+k+(j-1)*y+1,:)-volumeData((n-1)*t*y+1+k+(j-1)*y,:)),(volumeData((n-1)*t*y+k+j*y,:)-volumeData((n-1)*t*y+1+k+(j-1)*y,:)))))/6;
                        tetraV6(k+(j-1)*(y-1)+(n-1)*(y-1)*(t-1)) = abs(dot((volumeData(n*t*y+k+j*y+1,:)-volumeData((n-1)*t*y+k+j*y,:)),cross((volumeData(n*t*y+k+(j-1)*y+1,:)-volumeData((n-1)*t*y+k+j*y,:)),(volumeData((n-1)*t*y+k+j*y,:)-volumeData(n*t*y+k+(j-1)*y,:)))))/6;
                    end
                end
            end
            sumTetraV1 = sum(tetraV1);
            sumTetraV2 = sum(tetraV2);
            sumTetraV3 = sum(tetraV3);
            sumTetraV4 = sum(tetraV4);
            sumTetraV5 = sum(tetraV5);
            sumTetraV6 = sum(tetraV6);
            volume(i) = sumTetraV1 + sumTetraV2 + sumTetraV3 + sumTetraV4 + sumTetraV5 + sumTetraV6;
            volumeRatio(i) = volume(i)/volume(1);
            plot([1:i],volumeRatio);
            Frame(i) = getframe(1);
            if volumeRatio(i) > 2.0
                break
            end
            
        end
        toc
        volumeRatio(end-1:end)=[]
        figure1 = figure()
        plot([1:(i-2)],volumeRatio);
%         graphName = "springk=" + searchListK(searchNK) + "_c=" + searchListC(searchNC);
        graphName = "springk=" + searchNK + "_c=" + searchNC;
        saveas(figure1,graphName);
        VolumeHistroy(searchNC, searchNK) = volumeRatio(end);
        volumeRatio = [];
        close
    end
end
F_fib=springF(:,ul(1,2)+1:ul(1,2)+y*t*(h-1));

datai_x_name=strcat('output\', muscle_name, '_data_x.csv');
datai_y_name=strcat('output\', muscle_name, '_data_y.csv');
datai_z_name=strcat('output\', muscle_name, '_data_z.csv');
csvwrite(datai_x_name,datai_x);
csvwrite(datai_y_name,datai_y);
csvwrite(datai_z_name,datai_z);
csvwrite("VolumeHistory.csv",VolumeHistroy);

% v = VideoWriter('model_3d_7.avi');
%
% open(v);
% writeVideo(v,Frame);
% close(v);
%


