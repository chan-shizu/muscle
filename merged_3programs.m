clear
muscleNames = loadMuscleName();

% �e�����̍ő咸�_�� = Ni,Nj,Nk
%Ni=16,Nj=16��O��Ƃ��āAjava�R�[�h�������Ă���
%�����������߂�
Ni = 5;
Nj = 5;
Nk = 6;

%%surface2grid_circle_n
for m =1
    
    m
    %**********************
    % 1.  load data
    %**********************
    %muscleName = ['L_',muscleNames{m}];
    muscleName = [muscleNames{m}]
    
    %filename = ['traced\', muscleName, '.csv'];%defmuscle����o�͂��ꂽcsv����ёւ�������
    filename = ['C:\Users\kou_0\OneDrive\�h�L�������g\����matlab\traced\', muscleName, '.csv'];
%     filename = ['C:\Users\bubbl\Documents\shizuya_M1\DefMuscle_for_TUS_2018\matlab\surface2grid_0608\traced\', muscleName, '.csv'];
    [NUM, TXT, RAW] = xlsread(filename)%TXT��csv�t�@�C���̒��̃e�L�X�g�f�[�^����?
    
    gridInfo = TXT(:,1:2);%1//2, 0//0�Ƃ��̏��Ccsv��2, 3�s��
    XYZ = NUM(:, 4:end);
    [nFibers, numq] = size(XYZ);%nFibers�͎蓮�Œ��o����fiber���Cnumq�͊e���_��xyz���W�̐��C3�Ŋ���ƁC���ꂼ���fiber�őł����_�̐�
    %maxSlices = numq/3;
    
    %**********************
    % 2. �X���C�X�������낦��
    %**********************
    for iFiber = 1:nFibers
        maxq = 0;
        for q = 1:numq
            if isnan(XYZ(iFiber, q))
                break
            else
                maxq = q;
            end
        end
        nSlices = maxq/3;%���ꂼ���fiber�ł������_�̐�
        
        fiberCoordinate = zeros(nSlices, 3);%3=xyz
        for iSlice = 1:nSlices
            fiberCoordinate(iSlice, 1) = XYZ(iFiber, 3*(iSlice-1)+1);
            fiberCoordinate(iSlice, 2) = XYZ(iFiber, 3*(iSlice-1)+2);
            fiberCoordinate(iSlice, 3) = XYZ(iFiber, 3*(iSlice-1)+3);
            %fiberCoordinate(iSlice, 3) = XYZ(iFiber, 3*(iSlice-1)+3)+0.002;%��r�؂̉��������͏����O�Ɉړ������������S
        end %���ꂼ���fiber�ɂ����čs��̍s�Ɨ�����ւ�
        
        % fiberCoordinate =
        %     0.1360    1.3221    0.1013
        %     0.1397    1.3138    0.1011
        %     0.1428    1.3076    0.1001
        %     0.1446    1.3036    0.0993
        %     0.1470    1.2984    0.0980
        %     0.1507    1.2900    0.0957
        %     0.1542    1.2820    0.0934
        %     0.1574    1.2735    0.0896
        %     0.1614    1.2639    0.0863
        
        iFiber
        fiberCoordinate
        lengthCoordinate = calcLengthCoordinate(fiberCoordinate)%�ŏ��̓_��0�Ƃ��āC���O�̓_����̃L������ݐ�
        totalL = lengthCoordinate(end);%�ŏ��̓_����Ō�̓_�܂ł́C���O�̓_����̃L�����̗ݐ�
        newLengthCoordinate = [0:  totalL/(Nk-1)  :totalL]';%(Nk-1)�̃Z�O�����g�ɕ����@�w�肵���������œ��Ԋu�ɔz������
        %size(newLengthCoordinate) must = Nk
        
        %yy = splineAllColumn(x,y,xx)
        newFiberCoordinate = splineAllColumn(lengthCoordinate, fiberCoordinate, newLengthCoordinate);%���ꂼ���fiber�Ŏw�肵��z���̕�����(Nk)�̐��Ɏ��_���X�v���C����Ԃ���
        
        yy = newFiberCoordinate;
        plot3(yy(:,1), yy(:,2), yy(:,3), 'b.' );hold on
        axis equal;
        
        fiber{iFiber}.coordinate = newFiberCoordinate;%�\���̂��Ă�����ł�����?�I�u�W�F�N�g�w���I�ȁCfiber�ɂ́Ci,j,newFiberCoordinate�����ꂼ���fiber���ƂɊ܂܂�Ă�
        
        %i
        str = gridInfo{iFiber,1};
        k = strfind(str, '//');
        a = str(1:k-1);
        if  strcmp(a, '0')%0//0�Ȃ�C0
            fiber{iFiber}.i = '0';
        else
            b = str(k+2:end);%1//2�Ƃ��Ȃ�C1/2
            fiber{iFiber}.i = [a, '/', b];
        end
        
        %j
        str = gridInfo{iFiber,2};
        k = strfind(str, '//');
        a = str(1:k-1);
        if  strcmp(a, '0')
            fiber{iFiber}.j = '0';
        else
            b = str(k+2:end);
            fiber{iFiber}.j = [a, '/', b];
        end
        
        
    end % for iFiber = 1:nFibers
    
    
    %****************************************************
    % 3. Ni=16, Nj=7�̕\�ʃt�@�C�o�[���쐬
    %****************************************************
    clear k
    
    for k=1:Nk
        
        for side =1:4
            %��Ӗځij=0��i��������������j
            IJXYZ = [];
            for iFiber = 1:nFibers
                
                if side == 1 %west  J=0�̕�
                    index = str2num(fiber{iFiber}.j);
                    num01 = 0;%str 0 or 1
                    sort12 = 1;% 1(i) or 2(j)�Ń\�[�g
                elseif side == 2 %south  i=0�̕�
                    index = str2num(fiber{iFiber}.i);
                    num01 = 0;
                    sort12 = 2;
                elseif side == 3 %east  J=1�̕�
                    index = str2num(fiber{iFiber}.j);
                    num01 = 1;
                    sort12 = 1;
                elseif side == 4 %north  i=1�̕�
                    index = str2num(fiber{iFiber}.i);
                    num01 = 1;
                    sort12 = 2;
                end
                
                %�Y����fiber����������A�����~�ς��Ă���
                if (index == num01)
                    tmp = [str2num(fiber{iFiber}.i) ,... //������i
                        str2num(fiber{iFiber}.j),...  //������j
                        fiber{iFiber}.coordinate(k,1:3)  ];%k�X���C�X�ڂ�xyz���W
                    IJXYZ = [IJXYZ; tmp];
                    % ��jIJXYZ =
                    %          0            0    0.1871    1.1710    0.0719
                    %     1.0000         0    0.1861    1.1852    0.0674
                    %     0.5000         0    0.1861    1.1757    0.0684
                end
            end
            
            IJXYZ
            [Y, I] = sort(IJXYZ(:,sort12));%i�Ń\�[�g
            IJXYZ = IJXYZ(I,:)
            %     IJXYZ =
            %          0            0    0.1871    1.1710    0.0719
            %     0.5000         0    0.1861    1.1757    0.0684
            %     1.0000         0    0.1861    1.1852    0.0674
            coordinate = IJXYZ(:, 3:5);
            lengthCoordinate = calcLengthCoordinate(coordinate);
            totalL = lengthCoordinate(end);
            if (side==1 || side==3)
                newLengthCoordinate = [0:  totalL/(Ni-1)  :totalL]';%(Ni-1)�̃Z�O�����g�ɕ���
            else
                newLengthCoordinate = [0:  totalL/(Nj-1)  :totalL]';%(Nj-1)�̃Z�O�����g�ɕ���
            end
            
            newIJXYZ{side, k} = splineAllColumn(lengthCoordinate, IJXYZ, newLengthCoordinate)
            size(newIJXYZ{side, k})
            
        end % for side =1:4
        
    end %for k=1:Nk
    
    
    %************************************************
    %newIJXYZ{side, k}��csv�t�@�C���ɕۑ�
    %************************************************
    sfilename = ['C:\Users\kou_0\OneDrive\�h�L�������g\����matlab\aligned\', muscleName, '_aligned.csv'];
%     sfilename = ['C:\Users\bubbl\Documents\shizuya_M1\DefMuscle_for_TUS_2018\matlab\surface2grid_0608\aligned\', muscleName,  '_aligned.csv'];
    fid = fopen(sfilename,'wt'); % �������ݗp�Ƀt�@�C���I�[�v��
    for side =1:4
        
        if (side==1 || side==3)
            N = Ni;%16
        else
            N =  Nj;%7
        end
        
        for n = 1:N %1~16 or 1~7
            if (n==1 || n==N)
                fprintf(fid, '%s,',  'true');
            else
                fprintf(fid, '%s,',  'false');
            end
            
            if side == 1
                fprintf(fid, '%s,',   [num2str(n-1), '//',  num2str(N-1)] );
                fprintf(fid, '%s,',  '0//0');
            elseif side == 2
                fprintf(fid, '%s,',  '0//0');
                fprintf(fid, '%s,',   [num2str(n-1), '//',  num2str(N-1)] );
            elseif side == 3
                fprintf(fid, '%s,',   [num2str(n-1), '//',  num2str(N-1)] );
                fprintf(fid, '%s,',  '1//1');
            elseif side == 4
                fprintf(fid, '%s,',  '1//1');
                fprintf(fid, '%s,',   [num2str(n-1), '//',  num2str(N-1)] );
            end
            
            for k=1:Nk
                fprintf(fid, '%2.9f,',  newIJXYZ{side, k}(n,3) );%x���W
                fprintf(fid, '%2.9f,',  newIJXYZ{side, k}(n,4) );%y���W
                fprintf(fid, '%2.9f,',  newIJXYZ{side, k}(n,5) );%z���W
            end
            fprintf(fid, '\n');
            
        end %n = 1:N
        
    end
    
    fclose(fid);
    display(muscleName)
    
    
end

%% muscle_make & rotation
file_name_readcsv = ['C:\Users\kou_0\OneDrive\�h�L�������g\����matlab\aligned\',muscleName, '_aligned.csv']
% file_name_readcsv = strcat('C:\Users\bubbl\Documents\shizuya_M1\DefMuscle_for_TUS_2018\matlab\surface2grid_0608\aligned\', muscleName, '_aligned.csv')

file = csvread(file_name_readcsv,0,3);
fileNum=size(file)

y = Ni; 
t = Nj; 
h = Nk; 

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

file_n(1,:)=file(16,:);
file_n(2,:)=file(13,:);
file_n(3,:)=file(12,:);
file_n(4,:)=file(11,:);
file_n(5,:)=file(10,:);
file_n(6,:)=file(15,:);
file_n(7,:)=file(9,:);
file_n(8,:)=file(14,:);
file_n(9,:)=file(8,:);
file_n(10,:)=file(5,:);
file_n(11,:)=file(7,:);
file_n(12,:)=file(4,:);
file_n(13,:)=file(3,:);
file_n(14,:)=file(2,:);
file_n(15,:)=file(1,:);
file_n(16,:)=file(6,:);

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

% file_n(1,:)=file(1,:);
% file_n(2,:)=file(2,:);
% file_n(3,:)=file(3,:);
% file_n(4,:)=file(4,:);
% file_n(5,:)=file(8,:);
% file_n(6,:)=file(5,:);
% file_n(7,:)=file(6,:);
% file_n(8,:)=file(7,:);

%�����w�̂��ꂼ���fiber�̎��_��z�������̍����𕽋ςł��낦�Ă�
for i=1:h
    file_n(1:y*4-4,3*i)=mean(file_n(1:y*4-4,3*i));
end


%3��ɂȂ�悤�ɕ��ёւ�
for i=1:fileNum(2)/3
    mix(1+(i-1)*(y*4-4):(y*4-4)*i,1:3)=file_n(1:y*4-4,1+3*(i-1):3*i);
end

file_name_mix = ['C:\Users\kou_0\OneDrive\�h�L�������g\����matlab\muscle\',muscleName, '_min_n.csv']
% file_name_mix = strcat("C:\Users\bubbl\Documents\shizuya_M1\DefMuscle_for_TUS_2018\matlab\surface2grid_0608\muscle\", muscleName, "_min_n.csv")
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
if file_n(1,h*3+1) == 0
    file_n(:,h*3+1) = []
end

file_mix=file_n+file_p;

for i=1:fileNum(2)/3
    muscle(1+(i-1)*y*t:y*t*i,1:3)=file_mix(1:y*t,1+3*(i-1):3*i);
end

file_name_muscle = ['C:\Users\kou_0\OneDrive\�h�L�������g\����matlab\muscle\',muscleName, '_min_final.csv']
% file_name_muscle = strcat("C:\Users\bubbl\Documents\shizuya_M1\DefMuscle_for_TUS_2018\matlab\surface2grid_0608\muscle\", muscleName, "_min_final.csv")
csvwrite(file_name_muscle, muscle)

%% data�ړ�
%muscle = csvread("C:\Users\bubbl\OneDrive\�h�L�������g\shizuya_M1\DefMuscle_for_TUS_2018\matlab\surface2grid_0608\muscle\min\data0.csv",0,0);
for i=1:t*y
    data_move(1,1+3*(i-1):3*i)=muscle(t*y*h-t*y+i,1:3);
end

file_name_data_move = ['C:\Users\kou_0\OneDrive\�h�L�������g\����matlab\muscle\data_',muscleName, '.csv']
% file_name_data_move = strcat("C:\Users\bubbl\Documents\shizuya_M1\DefMuscle_for_TUS_2018\matlab\surface2grid_0608\muscle\data_", muscleName, ".csv")
csvwrite(file_name_data_move, data_move)

%% parameter
E=5*10^4; %file:///C:/Users/%E5%8D%A0%E9%83%A8%E3%80%80%E9%BA%BB%E9%87%8C%E5%AD%90/Downloads/nagano_05-02-05%20(1).pdf
v=0.49;
%k=8500;
%k=8237;

%ks=166.7558;
%ks=170;
Fmax=5;
PEsh = 4;
PExm = 0.4;

% file_name_data0 = strcat("C:\Users\bubbl\Documents\shizuya_M1\DefMuscle_for_TUS_2018\matlab\surface2grid_0608\muscle\min\", muscleName,"_min_final.csv")
% data0=csvread(file_name_data0, 0, 0)/1000; %mm����m�ɕϊ�
data0 = muscle

pointNum = size(data0);
dt=1*10^(-3);
%% se�쐬
for i=1:t
    yoko(1+(y-1)*(i-1):(y-1)*i,2)=1+y*(i-1):y*i-1;
    yoko(1+(y-1)*(i-1):(y-1)*i,3)=2+y*(i-1):y*i;
end

%��
for j=0:h-1
    for i=1:t
        yoko(1+(y-1)*(i-1)+j*(y-1)*t:(y-1)*i+j*(y-1)*t,2)=1+y*(i-1)+t*y*j:y*i-1+t*y*j;
        yoko(1+(y-1)*(i-1)+j*(y-1)*t:(y-1)*i+j*(y-1)*t,3)=2+y*(i-1)+t*y*j:y*i+t*y*j;
    end
end

%�c
for i=0:h-1
    tate(1+y*(t-1)*i:y*(t-1)*(i+1),2)=1+t*y*i:y*(t-1)+t*y*i;
    tate(1+y*(t-1)*i:y*(t-1)*(i+1),3)=1+y+t*y*i:t*y*(i+1);
end
tate(:,1)=1;

%����
takasa(1:(h-1)*t*y,2)=1:(h-1)*t*y;
takasa(1:(h-1)*t*y,3)=1+t*y:h*t*y;
takasa(:,1)=2;

%xy
for k=0:h-1
    for j=0:t-2
        for i=0:y-2
            xy(1+2*i+2*(y-1)*j+((y-1)*(t-1)*2*k),2)=1+i+y*j+t*y*k; xy(1+2*i+2*(y-1)*j+((y-1)*(t-1)*2*k),3)=2+y+i+y*j+t*y*k;
            xy(2+2*i+2*(y-1)*j+((y-1)*(t-1)*2*k),2)=2+i+y*j+t*y*k; xy(2+2*i+2*(y-1)*j+((y-1)*(t-1)*2*k),3)=1+y+i+y*j+t*y*k;
        end
    end
end
xy(:,1)=3;

%xz
for k=0:t-1
    for j=0:h-2
        for i=0:y-2
            xz(1+2*i+2*(y-1)*j+((y-1)*(h-1)*2*k),2)=1+i+y*t*j+y*k;    xz(1+2*i+2*(y-1)*j+((y-1)*(h-1)*2*k),3)=2+y*t+i+y*t*j+y*k;
            xz(2+2*i+2*(y-1)*j+((y-1)*(h-1)*2*k),2)=2+i+y*t*j+y*k;    xz(2+2*i+2*(y-1)*j+((y-1)*(h-1)*2*k),3)=1+y*t+i+y*t*j+y*k;
        end
    end
end
xz(:,1)=4;

%yz
for k=0:y-1
    for j=0:h-2
        for i=0:t-2
            yz(1+2*i+2*(t-1)*j+((t-1)*(h-1)*2*k),2)=1+t*i+t*y*j+k;      yz(1+2*i+2*(t-1)*j+((t-1)*(h-1)*2*k),3)=1+y*(t+1)+t*i+t*y*j+k;
            yz(2+2*i+2*(t-1)*j+((t-1)*(h-1)*2*k),2)=1+y+t*i+t*y*j+k;    yz(2+2*i+2*(t-1)*j+((t-1)*(h-1)*2*k),3)=1+y*t+t*i+t*y*j+k;
        end
    end
end
yz(:,1)=5;


se=vertcat(yoko,tate,takasa,xy,xz,yz);
seNum=size(se);
file_name_se = ['C:\Users\kou_0\OneDrive\�h�L�������g\����matlab\parameter\',muscleName, '_se.csv']
% file_name_se = strcat("C:\Users\bubbl\Documents\shizuya_M1\DefMuscle_for_TUS_2018\matlab\surface2grid_0608\parameter\", muscleName, "_se.csv")
csvwrite(file_name_se, se);
%% tetra
for j=1:(t-1)*(h-1)
    for i=1:y-1
        tetra(1+5*(j-1),1+10*(i-1):10*i)=i;
    end
end

%���ɍ���Ă���
for i=1:y-1
    tetra(2,1+10*(i-1))=1+(i-1);     tetra(3,1+10*(i-1))=1+y+y*t+(i-1);   tetra(4,1+10*(i-1))=1+y*t+(i-1);     tetra(5,1+10*(i-1))=2+y*t+(i-1);
    tetra(2,2+10*(i-1))=2+y+(i-1);   tetra(3,2+10*(i-1))=2+y+y*t+(i-1);   tetra(4,2+10*(i-1))=1+y+y*t+(i-1);   tetra(5,2+10*(i-1))=2+y*t+(i-1);
    tetra(2,3+10*(i-1))=1+(i-1);     tetra(3,3+10*(i-1))=2+(i-1);         tetra(4,3+10*(i-1))=2+y+(i-1);       tetra(5,3+10*(i-1))=2+y*t+(i-1);
    tetra(2,4+10*(i-1))=1+(i-1);     tetra(3,4+10*(i-1))=1+y+y*t+(i-1);   tetra(4,4+10*(i-1))=2+y+(i-1);       tetra(5,4+10*(i-1))=1+y+(i-1);
    tetra(2,5+10*(i-1))=1+(i-1);     tetra(3,5+10*(i-1))=2+y+(i-1);       tetra(4,5+10*(i-1))=1+y+y*t+(i-1);   tetra(5,5+10*(i-1))=2+y*t+(i-1);%���₵��tetra(2,5+10*(i-1))=1+(i-1);     tetra(3,5+10*(i-1))=2+y+(i-1);       tetra(4,5+10*(i-1))=1+y+y*t+(i-1);   tetra(5,5+10*(i-1))=2+y*t+(i-1);
    tetra(2,6+10*(i-1))=2+(i-1);     tetra(3,6+10*(i-1))=1+y*t+(i-1);     tetra(4,6+10*(i-1))=2+y*t+(i-1);     tetra(5,6+10*(i-1))=2+y+y*t+(i-1);
    tetra(2,7+10*(i-1))=1+y+(i-1);   tetra(3,7+10*(i-1))=1+y+y*t+(i-1);   tetra(4,7+10*(i-1))=1+y*t+(i-1);     tetra(5,7+10*(i-1))=2+y+y*t+(i-1);
    tetra(2,8+10*(i-1))=1+(i-1);     tetra(3,8+10*(i-1))=2+(i-1);         tetra(4,8+10*(i-1))=1+y+(i-1);       tetra(5,8+10*(i-1))=1+y*t+(i-1);
    tetra(2,9+10*(i-1))=2+(i-1);     tetra(3,9+10*(i-1))=2+y+(i-1);       tetra(4,9+10*(i-1))=1+y+(i-1);       tetra(5,9+10*(i-1))=2+y+y*t+(i-1);
    tetra(2,10+10*(i-1))=2+(i-1);    tetra(3,10+10*(i-1))=1+y+(i-1);      tetra(4,10+10*(i-1))=1+y*t+(i-1);    tetra(5,10+10*(i-1))=2+y+y*t+(i-1);
end


%�c
for i=1:t-2
    tetra(2+5*i:5+5*i,1:10*(y-1))=tetra(2+5*(i-1):5+5*(i-1),1:10*(y-1))+y;
end

%����
for i=1:h-2
    for k=1:t-1
        tetra(2+5*(k-1)+5*(t-1)*i:5+5*(k-1)+5*(t-1)*i,1:10*(y-1))=tetra(2+5*(k-1)+5*(t-1)*(i-1):5+5*(k-1)+5*(t-1)*(i-1),1:10*(y-1))+y*t;
    end
end
tetraNum=size(tetra);
file_name_tetra = ['C:\Users\kou_0\OneDrive\�h�L�������g\����matlab\parameter\',muscleName, '_tetra.csv']
% file_name_tetra = strcat("C:\Users\bubbl\Documents\shizuya_M1\DefMuscle_for_TUS_2018\matlab\surface2grid_0608\parameter\", muscleName, "_tetra.csv")
csvwrite(file_name_tetra, tetra);

%% 3D
tetraNum=size(tetra);
for j=1:tetraNum(1)/5;
    for i=1:tetraNum(2);
        x(i+(j-1)*tetraNum(2),1)=tetra(1+5*(j-1),i);  x(i+(j-1)*tetraNum(2),2)=data0(tetra(2+5*(j-1),i),1);  x(i+(j-1)*tetraNum(2),3)=data0(tetra(3+5*(j-1),i),1);  x(i+(j-1)*tetraNum(2),4)=data0(tetra(4+5*(j-1),i),1);  x(i+(j-1)*tetraNum(2),5)=data0(tetra(5+5*(j-1),i),1);
        y(i+(j-1)*tetraNum(2),1)=tetra(1+5*(j-1),i);  y(i+(j-1)*tetraNum(2),2)=data0(tetra(2+5*(j-1),i),2);  y(i+(j-1)*tetraNum(2),3)=data0(tetra(3+5*(j-1),i),2);  y(i+(j-1)*tetraNum(2),4)=data0(tetra(4+5*(j-1),i),2);  y(i+(j-1)*tetraNum(2),5)=data0(tetra(5+5*(j-1),i),2);
        z(i+(j-1)*tetraNum(2),1)=tetra(1+5*(j-1),i);  z(i+(j-1)*tetraNum(2),2)=data0(tetra(2+5*(j-1),i),3);  z(i+(j-1)*tetraNum(2),3)=data0(tetra(3+5*(j-1),i),3);  z(i+(j-1)*tetraNum(2),4)=data0(tetra(4+5*(j-1),i),3);  z(i+(j-1)*tetraNum(2),5)=data0(tetra(5+5*(j-1),i),3);
        V0(i+(j-1)*tetraNum(2),1)=1/6*dot(cross((data0(tetra(3+5*(j-1),i),:).'-data0(tetra(2+5*(j-1),i),:).'),(data0(tetra(4+5*(j-1),i),:).'-data0(tetra(2+5*(j-1),i),:).')),((data0(tetra(5+5*(j-1),i),:).'-data0(tetra(2+5*(j-1),i),:).')));
        
    end
end

%KCST�����߂�
D=E/(1+v)/(1-2*v)*[1-v v v 0 0 0;
    v 1-v v 0 0 0
    v v 1-v 0 0 0
    0 0 0 1/2-v 0 0
    0 0 0 0 1/2-v 0
    0 0 0 0 0 1/2-v];

for i=1:tetraNum(2)*tetraNum(1)/5;
    a1(i,1)=y(i,2)*(z(i,4)-z(i,3))-y(i,3)*(z(i,4)-z(i,2))+y(i,4)*(z(i,3)-z(i,2));
    a2(i,1)=-y(i,1)*(z(i,4)-z(i,3))+y(i,3)*(z(i,4)-z(i,1))-y(i,4)*(z(i,3)-z(i,1));
    a3(i,1)=y(i,1)*(z(i,4)-z(i,2))-y(i,2)*(z(i,4)-z(i,1))+y(i,4)*(z(i,2)-z(i,1));
    a4(i,1)=-y(i,1)*(z(i,3)-z(i,2))+y(i,2)*(z(i,3)-z(i,1))-y(i,3)*(z(i,2)-z(i,1));
    
    b1(i,1)=-x(i,2)*(z(i,4)-z(i,3))+x(i,3)*(z(i,4)-z(i,2))-x(i,4)*(z(i,3)-z(i,2));
    b2(i,1)=x(i,1)*(z(i,4)-z(i,3))-x(i,3)*(z(i,4)-z(i,1))+x(i,4)*(z(i,3)-z(i,1));
    b3(i,1)=-x(i,1)*(z(i,4)-z(i,2))+x(i,2)*(z(i,4)-z(i,1))-x(i,4)*(z(i,2)-z(i,1));
    b4(i,1)=x(i,1)*(z(i,3)-z(i,2))-x(i,2)*(z(i,3)-z(i,1))+x(i,3)*(z(i,2)-z(i,1));
    
    c1(i,1)=x(i,2)*(y(i,4)-y(i,3))-x(i,3)*(y(i,4)-y(i,2))+x(i,4)*(y(i,3)-y(i,2));
    c2(i,1)=-x(i,1)*(y(i,4)-y(i,3))+x(i,3)*(y(i,4)-y(i,1))-x(i,4)*(y(i,3)-y(i,1));
    c3(i,1)=x(i,1)*(y(i,4)-y(i,2))-x(i,2)*(y(i,4)-y(i,1))+x(i,4)*(y(i,2)-y(i,1));
    c4(i,1)=-x(i,1)*(y(i,3)-y(i,2))+x(i,2)*(y(i,3)-y(i,1))-x(i,3)*(y(i,2)-y(i,1));
    
    B(:,:,i)=1/(6*V0(i,1))*[a1(i,1) a2(i,1) a3(i,1) a4(i,1) 0 0 0 0 0 0 0 0
        0 0 0 0 b1(i,1) b2(i,1) b3(i,1) b4(i,1) 0 0 0 0
        0 0 0 0 0 0 0 0 c1(i,1) c2(i,1) c3(i,1) c4(i,1)
        b1(i,1) b2(i,1) b3(i,1) b4(i,1) a1(i,1) a2(i,1) a3(i,1) a4(i,1) 0 0 0 0
        0 0 0 0 c1(i,1) c2(i,1) c3(i,1) c4(i,1) b1(i,1) b2(i,1) b3(i,1) b4(i,1)
        c1(i,1) c2(i,1) c3(i,1) c4(i,1) 0 0 0 0 a1(i,1) a2(i,1) a3(i,1) a4(i,1)];
    
    KCST(:,:,i)=V0(i,1)*B(:,:,i).'*D*B(:,:,i);
end


%KMSM�����߂�
for i=1:tetraNum(2)*tetraNum(1)/5;
    H12(:,:,i)=1/(((x(i,1)-x(i,2))^2+(y(i,1)-y(i,2))^2+(z(i,1)-z(i,2))^2)^0.5)^2*[x(i,1)-x(i,2); y(i,1)-y(i,2); z(i,1)-z(i,2)]*[x(i,1)-x(i,2); y(i,1)-y(i,2); z(i,1)-z(i,2);].';
    H13(:,:,i)=1/(((x(i,1)-x(i,3))^2+(y(i,1)-y(i,3))^2+(z(i,1)-z(i,3))^2)^0.5)^2*[x(i,1)-x(i,3); y(i,1)-y(i,3); z(i,1)-z(i,3)]*[x(i,1)-x(i,3); y(i,1)-y(i,3); z(i,1)-z(i,3);].';
    H14(:,:,i)=1/(((x(i,1)-x(i,4))^2+(y(i,1)-y(i,4))^2+(z(i,1)-z(i,4))^2)^0.5)^2*[x(i,1)-x(i,4); y(i,1)-y(i,4); z(i,1)-z(i,4)]*[x(i,1)-x(i,4); y(i,1)-y(i,4); z(i,1)-z(i,4);].';
    H23(:,:,i)=1/(((x(i,2)-x(i,3))^2+(y(i,2)-y(i,3))^2+(z(i,2)-z(i,3))^2)^0.5)^2*[x(i,2)-x(i,3); y(i,2)-y(i,3); z(i,2)-z(i,3)]*[x(i,2)-x(i,3); y(i,2)-y(i,3); z(i,2)-z(i,3);].';
    H24(:,:,i)=1/(((x(i,2)-x(i,4))^2+(y(i,2)-y(i,4))^2+(z(i,2)-z(i,4))^2)^0.5)^2*[x(i,2)-x(i,4); y(i,2)-y(i,4); z(i,2)-z(i,4)]*[x(i,2)-x(i,4); y(i,2)-y(i,4); z(i,2)-z(i,4);].';
    H34(:,:,i)=1/(((x(i,3)-x(i,4))^2+(y(i,3)-y(i,4))^2+(z(i,3)-z(i,4))^2)^0.5)^2*[x(i,3)-x(i,4); y(i,3)-y(i,4); z(i,3)-z(i,4)]*[x(i,3)-x(i,4); y(i,3)-y(i,4); z(i,3)-z(i,4);].';
    
    KMSM(:,:,i)=[H12(:,:,i)+H13(:,:,i)+H14(:,:,i) -H12(:,:,i) -H13(:,:,i) -H14(:,:,i)
        -H12(:,:,i) H12(:,:,i)+H23(:,:,i)+H24(:,:,i) -H23(:,:,i) -H24(:,:,i)
        -H13(:,:,i) -H23(:,:,i) H13(:,:,i)+H23(:,:,i)+H34(:,:,i) -H34(:,:,i)
        -H14(:,:,i) -H24(:,:,i) -H34(:,:,i) H14(:,:,i)+H24(:,:,i)+H34(:,:,i)];
end


%KC�����߂�
for i=1:tetraNum(2)*tetraNum(1)/5;
    dr1(:,i)=1/(6*V0(i,1))*cross([x(i,2)-x(i,4); y(i,2)-y(i,4); z(i,2)-z(i,4)], [x(i,3)-x(i,4); y(i,3)-y(i,4); z(i,3)-z(i,4)]); %2�s1��
    dr2(:,i)=1/(6*V0(i,1))*cross([x(i,3)-x(i,1); y(i,3)-y(i,1); z(i,3)-z(i,1)], [x(i,4)-x(i,1); y(i,4)-y(i,1); z(i,4)-z(i,1)]);
    dr3(:,i)=1/(6*V0(i,1))*cross([x(i,4)-x(i,2); y(i,4)-y(i,2); z(i,4)-z(i,2)], [x(i,1)-x(i,2); y(i,1)-y(i,2); z(i,1)-z(i,2)]);
    dr4(:,i)=1/(6*V0(i,1))*cross([x(i,1)-x(i,3); y(i,1)-y(i,3); z(i,1)-z(i,3)], [x(i,2)-x(i,3); y(i,2)-y(i,3); z(i,2)-z(i,3)]);
    
    K11(:,:,i)=dr1(:,i)*dr1(:,i).';    K12(:,:,i)=dr1(:,i)*dr2(:,i).';    K13(:,:,i)=dr1(:,i)*dr3(:,i).';    K14(:,:,i)=dr1(:,i)*dr4(:,i).';
    K21(:,:,i)=dr2(:,i)*dr1(:,i).';    K22(:,:,i)=dr2(:,i)*dr2(:,i).';    K23(:,:,i)=dr2(:,i)*dr3(:,i).';    K24(:,:,i)=dr2(:,i)*dr4(:,i).';
    K31(:,:,i)=dr3(:,i)*dr1(:,i).';    K32(:,:,i)=dr3(:,i)*dr2(:,i).';    K33(:,:,i)=dr3(:,i)*dr3(:,i).';    K34(:,:,i)=dr3(:,i)*dr4(:,i).';
    K41(:,:,i)=dr4(:,i)*dr1(:,i).';    K42(:,:,i)=dr4(:,i)*dr2(:,i).';    K43(:,:,i)=dr4(:,i)*dr3(:,i).';    K44(:,:,i)=dr4(:,i)*dr4(:,i).';
    
    KC(:,:,i)=[K11(:,:,i) K12(:,:,i) K13(:,:,i) K14(:,:,i);K21(:,:,i) K22(:,:,i) K23(:,:,i) K24(:,:,i);K31(:,:,i) K32(:,:,i) K33(:,:,i) K34(:,:,i);K41(:,:,i) K42(:,:,i) K43(:,:,i) K44(:,:,i)];
    
end


%% �ŏ����@ �Q�l(http://www.eli.hokkai-s-u.ac.jp/~kikuchi/ma2/chap08.html)
% for i=1:tetraNum(2)*tetraNum(1)/5;
%     for j=1:12;
%         Y(12*(j-1)+1:12*j,1,i)=KCST(j,1:12,i);
%         X(12*(j-1)+1:12*j,1,i)=KMSM(j,1:12,i);   %k
%         X(12*(j-1)+1:12*j,2,i)=KC(j,1:12,i);   %ks
%     end
%     a(:,:,i)=X(:,:,i)\Y(:,:,i);
% end
% 
% ave(1,1)=mean(a(1,1,:)); %k
% ave(2,1)=mean(a(2,1,:)); %ks
% file_name_a = ['C:\Users\kou_0\OneDrive\�h�L�������g\����matlab\parameter\',muscleName, '_a.csv']
% % file_name_a = strcat("C:\Users\bubbl\Documents\shizuya_M1\DefMuscle_for_TUS_2018\matlab\surface2grid_0608\parameter\", muscleName, "_a.csv")
% csvwrite(file_name_a, ave);



