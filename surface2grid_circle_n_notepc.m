clear
muscleNames = loadMuscleName();

% �e�����̍ő咸�_�� = Ni,Nj,Nk
%Ni=16,Nj=16��O��Ƃ��āAjava�R�[�h�������Ă���
%�����������߂�
Ni = 5;
Nj = 5;
Nk = 6;


for m =1
    
m 
%**********************
% 1.  load data
%**********************
%muscleName = ['L_',muscleNames{m}];
muscleName = [muscleNames{m}]

%filename = ['traced\', muscleName, '.csv'];%defmuscle����o�͂��ꂽcsv����ёւ�������
% filename = "C:\Users\kou_0\OneDrive\�h�L�������g\����matlab\traced\L_Digastric_Posterior.csv"
filename = ['C:\Users\kou_0\OneDrive\�h�L�������g\����matlab\traced\', muscleName, '.csv'];
% filename = ['C:\Users\bubbl\Documents\shizuya_M1\DefMuscle_for_TUS_2018\matlab\surface2grid_0608\traced\', muscleName, '.csv'];
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
% sfilename = ['C:\Users\bubbl\Documents\shizuya_M1\DefMuscle_for_TUS_2018\matlab\surface2grid_0608\aligned\', muscleName,  '_aligned.csv'];
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





