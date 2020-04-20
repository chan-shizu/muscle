clear
muscleNames = loadMuscleName();

% 各方向の最大頂点数 = Ni,Nj,Nk
%Ni=16,Nj=16を前提として、javaコードを書いている
Ni = 16;
Nj = 16;
Nk = 30;


for m =1
    
m 
%**********************
% 1.  load data
%**********************
%muscleName = ['L_',muscleNames{m}];
muscleName = [muscleNames{m}]

filename = ['traced\', muscleName, '.csv'];
[NUM, TXT, RAW] = xlsread(filename)

gridInfo = TXT(:,1:2);
XYZ = NUM(:, 4:end);
[nFibers, numq] = size(XYZ);%nFibersは手動で抽出したfiber数
%maxSlices = numq/3;

%**********************
% 2. スライス数をそろえる
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
    nSlices = maxq/3;

    fiberCoordinate = zeros(nSlices, 3);%3=xyz
    for iSlice = 1:nSlices
        fiberCoordinate(iSlice, 1) = XYZ(iFiber, 3*(iSlice-1)+1);
        fiberCoordinate(iSlice, 2) = XYZ(iFiber, 3*(iSlice-1)+2);
        fiberCoordinate(iSlice, 3) = XYZ(iFiber, 3*(iSlice-1)+3);
        %fiberCoordinate(iSlice, 3) = XYZ(iFiber, 3*(iSlice-1)+3)+0.002;%上腕筋の下部だけは少し前に移動した方が安全
    end

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
    lengthCoordinate = calcLengthCoordinate(fiberCoordinate)
    totalL = lengthCoordinate(end);
    newLengthCoordinate = [0:  totalL/(Nk-1)  :totalL]';%(Nk-1)個のセグメントに分割
    %size(newLengthCoordinate) must = Nk

    %yy = splineAllColumn(x,y,xx)
    newFiberCoordinate = splineAllColumn(lengthCoordinate, fiberCoordinate, newLengthCoordinate);

    yy = newFiberCoordinate;
    plot3(yy(:,1), yy(:,2), yy(:,3), 'b.' );hold on
    axis equal;

    fiber{iFiber}.coordinate = newFiberCoordinate;

    %i
    str = gridInfo{iFiber,1};
    k = strfind(str, '//');
    a = str(1:k-1);
    if  strcmp(a, '0')
        fiber{iFiber}.i = '0';
    else
        b = str(k+2:end);
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
% 3. Ni=16, Nj=7の表面ファイバーを作成
%****************************************************
clear k

for k=1:Nk

    for side =1:4
        %一辺目（j=0でiが増加する方向）
        IJXYZ = [];
        for iFiber = 1:nFibers

            if side == 1 %west  J=0の辺
                index = str2num(fiber{iFiber}.j);
                num01 = 0;%str 0 or 1
                sort12 = 1;% 1(i) or 2(j)でソート
            elseif side == 2 %south  i=0の辺
                index = str2num(fiber{iFiber}.i);
                num01 = 0;
                sort12 = 2;
            elseif side == 3 %east  J=1の辺
                index = str2num(fiber{iFiber}.j);
                num01 = 1;
                sort12 = 1;
            elseif side == 4 %north  i=1の辺
                index = str2num(fiber{iFiber}.i);
                num01 = 1;
                sort12 = 2;
            end

            %該当のfiberがあったら、それを蓄積していく
            if (index == num01)
                tmp = [str2num(fiber{iFiber}.i) ,... //数字のi
                    str2num(fiber{iFiber}.j),...  //数字のj
                    fiber{iFiber}.coordinate(k,1:3)  ];%kスライス目のxyz座標
                IJXYZ = [IJXYZ; tmp];
                % 例）IJXYZ =
                %          0            0    0.1871    1.1710    0.0719
                %     1.0000         0    0.1861    1.1852    0.0674
                %     0.5000         0    0.1861    1.1757    0.0684
            end
        end

        IJXYZ
        [Y, I] = sort(IJXYZ(:,sort12));%iでソート
        IJXYZ = IJXYZ(I,:)
        %     IJXYZ =
        %          0            0    0.1871    1.1710    0.0719
        %     0.5000         0    0.1861    1.1757    0.0684
        %     1.0000         0    0.1861    1.1852    0.0674
        coordinate = IJXYZ(:, 3:5);
        lengthCoordinate = calcLengthCoordinate(coordinate);
        totalL = lengthCoordinate(end);
        if (side==1 || side==3)
            newLengthCoordinate = [0:  totalL/(Ni-1)  :totalL]';%(Ni-1)個のセグメントに分割
        else
            newLengthCoordinate = [0:  totalL/(Nj-1)  :totalL]';%(Nj-1)個のセグメントに分割
        end

        newIJXYZ{side, k} = splineAllColumn(lengthCoordinate, IJXYZ, newLengthCoordinate)
        size(newIJXYZ{side, k})

    end % for side =1:4

end %for k=1:Nk


%************************************************
%newIJXYZ{side, k}をcsvファイルに保存
%************************************************
sfilename = ['aligned\', muscleName,  '_aligned.csv'];
fid = fopen(sfilename,'wt'); % 書き込み用にファイルオープン
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
            fprintf(fid, '%2.9f,',  newIJXYZ{side, k}(n,3) );%x座標
            fprintf(fid, '%2.9f,',  newIJXYZ{side, k}(n,4) );%y座標
            fprintf(fid, '%2.9f,',  newIJXYZ{side, k}(n,5) );%z座標
        end
        fprintf(fid, '\n');

    end %n = 1:N

end

fclose(fid);
display(muscleName)


end





