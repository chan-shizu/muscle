function SwalInverseDynamic = SwalInverseDynamic();
% 粘弾性により逆動力学
clc;
clear;
warning('on','all');
warning;
weakenMuscle = 1;     % if weaken muscle then 1, else 0;

% input Data
Kht = 20;
Kct = 20;
Kx = 3;
Ky = 3;
Kzth = 0.001;
Cx = 1;
Cy = 1.0;
Czth = 0.00001;
%-----------------------------
csvwrite('basis.csv',[Kht;Kct;Kx;Ky;Kzth;Cx;Cy;Czth]);
Kxyz = [Kx Ky Kzth];
Cxyz = [Cx Cy Czth];

%% データの読み込み
Mscl = csvread('MuscleParam-0110.csv', 2, 1);             % 筋の始点・終点の座標，断面積(10列24行)，付着点情報
DynData = csvread('DynamicsData-100.csv', 2, 0);      % 舌骨の軌跡
head = csvread('cranium-100.csv',0,0);
Lce0 = csvread('Lce0-mix-0110.csv',0,0);
HTspring = csvread('spring coord-hyoid&thyroid.csv', 25, 1);
HTSNum = size(HTspring);
CTspring = csvread('spring coord-tracheal&cricoid.csv', 17, 1);
CTSNum = size(CTspring);
MsclNum = size(Mscl);                               % 筋の本数
DataNum = size(DynData);                       % データ数
BodyNum = 2;	dt = 1.60/(DataNum(1)-1);      % 時間間隔
% preparation of hill equation
Fmax = 0.7*10^6*Mscl(1:MsclNum(1), 10)*0.7;      % 最大筋力
Fmaxa=0.7*10^6*Mscl(1:MsclNum(1),10)*20*0.7;       %若年化させた際の筋長
Fmin(1:24) = 0;                        % 最小筋力（全て0）
Vsh = 0.3;
Vshl = 0.23;
Vml = 1.3;
Ver = 0.5;
alpha = 1;              % aは筋の活性化レベルを表す変数
PEsh = 4;
PExm = 0.34;
MsclName = table2cell(readtable('MuscleParam-0103.csv'));

syms xv xi yo yv yi zo zv AA BB CC z
eqn = [(yo-yv)/(zo-zv)==(yv-yi)/(xv-xi), xv == AA*yv^2+BB* yv+CC];
solx = solve(eqn,[xv yv]);
Func = (char(solx.yv(1)));

% 咽頭縫線の近似線を求める
Raphe = polyfit([Mscl(1:16,8);Mscl(37:44,8)],[Mscl(1:16,7);Mscl(37:44,7)],1);

% 最適化計算の準備
H = 2*diag(Mscl(1:MsclNum(1), (BodyNum+1)*3+1))^-2;  % 2次項（対称行列）
fl(1:MsclNum(1)) = 0;                           % 線形項は0
Fmin(1:MsclNum(1)) = 0;               % 最小筋力 （全て0）
errortime= [];
%% 計算の準備
% 舌骨と甲状軟骨を分別
for i = 1:BodyNum;
    TraData(1:DataNum(1),1,i) = DynData(1:DataNum(1),1);       TraData(1:DataNum(1),2,i) = DynData(1:DataNum(1),(i-1)*6+2);      TraData(1:DataNum(1),3,i) = DynData(1:DataNum(1),(i-1)*6+3);
end
% バネの行列を作る
Spring(1:HTSNum(1)+CTSNum,1:9) = 0;
Spring(1:HTSNum(1),1:HTSNum(2)) = HTspring;
Spring(HTSNum(1)+1:HTSNum(1)+CTSNum(1),1+3:3+CTSNum(2)) = CTspring;
SNum = size(Spring);

% ボディの切り替え位置を探索
for i = 1:BodyNum;            % ボディの切り替え位置を探索
    BodyTrans(i) = 1;
    for j = 1:MsclNum(1);
        if Mscl(j, (BodyNum*3+3+1)+i) == 1;
            BodyTrans(i) = j;
            break;
        end;
    end;
end;
% バネの切り替え位置を入力
SprgTrans = [1, HTSNum(1)+1];
SprgTrans2 = [HTSNum(1) , SNum(1)];   % the last spring position
K(1:HTSNum(1)) = Kht;
K(HTSNum(1)+1:SNum(1)) = Kct;

%% Hill'sの計算、逆動力学、最適化計算
for i = 1:DataNum(1);
    % 膜バネの座標変換
    STrans(1:SNum(1), 1:SNum(2)) = 0;       % 並進の変換行列(0)
    STrans(1:HTSNum(1), 1) = TraData(i,2,1);  % 並進の変換行列(1)
    STrans(1:HTSNum(1), 2) = TraData(i,3,1);  % 並進の変換行列(2)
    STrans(:, 4) = TraData(i,2,2);             % 並進の変換行列(4)
    STrans(:, 5) = TraData(i,3,2);             % 並進の変換行列(5)
    SRotMat(:,:) = eye(SNum(2));
    for k = 1:BodyNum;
        Xrot = 0;
        Yrot = 0;
        Zrot = TraData(i,4,k);	SRotMat(3*(k-1)+1, 3*(k-1)+1) = cos(Zrot)*cos(Yrot)
        SRotMat(3*(k-1)+1, 3*(k-1)+2) = sin(Zrot)*cos(Yrot);
        SRotMat(3*(k-1)+1, 3*(k-1)+3) = -sin(Yrot);
        SRotMat(3*(k-1)+2,3*(k-1)+1)=cos(Zrot)*sin(Yrot)*sin(Xrot) - sin(Zrot)*cos(Xrot);
        SRotMat(3*(k-1)+2,3*(k-1)+2)= sin(Zrot)*sin(Yrot)*sin(Xrot)+cos(Zrot)*cos(Xrot);
        SRotMat(3*(k-1)+2,3*(k-1)+3)=cos(Yrot)*sin(Xrot);
        SRotMat(3*(k-1)+3, 3*(k-1)+1) = cos(Zrot)*sin(Yrot)*cos(Xrot) + sin(Zrot)*sin(Xrot);
        SRotMat(3*(k-1)+3, 3*(k-1)+2) = sin(Zrot)*sin(Yrot)*cos(Xrot) - cos(Zrot)*sin(Xrot);
        SRotMat(3*(k-1)+3, 3*(k-1)+3) = cos(Yrot)*cos(Xrot);
    end;
    tmpGSpring(:,:)=Spring*SRotMat+STrans(:,:);        %ローカル座標をグローバル座標に変換
    for k = 1:BodyNum;
        % 付着位置毎にバネを分類
        GSpring(1:SNum(1),1:6,k) = 0;
        GSpring(1:SprgTrans2(k),1:3,k) = tmpGSpring(1:SprgTrans2(k),(k-1)*3+1:(k-1)*3+3);
        GSpring(1:SprgTrans2(k),4:6,k)=tmpGSpring(1:SprgTrans2(k),(k)*3+1:(k)*3+3);
        if k == 2;
            GSpring(1:HTSNum(1),4:6,k) = tmpGSpring(1:HTSNum(1),1:3);
        end;
    end;
    for k = 1:BodyNum;
        % バネの角度を計算
        for j = 1:SNum(1);
            SxyAngle(j,1,k)=atan((GSpring(j,5,k)-(GSpring(j,2,k)))/(GSpring(j,4,k)-(GSpring(j,1,k))))
            if GSpring(j,4,k)-(GSpring(j,1,k)) < 0;    % (x,y) = (-, +) or (-, -)
                SxyAngle(j,1,k) = SxyAngle(j,1,k) + pi();
            elseif GSpring(j,5,k)-(GSpring(j,2,k)) < 0;         % (x,y) = (+, -)
                SxyAngle(j,1,k) = SxyAngle(j,1,k) + 2*pi();
            end;
            if GSpring(j,5,k)-(GSpring(j,2,k)) == 0 && (GSpring(j,4,k)-(GSpring(j,1,k))) == 0;
                SxyAngle(j,1,k) = 0;
            end;
            SzAngle(j,1,k) = acos((GSpring(j,6,k)-GSpring(j,3,k))./sqrt((GSpring(j,4,k)-GSpring(j,1,k)).^2 + (GSpring(j,5,k)-GSpring(j,2,k)).^2+(GSpring(j,6,k)-GSpring(j,3,k)).^2));
            if GSpring(j,1:6,k) == 0;
                SzAngle(j,1,k) = 0;
            end;
            
            % spring vector
            SXdir(i,j,k) = cos(SxyAngle(j,1,k)) * sin(SzAngle(j,1,k));   %筋のx方向成分
            SYdir(i,j,k) = sin(SxyAngle(j,1,k)) * sin(SzAngle(j,1,k));    %筋のy方向成分
            SZdir(i,j,k) = cos(SzAngle(j,1,k));
            if GSpring(j,5,k)-(GSpring(j,2,k)) == 0 && (GSpring(j,4,k)-(GSpring(j,1,k))) == 0;
                SXdir(j,1,k) = 0;
            end;
            
            SprgVec = [SXdir(i,j,k) SYdir(i,j,k) SZdir(i,j,k)];
            SHydVec = [(GSpring(j,1,k) - TraData(i,2,k)) (GSpring(j,2,k)-TraData(i,3,k)) (GSpring(j,3,k)-0)];     % 剛体上の筋肉の付着点の向き
            SCrossVec = cross(SHydVec, SprgVec);      % 外積
            SMomArm(i,j,k) = SCrossVec(3);
            
            %　バネの長さと伸びの計算
            Length(i,j,k) = sqrt((GSpring(j,4,k)-GSpring(j,1,k)).^2 + (GSpring(j,5,k)-GSpring(j,2,k)).^2+(GSpring(j,6,k)-GSpring(j,3,k)).^2);
            if i == 1;
                Lengthen(1,1:SNum(1),k) = 0;
            else
                Lengthen(i,j,k) = Length(i,j,k) - Length(1,j,k);	end;
            if Lengthen(i,j,k) < 0;
                Lengthen(i ,j ,k) = 0;
            end;
        end;
    end;
    %     cla;
    tmp_Mscl = HeadMov(Mscl,head(i,:),MsclNum);     %頭の回転と変位の計算
    %筋の座標変換-------------------------------------------------------------
    for k = 1:BodyNum;
        mTrans(MsclNum(1), MsclNum(2),k) = 0;   % 並進の変換行列(0)
        mTrans(BodyTrans(k):MsclNum(1), (k-1)*3+1 ,k) = TraData(i,2,k);           % 並進の変換行列(1)
        mTrans(BodyTrans(k):MsclNum(1), (k-1)*3+2 ,k) = TraData(i,3,k);           % 並進の変換行列(2)
        
        mRotMat(:,:,k) = eye(MsclNum(2));
        Xrot = 0;
        Yrot = 0;
        Zrot = TraData(i,4,k);
        mRotMat(3*(k-1)+1, 3*(k-1)+1,k) = cos(Zrot)*cos(Yrot);  % 回転行列の(1, 1)
        mRotMat(3*(k-1)+1, 3*(k-1)+2,k) = sin(Zrot)*cos(Yrot);  % 回転行列の(1, 2)
        mRotMat(3*(k-1)+1, 3*(k-1)+3,k) = -sin(Yrot);  % 回転行列の(1, 3)
        mRotMat(3*(k-1)+2, 3*(k-1)+1,k) = cos(Zrot)*sin(Yrot)*sin(Xrot) - sin(Zrot)*cos(Xrot); % 回転行列の(2, 1)
        mRotMat(3*(k-1)+2, 3*(k-1)+2,k) = sin(Zrot)*sin(Yrot)*sin(Xrot) + cos(Zrot)*cos(Xrot);	% 回転行列の(2, 2)
        mRotMat(3*(k-1)+2, 3*(k-1)+3,k) = cos(Yrot)*sin(Xrot);  % 回転行列の(2, 3)
        mRotMat(3*(k-1)+3, 3*(k-1)+1,k) = cos(Zrot)*sin(Yrot)*cos(Xrot) + sin(Zrot)*sin(Xrot);  % 回転行列の(3, 1)
        mRotMat(3*(k-1)+3, 3*(k-1)+2,k) = sin(Zrot)*sin(Yrot)*cos(Xrot) - cos(Zrot)*sin(Xrot);  % 回転行列の(3, 2)
        mRotMat(3*(k-1)+3, 3*(k-1)+3,k) = cos(Yrot)*cos(Xrot);  % 回転行列の(3, 3      
        tmpGMscl(:,:,k) = tmp_Mscl*mRotMat(1:MsclNum(2),1:MsclNum(2),k)+mTrans(:,:,k);                      % ローカル座標をグローバル座標に変換
    end;
    
    % 付着位置毎に筋を分類
    GMscl(1:MsclNum(1), 1:6, 1:2) = 0;
    for k = 1:BodyNum;
        for j = 1:MsclNum(1);
            if tmpGMscl(j,11,k) == 1 && tmpGMscl(j,12,k) == 1;
                if k == 1;
                    GMscl(j,1:3,k) = tmpGMscl(j,(k-1)*3+1:(k-1)*3+3,k);
                    GMscl(j,4:6,k) = tmpGMscl(j,4:6,k+1);
                else
                    GMscl(j,1:3,k) = tmpGMscl(j,(k-1)*3+1:(k-1)*3+3,k);
                    GMscl(j,4:6,k) = tmpGMscl(j,1:3,1);
                end;
            elseif tmpGMscl(j,k+10,k) == 1;
                GMscl(j,1:3,k) = tmpGMscl(j,(k-1)*3+1:(k-1)*3+3,k);
                GMscl(j,4:6,k) = tmpGMscl(j,7:9,k);
            end;
        end;
    end;
    
    % 筋の長さを計算
    Lce(i,:) = MLength(GMscl(:,1:3,:),GMscl(:,4:6,:),MsclNum);
    %------- PG bend------------------------
    PGTransit(1:8,3) = GMscl(1:8,3,1);
    for j = 1:8;
        a = Raphe(1);
        b = Raphe(2) + (TraData(i,2,1)-TraData(1,2,1));
        c = GMscl(j,2,1);
        d = GMscl(j,1,1);
        e = GMscl(j,5,1);	f = PGTransit(j,3);
        g = GMscl(j,6,1);
        PGTransit(j,2) = (-(b-d-a*e-abs(f-g))-sqrt((b-d-a*e-abs(f-g))^2-4*a*(c*abs(f-g)-(b-d)*e)))/(2*a);
        PGTransit(j,1) = PGTransit(j,2) * a + b;
    end;
    PG_length1 = sqrt((GMscl(1:8,1,1)-PGTransit(:,1)).^2 + (GMscl(1:8,2,1)-PGTransit(:,2)).^2 + (GMscl(1:8,3,1)-PGTransit(:,3)).^2);
    PG_length2 = sqrt((GMscl(1:8,4,1)-PGTransit(:,1)).^2 + (GMscl(1:8,5,1)-PGTransit(:,2)).^2 + (GMscl(1:8,6,1)-PGTransit(:,3)).^2);
    GMscl(1:8,4:6,1) = PGTransit;
    Lce(i,1:8) = PG_length1 + PG_length2;
    %------- PL bend------------------------
    PLTransit(1:8,3) = GMscl(1:8,3,1);
    for j = 1:8;
        a = Raphe(1);
        b = Raphe(2) + (TraData(i,2,1)-TraData(1,2,1));
        c = GMscl(8+j,2,1);
        d = GMscl(8+j,1,1);
        e = GMscl(8+j,5,1);
        f = PLTransit(j,3);
        g = GMscl(8+j,6,1);
        PLTransit(j,2) = (-(b-d-a*e-abs(f-g))-sqrt((b-d-a*e-abs(f-g))^2-4*a*(c*abs(f-g)-(b-d)*e)))/(2*a);
        PLTransit(j,1) = PLTransit(j,2) * a + b;
    end;
    PL_length1 = sqrt((GMscl(9:16,1,1)-PLTransit(:,1)).^2 + (GMscl(9:16,2,1)-PLTransit(:,2)).^2 + (GMscl(9:16,3,1)-PLTransit(:,3)).^2);
    PL_length2 = sqrt((GMscl(9:16,4,1)-PLTransit(:,1)).^2 + (GMscl(9:16,5,1)-PLTransit(:,2)).^2 + (GMscl(9:16,6,1)-PLTransit(:,3)).^2);
    GMscl(9:16,4:6,1) = PLTransit;
    Lce(i,9:16) = PL_length1 + PL_length2;
    
    %----------PT bend---------------------
    PTTransit(1:8,3) = GMscl(1:8,3,1);
    for j = 1:8;
        a = Raphe(1);
        b = Raphe(2) + (TraData(i,2,1)-TraData(1,2,1));
        c = GMscl(36+j,2,2);
        d = GMscl(36+j,1,2);
        e = GMscl(36+j,5,2);
        f = PTTransit(j,3);
        g = GMscl(36+j,6,2);
        PTTransit(j,2) = (-(b-d-a*e-abs(f-g))-sqrt((b-d-a*e-abs(f-g))^2-4*a*(c*abs(f-g)-(b-d)*e)))/(2*a);
        PTTransit(j,1) = PTTransit(j,2) * a + b;
    end;
    PT_Length1 =  sqrt((GMscl(37:44,1,2)-PTTransit(:,1)).^2 + (GMscl(37:44,2,2)-PTTransit(:,2)).^2 + (GMscl(37:44,3,2)-PTTransit(:,3)).^2);
    PT_Length2 = sqrt((GMscl(37:44,4,2)-PTTransit(:,1)).^2 + (GMscl(37:44,5,2)-PTTransit(:,2)).^2 + (GMscl(37:44,6,2)-PTTransit(:,3)).^2);
    GMscl(37:44,4:6,2) = PTTransit;
    Lce(i,37:44) = PT_Length1 + PT_Length2;
    
    for k = 1:BodyNum;
        for j = 1:MsclNum(1);
            mxyAngle(j,1,k) = atan((GMscl(j,5,k)-GMscl(j,2,k))/(GMscl(j,4,k)-GMscl(j,1,k)));        %角度を計算
            if (GMscl(j,5,k)-GMscl(j,2,k)) == 0  &&  (GMscl(j,4,k)-GMscl(j,1,k)) == 0;
                mxyAngle(j,1,k) = 0;
            end;
            if (GMscl(j,4,k)-GMscl(j,1,k)) < 0;   % (x,y) = (-, +) or (-, -)
                mxyAngle(j,1,k) = mxyAngle(j,1,k) + pi();
            elseif (GMscl(j,5,k)-GMscl(j,2,k)) < 0;         % (x,y) = (+, -)
                mxyAngle(j,1,k) = mxyAngle(j,1,k) + 2*pi();
            end;
            mzAngle(j,1,k) = acos((GMscl(j,6,k)	GMscl(j,3,k))/(sqrt((GMscl(j,4,k)-GMscl(j,1,k))^2 + (GMscl(j,5,k)-GMscl(j,2,k))^2 + (GMscl(j,6,k)-GMscl(j,3,k))^2)));
            if GMscl(j,:,k) == 0;
                mzAngle(j,1,k) = 0;
            end;
            
            % muscle vector
            mXdir(i,j,k) = cos(mxyAngle(j,1,k)) * sin(mzAngle(j,1,k));      % 筋のx方向成分
            mYdir(i,j,k) = sin(mxyAngle(j,1,k)) * sin(mzAngle(j,1,k));      % 筋のy方向成分
            mZdir(i,j,k) = cos(mzAngle(j,1,k));
            mXdir(i,BodyTrans(2)+2:MsclNum(1),1) = 0;
            mXdir(i,1:BodyTrans(2)-1,2) = 0;
            mZdir(i,BodyTrans(2)+2:MsclNum(1),1) = 0;
            mZdir(i,1:BodyTrans(2)-1,2) = 0;
            mMsclVec = [mXdir(i,j,k) mYdir(i,j,k) mZdir(i,j,k)];
            mHydVec = [(GMscl(j,1,k) - TraData(i,2,k)) (GMscl(j,2,k)-TraData(i,3,k)) (GMscl(j,3,k))];   % 剛体上の筋肉の付着点の向き
            mCrossVec = cross(mHydVec, mMsclVec);          % 外積
            mxMomArm(i,j,k) = mCrossVec(1);
            myMomArm(i,j,k) = mCrossVec(2);
            mzMomArm(i,j,k) = mCrossVec(3);
        end;
    end;
    
    %% Hill'sの計算
    if i == 1;
        Vce(i,1:MsclNum(1)) = 0;
    else;
        Vce(i,:) = -(Lce(i,:) - Lce(i-1,:))/dt;
    end;
    
    Lmus0 = Lce0;
    Lcesh = 0.5 * Lce0;
    f_Lce(i,:) = exp(-( (Lce(i,:) - Lce0)./Lcesh).^2);
    % calculation of f(Vce)
    Vvm = 6 * Lce0;
    Vmax(i,:) = Vvm .* ( 1 - Ver .* (1 - alpha * f_Lce(i,:)));
    for j = 1:MsclNum(1)
        if Vce(i,j) <= -Vmax(i,j);
            f_Vce(i,j) = 0;
        elseif Vce(i,j) > -Vmax(i,j)  &&  Vce(i,j) < 0;
            f_Vce(i,j) = (Vsh * Vmax(i,j) + Vsh * Vce(i,j)) / (Vsh * Vmax(i,j) - Vce(i,j));
        elseif Vce(i,j) >= 0;
            f_Vce(i,j) = (Vsh * Vshl * Vmax(i,j) + Vml * Vce(i,j)) / (Vsh * Vshl * Vmax(i,j) + Vce(i,j));
        end;
    end;
    HillPassive(i,:) = Fmax'./exp(PEsh) .* (exp((Lce(i,:)-Lmus0) ./ (Lmus0) .* PEsh/PExm) -1);
    if Lce(i,:)<= Lmus0;
        HillPassive(i,j) =0;
    end;
    HillActive(i,:) = alpha .* f_Lce(i,:) .* f_Vce(i,:) .* Fmaxa';
    Lpercentage(i,:) = Lce(i,:) ./ Lce0;
    %% 逆動力学計算
    for k = 1:BodyNum;
        if i == 1;
            SpringForce(1,1:SNum(1),1:2) = 0;
        else;
            SpringForce(i,:,k) = K.*Lengthen(i,:,k);
        end
        
        % 舌骨と甲状軟骨の変位と速度とそのpassive力
        Displacement(i,1:3,k) = TraData(i,2:4,k) - TraData(1,2:4,k);
        if i > 1;
            Bvel(i,1:3,k) = (TraData(i,2:4,k) - TraData(i-1,2:4,k)) / dt;
        else
            Bvel(i,1:3,k) =  0;
        end;
        Fpb(i,1:3,k) = Displacement(i,1:3,k) .* Kxyz + Bvel(i,1:3,k) .* Cxyz;  	% Spring passive力
        Fx(i,1,k) = SpringForce(i, 1:SNum(1),k)*SXdir(i,:,k)';      %x方向のpassive力
        Fy(i,1,k) = SpringForce(i, 1:SNum(1),k)*SYdir(i,:,k)';           %　y方向のpassive力
        Tr(i,1,k) = SpringForce(i, 1:SNum(1),k)*SMomArm(i,:,k)';   %回転のpassive力
        
        Fhill(i,1,k) = HillPassive(i, 1:MsclNum(1))*mXdir(i,:,k)';   %x方向のpassive力
        Fhill(i,2,k) = HillPassive(i, 1:MsclNum(1))*mYdir(i,:,k)';   %y方向のpassive力
        Fhill(i,3,k) = HillPassive(i, 1:MsclNum(1))*mzMomArm(i,:,k)';    %回転のpassive力
        
        HydF0(i,1,k) = DynData(i,(k-1)*3+14);  %挙動から求めた力F
        HydF0(i,2,k) = DynData(i,(k-1)*3+15);
        HydF0(i,3,k) = DynData(i,(k-1)*3+16);
        
        tmpHydF(i,1,k) = DynData(i,(k-1)*3+14) - Fhill(i,1,k);
        tmpHydF(i,2,k) = DynData(i,(k-1)*3+15) - Fhill(i,2,k);
        tmpHydF(i,3,k) = DynData(i,(k-1)*3+16) - Fhill(i,3,k);
        
        HydF(i,1,k) = tmpHydF(i,1,k) - Fx(i,1,k) + Fpb(i,1,k);
        HydF(i,2,k) = tmpHydF(i,2,k) - Fy(i,1,k) + Fpb(i,2,k);
        HydF(i,3,k) = tmpHydF(i,3,k) - Tr(i,1,k) + Fpb(i,3,k);
    end;
    %% 最適化計算
    Aeq(1:BodyNum*6, 1:MsclNum(1)) = 0;
    for j = 1:BodyNum;
        Aeq(j*6-5, :) = mXdir(i,:, j)';
        Aeq(j*6-4, :) = mYdir(i,:, j)';
        Aeq(j*6-3, :) = mZdir(i,:, j)';
        Aeq(j*6-2, :) = mxMomArm(i,:, j);
        Aeq(j*6-1, :) = myMomArm(i,:, j);
        Aeq(j*6, :) = mzMomArm(i,:, j);
    end;
    
    beq(1:BodyNum*6) = [HydF(i,1,1) HydF(i,2,1) 0 0 0 HydF(i,3,1) HydF(i,1,2) HydF(i,2,2) 0 0 0 HydF(i,3,2)];
    
    % 二次計画法
    [x, fval, exitflag, output, lambda] = quadprog(H,fl, [], [], Aeq, beq, Fmin, HillActive(i,:));
    MForce(i, 1) = (i-1) * (dt);
    MForce(i, 2:MsclNum(1)+1) = x;
    if any([-2 -1 0]  == exitflag);
        errortime = [errortime
            i exitflag];
    end;
    
    xydir(i,:) = mXdir(i,:,1);
    xydir(i,35:46) = mXdir(i,35:46,2);
    xydir(i,[2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34]) = mYdir(i,[2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34],1);
    xydir(i,[36 38 40 42 44 46]) = mYdir(i,[36 38 40 42 44 46],2);
    Fxy(i,:) = x'.*xydir(i,:);
    Fpexy(i,:) = HillPassive(i,:).*xydir(i,:);
    disp(i);
end;
csvwrite('SwalOptResult.csv', MForce);
end

function Ans = HeadMov(Mscl,x,Num)
RotMat(:,:) = eye(Num(2));
RotMat(7,7) = cos(x(7)*pi()/180);
RotMat(7,8) = sin(x(7)*pi()/180);
RotMat(8,7) = -sin(x(7)*pi()/180);
RotMat(8,8) = cos(x(7)*pi()/180);       % 回転の変換行列
Head = Mscl*RotMat;
Ans = Mscl;
for j=1:Num(1);
    if Mscl(j,13) == 1 || Mscl(j,13) == 2;
        Ans(j,7) = Head(j,7) + x(2);
        Ans(j,8) = Head(j,8) + x(3);	  end;
end;
end

% 咽頭縫線の移動
function Ans = RapheTrans(Mscl,s,Num)
formula = RapheFormula(Mscl,s);
Ans = Mscl;
for j=1:Num(1);
    if Mscl(j,13) == 4;
        Ans(j,7) = polyval(formula, Ans(j,8));
    end;
end;
end

% 咽頭縫線の方程式
function Ans = RapheFormula(Mscl,s);
for l = 1:4;
    Raphe(l,:) = s(l*2:l*2+1);
end;
Ans = polyfix(Raphe(:,2),Raphe(:,1),2,[s(3)],[s(2)]);
end

%　筋の長さ
function Ans = MLength(P1,P2,n);
for  j = 1:n;
    Ans(j) = sqrt((P1(j,1,1)-P2(j,1,1))^2+(P1(j,2,1)-P2(j,2,1))^2+(P1(j,3,1)-P2(j,3,1))^2);
    if Ans(j) == 0;
        Ans(j) = sqrt((P1(j,1,2)-P2(j,1,2))^2+(P1(j,2,2)-P2(j,2,2))^2+(P1(j,3,2)-P2(j,3,2))^2);
    end;
end;
end % 咽頭筋を曲げる

function Ans = ViaPoint(PG1,P,raphe,displacement,solx,Func)
AA = raphe(1);
BB = raphe(2);
CC = raphe(3);
xi = P(1);
yi = P(2);
yo = P(5);
zo = P(6);
if P(3) >0;
    Ans(3) = PG1(1,3);
else;
    Ans(3) = PG1(2,3);
end;
zv = (Ans(3));

YV = vpa(subs(solx.yv(1)));
if (YV(1) <=P(2) && YV(1) >=P(5)) || (YV(1) >=P(2) && YV(1) <=P(5)) && isreal(YV(1)) == 1;
    Ans(2) = YV(1);
else;
    zv = -zv;
    zo = -zo;
    YV = vpa(subs(solx.yv));
    if (YV(1) <=P(2) && YV(1) >=P(5)) || (YV(1) >=P(2) && YV(1) <=P(5)) && isreal(YV(1)) == 1;
        Ans(2) = YV(1);
    else;
        disp('error_YV');
    end;
end;
XV = AA*Ans(2)^2 + BB*Ans(2) + CC;
Ans(1) = XV;
Ans(4) = sqrt((Ans(1)-P(1))^2+(Ans(2)-P(2))^2+(Ans(3)-P(3))^2) +  sqrt((Ans(1)-P(4))^2+(Ans(2)-P(5))^2+(Ans(3)-P(6))^2);
end
