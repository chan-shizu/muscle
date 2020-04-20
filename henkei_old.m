%function area_mive=2Darea():

clc;
clear;
warning('on','all');
warning;

%入力
% a=csvread('a_r.csv', 0, 0); 
% springk=a(1,1)/10000;    %1900
springk=100;  %100

% kv=a(2,1)*50000;   %85
kv=0.25;    %0.25
muscleNames = loadMuscleName();
muscle_name = [muscleNames{1}]

y=5;
t=5;
h=6;
mass=5*10^(-3);
conk=60;%60
conc=60;
bNum=(y-1)*(t-1)*(h-1)/8;   %blockNumber

%-9.98031度回転させる
file_name_data = strcat("C:\Users\bubbl\Documents\shizuya_M1\DefMuscle_for_TUS_2018\matlab\surface2grid_0608\muscle\min\urabe_arm_data.csv")%"data_rot_h.csv"
file_name_data0 = strcat("C:\Users\bubbl\Documents\shizuya_M1\DefMuscle_for_TUS_2018\matlab\surface2grid_0608\muscle\min\urabe_arm_data0.csv")%"data0_rot_h.csv"
file_name_se = strcat("C:\Users\bubbl\Documents\shizuya_M1\DefMuscle_for_TUS_2018\matlab\surface2grid_0608\parameter\", muscle_name, "_se.csv")
file_name_tetra = strcat("C:\Users\bubbl\Documents\shizuya_M1\DefMuscle_for_TUS_2018\matlab\surface2grid_0608\parameter\", muscle_name, "_tetra.csv")
data= csvread(file_name_data, 0, 0)/1000;
data0=csvread(file_name_data0, 0, 0)/1000;
se=csvread(file_name_se, 0, 0);
tetra=csvread(file_name_tetra, 0, 0);

% data= csvread('data_rot_h.csv', 0, 0)/1000; 
% data0=csvread('data0_rot_h.csv', 0, 0)/1000; 
% tetra=csvread('tetra_arm.csv', 0, 0);
% se=csvread('se_arm.csv', 0, 0);

timeNum = size(data); 
pointNum = size(data0); 
seNum = size(se); 
tetraNum = size(tetra); 
dt=1*10^(-3); 




% Hillモデル計算の準備
Fmax = 2940/y/t/(h-1);                 % 筋線維1本当たりの最大筋力(筋肉ごとに算出)
Vsh = 0.3;
Vshl = 0.23;         
Vml = 1.3;        
Ver = 0.5;
alpha = 1*0.1;                                              % 筋の最大活性化レベルを表す変数
PEsh = 10;                  %筋肉の部位ごとの特性値
PExm = 0.4;                 %筋肉の部位ごとの特性値

tic
for i=1:timeNum(1)
if i==1
datai_x(i,1:pointNum)=data0(:,1).';  %i時における各点の座標
datai_y(i,1:pointNum)=data0(:,2).';
datai_z(i,1:pointNum)=data0(:,3).';
else
     %固定面
    datai_x(i,1:y*t)=data0(1:y*t,1);   
    datai_y(i,1:y*t)=data0(1:y*t,2);
    datai_z(i,1:y*t)=data0(1:y*t,3);
    for k=1:y*t
    datai_x(i,(pointNum(1)-k+1))=data(i,3*(y*t-k)+1);  %x座標
    datai_y(i,(pointNum(1)-k+1))=data(i,3*(y*t-k)+2);
    datai_z(i,(pointNum(1)-k+1))=data(i,3*(y*t-k)+3);
    end
%     
end

if i==1
    vn_x(i,1:pointNum(1))=0;   %1列目:node4 2列目:node5 
    vn_y(i,1:pointNum(1))=0;
    vn_z(i,1:pointNum(1))=0;

else
     
%平等に変位させる
    vn_x(i,(y*t+1):(pointNum(1)-y*t))=vn_x(i-1,(y*t+1):(pointNum(1)-y*t))+Fn_x(i-1,(y*t+1):(pointNum(1)-y*t))*dt/2;
    vn_y(i,(y*t+1):(pointNum(1)-y*t))=vn_y(i-1,(y*t+1):(pointNum(1)-y*t))+Fn_y(i-1,(y*t+1):(pointNum(1)-y*t))*dt/2;
    vn_z(i,(y*t+1):(pointNum(1)-y*t))=vn_z(i-1,(y*t+1):(pointNum(1)-y*t))+Fn_z(i-1,(y*t+1):(pointNum(1)-y*t))*dt/2;
    datai_x(i,(y*t+1):(pointNum(1)-y*t))=datai_x(i-1,(y*t+1):(pointNum(1)-y*t))+dt*vn_x(i,(y*t+1):(pointNum(1)-y*t));
    datai_y(i,(y*t+1):(pointNum(1)-y*t))=datai_y(i-1,(y*t+1):(pointNum(1)-y*t))+dt*vn_y(i,(y*t+1):(pointNum(1)-y*t));
    datai_z(i,(y*t+1):(pointNum(1)-y*t))=datai_z(i-1,(y*t+1):(pointNum(1)-y*t))+dt*vn_z(i,(y*t+1):(pointNum(1)-y*t));

% 骨
    for j=1:h-1
        for k=1:y
            if datai_y(i,y*(t-1)+k+j*y*t)>data0(y*(t-1)+k+j*y*t,2)
        datai_y(i,y*(t-1)+k+j*y*t)=data0(y*(t-1)+k+j*y*t,2);
            end
        end
    end

    %4層目まで骨
% for j=1:h-3
%         for k=1:y
%             if datai_y(i,y*(t-1)+k+j*y*t)>data0(y*(t-1)+k+j*y*t,2)
%         datai_y(i,y*(t-1)+k+j*y*t)=data0(y*(t-1)+k+j*y*t,2);
%             end
%         end
% end

    
end

%ばね要素長さ
for s=1:seNum(1);  
  length_s(i,s)=sqrt((datai_x(i,se(s,3))-datai_x(i,se(s,2)))^2+(datai_y(i,se(s,3))-datai_y(i,se(s,2)))^2+(datai_z(i,se(s,3))-datai_z(i,se(s,2)))^2);      
  lengthen_s(i,s)=length_s(i,s)-length_s(1,s);
  
%角度計算 
  if datai_y(i,se(s,3))-datai_y(i,se(s,2))==0&&datai_x(i,se(s,3))-datai_x(i,se(s,2))>=0
      se(s,5)=0;
  elseif datai_y(i,se(s,3))-datai_y(i,se(s,2))==0&&datai_x(i,se(s,3))-datai_x(i,se(s,2))<=0
      se(s,5)=pi;
  elseif datai_x(i,se(s,3))-datai_x(i,se(s,2))==0&&datai_y(i,se(s,3))-datai_y(i,se(s,2))>=0;
      se(s,5)=pi/2;
  elseif datai_x(i,se(s,3))-datai_x(i,se(s,2))==0&&datai_y(i,se(s,3))-datai_y(i,se(s,2))<=0
      se(s,5)=-pi/2;
  elseif datai_x(i,se(s,3))-datai_x(i,se(s,2))<=0&&datai_y(i,se(s,3))-datai_y(i,se(s,2))>=0;
      se(s,5)=atan((datai_y(i,se(s,3))-datai_y(i,se(s,2)))/(datai_x(i,se(s,3))-datai_x(i,se(s,2))))+pi;   %角度β(xy平面)
  elseif datai_x(i,se(s,3))-datai_x(i,se(s,2))<=0&&datai_y(i,se(s,3))-datai_y(i,se(s,2))<=0;
      se(s,5)=atan((datai_y(i,se(s,3))-datai_y(i,se(s,2)))/(datai_x(i,se(s,3))-datai_x(i,se(s,2))))+pi;
  else
      se(s,5)=atan((datai_y(i,se(s,3))-datai_y(i,se(s,2)))/(datai_x(i,se(s,3))-datai_x(i,se(s,2))));
  end
  
  se(s,4)=atan((datai_z(i,se(s,3))-datai_z(i,se(s,2)))/sqrt((datai_x(i,se(s,3))-datai_x(i,se(s,2)))^2+(datai_y(i,se(s,3))-datai_y(i,se(s,2)))^2));   %角度α
end


%% 慣性力
if i==1;
    Mr_x(i,1:pointNum(1))=0;
    Mr_y(i,1:pointNum(1))=0;
    Mr_z(i,1:pointNum(1))=0;
else
Mr_x(i,1:pointNum(1))=mass*(datai_x(i,1:pointNum(1))-datai_x(i-1,1:pointNum(1)))/dt/dt;
Mr_y(i,1:pointNum(1))=mass*(datai_y(i,1:pointNum(1))-datai_y(i-1,1:pointNum(1)))/dt/dt;
Mr_z(i,1:pointNum(1))=mass*(datai_z(i,1:pointNum(1))-datai_z(i-1,1:pointNum(1)))/dt/dt;
end

%% 収縮力(筋力)

Fm_x(i,1:pointNum(1))=0;    Fm_y(i,1:pointNum(1))=0;    Fm_z(i,1:pointNum(1))=0;


%% ばね力(横)
for n=1:6;
for j=1:seNum(1);
    if n==6;
        ul(1,n)=seNum(1);
    end
    if se(j,1)>n-1;
        ul(1,n)=j-1;  %upper limit
        break
    end
end
end

for s=1:seNum(1);
    springF(i,s)=springk*lengthen_s(i,s);
    
    
if se(s,1)==2  %筋線維だけHillのモデル計算
    
    Lcesh(1,s-ul(1,2)) = 0.5 * length_s(1,s);
    f_Lce(i,s-ul(1,2)) = exp(-( lengthen_s(i,s)/Lcesh(1,s-ul(1,2))).^2);
    if i == 1;
        Vce(i,s-ul(1,2)) = 0;
    else;
        Vce(i,s-ul(1,2)) = -(length_s(i,s) - length_s(i-1,s))/0.1;
    end;
    Vvm = 6 * length_s(1,s);
    Vmax(1,s-ul(1,2)) = Vvm .* ( 1 - Ver .* (1 - alpha * f_Lce(i,s-ul(1,2))));
 
    if Vce(i,s-ul(1,2)) <= -Vmax(1,s-ul(1,2));
        f_Vce(i,s-ul(1,2)) = 0;
    elseif Vce(i,s-ul(1,2)) > -Vmax(1,s-ul(1,2))  &&  Vce(i,s-ul(1,2)) < 0;
       f_Vce(i,s-ul(1,2)) = (Vsh * Vmax(1,s-ul(1,2)) + Vsh * Vce(i,s-ul(1,2))) / (Vsh *Vmax(1,s-ul(1,2)) - Vce(i,s-ul(1,2)));
    elseif Vce(i,s-ul(1,2)) >= 0;
        f_Vce(i,s-ul(1,2)) = (Vsh * Vshl * Vmax(1,s-ul(1,2)) + Vml * Vce(i,s-ul(1,2))) / (Vsh * Vshl * Vmax(1,s-ul(1,2)) + Vce(i,s-ul(1,2)));
    end;
    
    length_per(i,s-ul(1,2)) =lengthen_s(i,s)/ length_s(1,s);
    if lengthen_s(i,s)<=0;
    HillPassive(i,s-ul(1,2)) =0;
    else
    HillPassive(i,s-ul(1,2)) = Fmax/exp(PEsh)* (exp(lengthen_s(i,s)/ length_s(1,s)* PEsh/PExm) -1); 
    end
    HillActive(i,s-ul(1,2)) = alpha .* f_Lce(i,s-ul(1,2)) .* f_Vce(i,s-ul(1,2)) .* Fmax'; 
    if i==1
    springF(i,s)=0;
    else
    springF(i,s)=HillActive(i,s-ul(1,2))+HillPassive(i,s-ul(1,2));
    end
end

end

%%質点ごとのspringforce
%分布範囲
for n=1:6;
Fs_x(i,1:pointNum(1),n)=0;
Fs_y(i,1:pointNum(1),n)=0;
Fs_z(i,1:pointNum(1),n)=0;
senum=(1:seNum);

if n==1; 
    Fsl_x(i,1:pointNum(1))=0;
    Fsl_y(i,1:pointNum(1))=0;
    Fsl_z(i,1:pointNum(1))=0;
    Fsr_x(i,1:pointNum(1))=0;   
    Fsr_y(i,1:pointNum(1))=0; 
    Fsr_z(i,1:pointNum(1))=0; 
    
    %左だけ検索
    for k=1:ul(1,n)
    for j=1:ul(1,n) 
        match_l(j,1)=se(j,2)==se(k,2);
    end
    match_lNum=senum(match_l);   %同じ質点を有する三角Noを抽出
    sN=size(match_lNum);
    if sN(2)>1;   
    matchse_l=[];
    for s=1:sN(2);   %matchnumに入っている三角Noのxy方向力をmatchf_Aに表示
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
    
    %右だけ検索
    for k=1:ul(1,n)
    for j=1:ul(1,n) 
        match_r(j,1)=se(j,3)==se(k,3);
    end
    match_rNum=senum(match_r);   %同じ質点を有する三角Noを抽出
    sN=size(match_rNum);
    if sN(2)>1;   
    matchse_r=[];
    for s=1:sN(2);   %matchnumに入っている三角Noのxy方向力をmatchf_Aに表示
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
    
    
    %左だけ検索
    for k=ul(1,n-1)+1:ul(1,n)
    for j=ul(1,n-1)+1:ul(1,n) 
        match_l(j,1)=se(j,2)==se(k,2);
    end
    match_lNum=senum(match_l);   %同じ質点を有する三角Noを抽出
    sN=size(match_lNum);
    if sN(2)>1;   
    matchse_l=[];
    for s=1:sN(2);   %matchnumに入っている三角Noのxy方向力をmatchf_Aに表示
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
    
    %右だけ検索
    for k=ul(1,n-1)+1:ul(1,n)
    for j=ul(1,n-1)+1:ul(1,n) 
        match_r(j,1)=se(j,3)==se(k,3);
    end
    match_rNum=senum(match_r);   %同じ質点を有する三角Noを抽出
    sN=size(match_rNum);
    if sN(2)>1;   
    matchse_r=[];
    for s=1:sN(2);   %matchnumに入っている三角Noのxy方向力をmatchf_Aに表示
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


%% 拘束条件

    con_x(i,1:pointNum(1))=-1*conk*(datai_x(i,1:pointNum(1))-datai_x(1,1:pointNum(1)));
    con_y(i,1:pointNum(1))=-1*conk*(datai_y(i,1:pointNum(1))-datai_y(1,1:pointNum(1)));
    con_z(i,1:pointNum(1))=-1*conk*(datai_z(i,1:pointNum(1))-datai_z(1,1:pointNum(1)));


    if i==1;
        vis_x(i,1:pointNum(1))=0;
        vis_y(i,1:pointNum(1))=0;
        vis_z(i,1:pointNum(1))=0;
    else
        vis_x(i,1:pointNum(1))=-1*conc*(datai_x(i,1:pointNum(1))-datai_x(i-1,1:pointNum(1)))/dt;
        vis_y(i,1:pointNum(1))=-1*conc*(datai_y(i,1:pointNum(1))-datai_y(i-1,1:pointNum(1)))/dt;
        vis_z(i,1:pointNum(1))=-1*conc*(datai_z(i,1:pointNum(1))-datai_z(i-1,1:pointNum(1)))/dt;
    end



%% 体積保存力

for j=1:tetraNum(1)/5;
for k=1:tetraNum(2);
    if i==1;
    V0(k+(j-1)*tetraNum(2),1)=1/6*abs(dot(cross((data0(tetra(3+5*(j-1),k),:).'-data0(tetra(2+5*(j-1),k),:).'),(data0(tetra(4+5*(j-1),k),:).'-data0(tetra(2+5*(j-1),k),:).')),((data0(tetra(5+5*(j-1),k),:).'-data0(tetra(2+5*(j-1),k),:).'))));
    end
    %fvAを求めるのに必要なﾍﾞｸﾄﾙ
    Vec1(:,1)=[datai_x(i,tetra(3+5*(j-1),k));datai_y(i,tetra(3+5*(j-1),k));datai_z(i,tetra(3+5*(j-1),k))]-[datai_x(i,tetra(5+5*(j-1),k));datai_y(i,tetra(5+5*(j-1),k));datai_z(i,tetra(5+5*(j-1),k))];   %B-D
    Vec2(:,1)=[datai_x(i,tetra(4+5*(j-1),k));datai_y(i,tetra(4+5*(j-1),k));datai_z(i,tetra(4+5*(j-1),k))]-[datai_x(i,tetra(5+5*(j-1),k));datai_y(i,tetra(5+5*(j-1),k));datai_z(i,tetra(5+5*(j-1),k))];   %C-D
    
    %fvB
    Vec1(:,2)=[datai_x(i,tetra(4+5*(j-1),k));datai_y(i,tetra(4+5*(j-1),k));datai_z(i,tetra(4+5*(j-1),k))]-[datai_x(i,tetra(2+5*(j-1),k));datai_y(i,tetra(2+5*(j-1),k));datai_z(i,tetra(2+5*(j-1),k))];   %C-A
    Vec2(:,2)=[datai_x(i,tetra(5+5*(j-1),k));datai_y(i,tetra(5+5*(j-1),k));datai_z(i,tetra(5+5*(j-1),k))]-[datai_x(i,tetra(2+5*(j-1),k));datai_y(i,tetra(2+5*(j-1),k));datai_z(i,tetra(2+5*(j-1),k))];   %D-A
    
    %fvC
    Vec1(:,3)=[datai_x(i,tetra(5+5*(j-1),k));datai_y(i,tetra(5+5*(j-1),k));datai_z(i,tetra(5+5*(j-1),k))]-[datai_x(i,tetra(3+5*(j-1),k));datai_y(i,tetra(3+5*(j-1),k));datai_z(i,tetra(3+5*(j-1),k))];   %D-B
    Vec2(:,3)=[datai_x(i,tetra(2+5*(j-1),k));datai_y(i,tetra(2+5*(j-1),k));datai_z(i,tetra(2+5*(j-1),k))]-[datai_x(i,tetra(3+5*(j-1),k));datai_y(i,tetra(3+5*(j-1),k));datai_z(i,tetra(3+5*(j-1),k))];   %A-B
    
    %fvD
    Vec1(:,4)=[datai_x(i,tetra(2+5*(j-1),k));datai_y(i,tetra(2+5*(j-1),k));datai_z(i,tetra(2+5*(j-1),k))]-[datai_x(i,tetra(4+5*(j-1),k));datai_y(i,tetra(4+5*(j-1),k));datai_z(i,tetra(4+5*(j-1),k))];   %A-C
    Vec2(:,4)=[datai_x(i,tetra(3+5*(j-1),k));datai_y(i,tetra(3+5*(j-1),k));datai_z(i,tetra(3+5*(j-1),k))]-[datai_x(i,tetra(4+5*(j-1),k));datai_y(i,tetra(4+5*(j-1),k));datai_z(i,tetra(4+5*(j-1),k))];   %B-C
    
    V(k+(j-1)*tetraNum(2),1)=1/6*abs(dot(cross((-1)*Vec2(:,3),Vec1(:,2)),Vec2(:,2)));

%体積時系列
if data(i,3)-data(1,3)>0;    %こっちが標準
         fvA(:,k+(j-1)*tetraNum(2))=1/6*kv*(V(k+(j-1)*tetraNum(2),1)-V0(k+(j-1)*tetraNum(2),1))/V0(k+(j-1)*tetraNum(2),1)^2*cross(Vec2(:,1),Vec1(:,1));
         fvB(:,k+(j-1)*tetraNum(2))=1/6*kv*(V(k+(j-1)*tetraNum(2),1)-V0(k+(j-1)*tetraNum(2),1))/V0(k+(j-1)*tetraNum(2),1)^2*cross(Vec1(:,2),Vec2(:,2));
         fvC(:,k+(j-1)*tetraNum(2))=1/6*kv*(V(k+(j-1)*tetraNum(2),1)-V0(k+(j-1)*tetraNum(2),1))/V0(k+(j-1)*tetraNum(2),1)^2*cross(Vec2(:,3),Vec1(:,3));
         fvD(:,k+(j-1)*tetraNum(2))=1/6*kv*(V(k+(j-1)*tetraNum(2),1)-V0(k+(j-1)*tetraNum(2),1))/V0(k+(j-1)*tetraNum(2),1)^2*cross(Vec1(:,4),Vec2(:,4));
else
         fvA(:,k+(j-1)*tetraNum(2))=1/6*kv*(V(k+(j-1)*tetraNum(2),1)-V0(k+(j-1)*tetraNum(2),1))/V0(k+(j-1)*tetraNum(2),1)^2*cross(Vec1(:,1),Vec2(:,1));
         fvB(:,k+(j-1)*tetraNum(2))=1/6*kv*(V(k+(j-1)*tetraNum(2),1)-V0(k+(j-1)*tetraNum(2),1))/V0(k+(j-1)*tetraNum(2),1)^2*cross(Vec2(:,2),Vec1(:,2));
         fvC(:,k+(j-1)*tetraNum(2))=1/6*kv*(V(k+(j-1)*tetraNum(2),1)-V0(k+(j-1)*tetraNum(2),1))/V0(k+(j-1)*tetraNum(2),1)^2*cross(Vec1(:,3),Vec2(:,3));
         fvD(:,k+(j-1)*tetraNum(2))=1/6*kv*(V(k+(j-1)*tetraNum(2),1)-V0(k+(j-1)*tetraNum(2),1))/V0(k+(j-1)*tetraNum(2),1)^2*cross(Vec2(:,4),Vec1(:,4));
end
end
end

FvA(1:pointNum(1),1:3)=0;
tetranum=(1:tetraNum(1)*tetraNum(2)/5);
%a行b列目のnode番号を有するtetraNumに1 
for m=1:tetraNum(1)/5;
for k=1:tetraNum(2);
for j=1:tetraNum(1)/5;
matchA(1+(j-1)*tetraNum(2):tetraNum(2)*j,1)=tetra(2+5*(j-1),:).'==tetra(2+5*(m-1),k);    %同じ質点のある三角Noが　ある場合：1 ない：0
end
matchAnum=tetranum(matchA);   %同じ質点を有する三角Noを抽出
mN=size(matchAnum);      %1つ以上なら同じ質点を有する三角Noがほかにあるということ
if mN(2)>1;   
    matchf_A=[];
for n=1:mN(2);   %matchnumに入っている三角Noのxy方向力をmatchf_Aに表示
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
%a行b列目のnode番号を有するtetraNumに1 
for m=1:tetraNum(1)/5;
for k=1:tetraNum(2);
for j=1:tetraNum(1)/5;
matchB(1+(j-1)*tetraNum(2):tetraNum(2)*j,1)=tetra(3+5*(j-1),:).'==tetra(3+5*(m-1),k);    %同じ質点のある三角Noが　ある場合：1 ない：0
end
matchBnum=tetranum(matchB);   %同じ質点を有する三角Noを抽出
mN=size(matchBnum);      %1つ以上なら同じ質点を有する三角Noがほかにあるということ
if mN(2)>1;   
    matchf_B=[];
for n=1:mN(2);   %matchnumに入っている三角Noのxy方向力をmatchf_Bに表示
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
%a行b列目のnode番号を有するtetraNumに1 
for m=1:tetraNum(1)/5;
for k=1:tetraNum(2);
for j=1:tetraNum(1)/5;
matchC(1+(j-1)*tetraNum(2):tetraNum(2)*j,1)=tetra(4+5*(j-1),:).'==tetra(4+5*(m-1),k);    %同じ質点のある三角Noが　ある場合：1 ない：0
end
matchCnum=tetranum(matchC);   %同じ質点を有する三角Noを抽出
mN=size(matchCnum);      %1つ以上なら同じ質点を有する三角Noがほかにあるということ
if mN(2)>1;   
    matchf_C=[];
for n=1:mN(2);   %matchnumに入っている三角Noのxy方向力をmatchf_Cに表示
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
%a行b列目のnode番号を有するtetraNumに1 
for m=1:tetraNum(1)/5;
for k=1:tetraNum(2);
for j=1:tetraNum(1)/5;
matchD(1+(j-1)*tetraNum(2):tetraNum(2)*j,1)=tetra(5+5*(j-1),:).'==tetra(5+5*(m-1),k);    %同じ質点のある三角Noが　ある場合：1 ない：0
end
matchDnum=tetranum(matchD);   %同じ質点を有する三角Noを抽出
mN=size(matchDnum);      %1つ以上なら同じ質点を有する三角Noがほかにあるということ
if mN(2)>1;   
    matchf_D=[];
for n=1:mN(2);   %matchnumに入っている三角Noのxy方向力をmatchf_Dに表示
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


if i==1;
    Fv_x(i,1:pointNum(1))=0;
    Fv_y(i,1:pointNum(1))=0;
    Fv_z(i,1:pointNum(1))=0;
else

Fv_x(i,1:pointNum(1))=FvA(1:pointNum(1),1).'+FvB(1:pointNum(1),1).'+FvC(1:pointNum(1),1).'+FvD(1:pointNum(1),1).';
Fv_y(i,1:pointNum(1))=FvA(1:pointNum(1),2).'+FvB(1:pointNum(1),2).'+FvC(1:pointNum(1),2).'+FvD(1:pointNum(1),2).';
Fv_z(i,1:pointNum(1))=FvA(1:pointNum(1),3).'+FvB(1:pointNum(1),3).'+FvC(1:pointNum(1),3).'+FvD(1:pointNum(1),3).';

end

%0：横　1：縦　2：高さ　3：xy平面斜め　4：xz平面斜め　5：yz平面斜め　6：空間斜め

Fn_x(i,:)=FsSum_x(i,:)+Fv_x(i,:)+Mr_x(i,:)+con_x(i,:)+vis_x(i,:);
Fn_y(i,:)=FsSum_y(i,:)+Fv_y(i,:)+Mr_y(i,:)+con_y(i,:)+vis_y(i,:);
Fn_z(i,:)=FsSum_z(i,:)+Fv_z(i,:)+Mr_z(i,:)+con_z(i,:)+vis_z(i,:);



vn_x(i,1:pointNum(1))=vn_x(i,1:pointNum(1))+Fn_x(i,1:pointNum(1))/2*dt;
vn_y(i,1:pointNum(1))=vn_y(i,1:pointNum(1))+Fn_y(i,1:pointNum(1))/2*dt;
vn_z(i,1:pointNum(1))=vn_z(i,1:pointNum(1))+Fn_z(i,1:pointNum(1))/2*dt;

%% 描画
% if i==1||mod(i,10)==0;
%     if i==1;
%         l=1;
%     else
%     l=i/10;
%     end
% for j=1:(t-1);  %middle point
% mp(j,:)=y+1+(y*2-1)*(j-1):(y*2-1+(y*2-1)*(j-1));
% end
% 
% 
% %%横線
% %伸ばす線がないnode
% for h=1:h;
%     for t=1:t;
%         limit(h,t)=y*t+y*t*(h-1);
%     end
% end
% 
% for j=1:ul(1,1);
% x_yoko(j,1)=datai_x(i,se(j,2));
% x_yoko(j,2)=datai_x(i,se(j,3));
% y_yoko(j,1)=datai_y(i,se(j,2));
% y_yoko(j,2)=datai_y(i,se(j,3));
% z_yoko(j,1)=datai_z(i,se(j,2));
% z_yoko(j,2)=datai_z(i,se(j,3));
% plot3(x_yoko(j,:)*1000,y_yoko(j,:)*1000,z_yoko(j,:)*1000,'k');
% hold on
% end
% 
% 
% 
% %%縦線
% %伸ばす線がないnode
% for j=1+ul(1,1):ul(1,2);
% x_tate(j,1)=datai_x(i,se(j,2));
% x_tate(j,2)=datai_x(i,se(j,3));
% y_tate(j,1)=datai_y(i,se(j,2));
% y_tate(j,2)=datai_y(i,se(j,3));
% z_tate(j,1)=datai_z(i,se(j,2));
% z_tate(j,2)=datai_z(i,se(j,3));
% plot3(x_tate(j,:)*1000,y_tate(j,:)*1000,z_tate(j,:)*1000,'k');
% hold on
% end
% 
% 
% %%高さ線
% for j=1+ul(1,2):ul(1,3);
% x_high(j,1)=datai_x(i,se(j,2));
% x_high(j,2)=datai_x(i,se(j,3));
% y_high(j,1)=datai_y(i,se(j,2));
% y_high(j,2)=datai_y(i,se(j,3));
% z_high(j,1)=datai_z(i,se(j,2));
% z_high(j,2)=datai_z(i,se(j,3));
% plot3(x_high(j,:)*1000,y_high(j,:)*1000,z_high(j,:)*1000,'k');
% hold on
% end
% 
% %%xy平面斜め線(3)
% for j=1+ul(1,3):ul(1,4);
% x_xy(j,1)=datai_x(i,se(j,2));
% x_xy(j,2)=datai_x(i,se(j,3));
% y_xy(j,1)=datai_y(i,se(j,2));
% y_xy(j,2)=datai_y(i,se(j,3));
% z_xy(j,1)=datai_z(i,se(j,2));
% z_xy(j,2)=datai_z(i,se(j,3));
% plot3(x_xy(j,:)*1000,y_xy(j,:)*1000,z_xy(j,:)*1000,'k');
% hold on
% end
% 
% 
% %%xz平面斜め線(4)
% for j=1+ul(1,4):ul(1,5);
% x_xz(j,1)=datai_x(i,se(j,2));
% x_xz(j,2)=datai_x(i,se(j,3));
% y_xz(j,1)=datai_y(i,se(j,2));
% y_xz(j,2)=datai_y(i,se(j,3));
% z_xz(j,1)=datai_z(i,se(j,2));
% z_xz(j,2)=datai_z(i,se(j,3));
% plot3(x_xz(j,:)*1000,y_xz(j,:)*1000,z_xz(j,:)*1000,'k');
% hold on
% end
% 
% 
% %%yz平面斜め線(5)
% for j=1+ul(1,5):ul(1,6);
% x_yz(j,1)=datai_x(i,se(j,2));
% x_yz(j,2)=datai_x(i,se(j,3));
% y_yz(j,1)=datai_y(i,se(j,2));
% y_yz(j,2)=datai_y(i,se(j,3));
% z_yz(j,1)=datai_z(i,se(j,2));
% z_yz(j,2)=datai_z(i,se(j,3));
% plot3(x_yz(j,:)*1000,y_yz(j,:)*1000,z_yz(j,:)*1000,'k');
% hold on
% end
% 
% 
% 
% for k=1:pointNum(1);
% %if k==1;
% 
% plot3(datai_x(i,k)*1000,datai_y(i,k)*1000,datai_z(i,k)*1000,'.','MarkerSize',30,'MarkerEdgeColor','k'); 
% hold on
% % else
% %     hold on
% %     plot(datai_x(i,k)*1000,datai_y(i,k)*1000,'.','MarkerSize',10,'MarkerEdgeColor','k');   
% 
% % 
% % hold off
% % end
% end
% 
% % xlim([67.5,367.5]);
% % ylim([-100,200]);
% % zlim([1000,1300]);
% Frame(l) = getframe(1);
% 
% 
% hold off
% end
end
toc
F_fib=springF(:,ul(1,2)+1:ul(1,2)+y*t*(h-1));

csvwrite('output\data_x.csv',datai_x)
csvwrite('output\data_y.csv',datai_y)
csvwrite('output\data_z.csv',datai_z)


% v = VideoWriter('model_3d_7.avi');
% 
% open(v);
% writeVideo(v,Frame);
% close(v);
% 


