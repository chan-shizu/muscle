function parameter
%ばね+筋線維のモデル
% clc; %コマンドウェイドウからすべてのテキストをクリア
% clear; %ワークスペースからアイテムを削除し，システムメモリ開放
warning('on','all');%すべての警告を有効にする
warning;
muscleNames = loadMuscleName();
muscle_name = [muscleNames{1}];

E=5*10^4; %file:///C:/Users/%E5%8D%A0%E9%83%A8%E3%80%80%E9%BA%BB%E9%87%8C%E5%AD%90/Downloads/nagano_05-02-05%20(1).pdf
v=0.49;
%k=8500;
%k=8237;
%ks=166.7558;
%ks=170;
Fmax=5;
PEsh = 4;
PExm = 0.4;

divisionMuscle = readmatrix("divisionMuscle.csv");
y=divisionMuscle(1);
t=divisionMuscle(2);
h=divisionMuscle(3);

file_name_data0 = ['muscle\',muscle_name, '_min_final.csv'];
% file_name_data0 = strcat("C:\Users\bubbl\Documents\shizuya_M1\DefMuscle_for_TUS_2018\matlab\surface2grid_0608\muscle\min\", muscle_name,"_min_final.csv")
% file_name_data0 = strcat("C:\Users\bubbl\Documents\shizuya_M1\DefMuscle_for_TUS_2018\matlab\surface2grid_0608\muscle\min\data0_arm_rot.csv")
data0=csvread(file_name_data0, 0, 0)/1000; %mmからmに変換

pointNum = size(data0); 
dt=1*10^(-3); 
%% se作成
for i=1:t
yoko(1+(y-1)*(i-1):(y-1)*i,2)=1+y*(i-1):y*i-1;
yoko(1+(y-1)*(i-1):(y-1)*i,3)=2+y*(i-1):y*i;
end

%横
for j=0:h-1
for i=1:t
yoko(1+(y-1)*(i-1)+j*(y-1)*t:(y-1)*i+j*(y-1)*t,2)=1+y*(i-1)+t*y*j:y*i-1+t*y*j;
yoko(1+(y-1)*(i-1)+j*(y-1)*t:(y-1)*i+j*(y-1)*t,3)=2+y*(i-1)+t*y*j:y*i+t*y*j;
end
end

%縦
for i=0:h-1
tate(1+y*(t-1)*i:y*(t-1)*(i+1),2)=1+t*y*i:y*(t-1)+t*y*i;
tate(1+y*(t-1)*i:y*(t-1)*(i+1),3)=1+y+t*y*i:t*y*(i+1);
end
tate(:,1)=1;

%高さ
takasa(1:(h-1)*t*y,2)=1:(h-1)*t*y;
takasa(1:(h-1)*t*y,3)=1+t*y:h*t*y;
takasa(:,1)=2;

%xy
for k=0:h-1
for j=0:t-2
for i=0:y-2
xy(1+2*i+2*(y-1)*j+((y-1)*(t-1)*2*k),2)=1+i+y*j+t*y*k;    xy(1+2*i+2*(y-1)*j+((y-1)*(t-1)*2*k),3)=2+y+i+y*j+t*y*k;
xy(2+2*i+2*(y-1)*j+((y-1)*(t-1)*2*k),2)=2+i+y*j+t*y*k;    xy(2+2*i+2*(y-1)*j+((y-1)*(t-1)*2*k),3)=1+y+i+y*j+t*y*k;
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
file_name_se = ['parameter\',muscle_name, '_se.csv'];
% file_name_se = strcat("C:\Users\bubbl\Documents\shizuya_M1\DefMuscle_for_TUS_2018\matlab\surface2grid_0608\parameter\", muscle_name, "_se.csv")
csvwrite(file_name_se, se);
%% tetra
for j=1:(t-1)*(h-1)
for i=1:y-1
tetra(1+5*(j-1),1+10*(i-1):10*i)=i;
end
end

%横に作っていく
for i=1:y-1
tetra(2,1+10*(i-1))=1+(i-1);     tetra(3,1+10*(i-1))=1+y+y*t+(i-1);   tetra(4,1+10*(i-1))=1+y*t+(i-1);     tetra(5,1+10*(i-1))=2+y*t+(i-1);
tetra(2,2+10*(i-1))=2+y+(i-1);   tetra(3,2+10*(i-1))=2+y+y*t+(i-1);   tetra(4,2+10*(i-1))=1+y+y*t+(i-1);   tetra(5,2+10*(i-1))=2+y*t+(i-1);
tetra(2,3+10*(i-1))=1+(i-1);     tetra(3,3+10*(i-1))=2+(i-1);         tetra(4,3+10*(i-1))=2+y+(i-1);       tetra(5,3+10*(i-1))=2+y*t+(i-1);
tetra(2,4+10*(i-1))=1+(i-1);     tetra(3,4+10*(i-1))=1+y+y*t+(i-1);   tetra(4,4+10*(i-1))=2+y+(i-1);       tetra(5,4+10*(i-1))=1+y+(i-1);
tetra(2,5+10*(i-1))=1+(i-1);     tetra(3,5+10*(i-1))=2+y+(i-1);       tetra(4,5+10*(i-1))=1+y+y*t+(i-1);   tetra(5,5+10*(i-1))=2+y*t+(i-1);
tetra(2,6+10*(i-1))=2+(i-1);     tetra(3,6+10*(i-1))=1+y*t+(i-1);     tetra(4,6+10*(i-1))=2+y*t+(i-1);     tetra(5,6+10*(i-1))=2+y+y*t+(i-1);
tetra(2,7+10*(i-1))=1+y+(i-1);   tetra(3,7+10*(i-1))=1+y+y*t+(i-1);   tetra(4,7+10*(i-1))=1+y*t+(i-1);     tetra(5,7+10*(i-1))=2+y+y*t+(i-1);
tetra(2,8+10*(i-1))=1+(i-1);     tetra(3,8+10*(i-1))=2+(i-1);         tetra(4,8+10*(i-1))=1+y+(i-1);       tetra(5,8+10*(i-1))=1+y*t+(i-1);
tetra(2,9+10*(i-1))=2+(i-1);     tetra(3,9+10*(i-1))=2+y+(i-1);       tetra(4,9+10*(i-1))=1+y+(i-1);       tetra(5,9+10*(i-1))=2+y+y*t+(i-1);
tetra(2,10+10*(i-1))=2+(i-1);    tetra(3,10+10*(i-1))=1+y+(i-1);      tetra(4,10+10*(i-1))=1+y*t+(i-1);    tetra(5,10+10*(i-1))=2+y+y*t+(i-1);
end


%縦
for i=1:t-2
tetra(2+5*i:5+5*i,1:10*(y-1))=tetra(2+5*(i-1):5+5*(i-1),1:10*(y-1))+y;     
end

%高さ
for i=1:h-2
    for k=1:t-1
        tetra(2+5*(k-1)+5*(t-1)*i:5+5*(k-1)+5*(t-1)*i,1:10*(y-1))=tetra(2+5*(k-1)+5*(t-1)*(i-1):5+5*(k-1)+5*(t-1)*(i-1),1:10*(y-1))+y*t;
    end
end
tetraNum=size(tetra);
file_name_tetra = ['parameter\',muscle_name, '_tetra.csv'];
% file_name_tetra = strcat("C:\Users\bubbl\Documents\shizuya_M1\DefMuscle_for_TUS_2018\matlab\surface2grid_0608\parameter\", muscle_name, "_tetra.csv")
csvwrite(file_name_tetra, tetra);

%% 3D
tetraNum=size(tetra);
for j=1:tetraNum(1)/5
    for i=1:tetraNum(2)
        x(i+(j-1)*tetraNum(2),1)=tetra(1+5*(j-1),i);  x(i+(j-1)*tetraNum(2),2)=data0(tetra(2+5*(j-1),i),1);  x(i+(j-1)*tetraNum(2),3)=data0(tetra(3+5*(j-1),i),1);  x(i+(j-1)*tetraNum(2),4)=data0(tetra(4+5*(j-1),i),1);  x(i+(j-1)*tetraNum(2),5)=data0(tetra(5+5*(j-1),i),1);
        y(i+(j-1)*tetraNum(2),1)=tetra(1+5*(j-1),i);  y(i+(j-1)*tetraNum(2),2)=data0(tetra(2+5*(j-1),i),2);  y(i+(j-1)*tetraNum(2),3)=data0(tetra(3+5*(j-1),i),2);  y(i+(j-1)*tetraNum(2),4)=data0(tetra(4+5*(j-1),i),2);  y(i+(j-1)*tetraNum(2),5)=data0(tetra(5+5*(j-1),i),2);
        z(i+(j-1)*tetraNum(2),1)=tetra(1+5*(j-1),i);  z(i+(j-1)*tetraNum(2),2)=data0(tetra(2+5*(j-1),i),3);  z(i+(j-1)*tetraNum(2),3)=data0(tetra(3+5*(j-1),i),3);  z(i+(j-1)*tetraNum(2),4)=data0(tetra(4+5*(j-1),i),3);  z(i+(j-1)*tetraNum(2),5)=data0(tetra(5+5*(j-1),i),3);
        V0(i+(j-1)*tetraNum(2),1)=1/6*dot(cross((data0(tetra(3+5*(j-1),i),:).'-data0(tetra(2+5*(j-1),i),:).'),(data0(tetra(4+5*(j-1),i),:).'-data0(tetra(2+5*(j-1),i),:).')),((data0(tetra(5+5*(j-1),i),:).'-data0(tetra(2+5*(j-1),i),:).')));

    end
end

%KCSTを求める
D=E/(1+v)/(1-2*v)*[1-v v v 0 0 0; 
                   v 1-v v 0 0 0
                   v v 1-v 0 0 0
                   0 0 0 1/2-v 0 0 
                   0 0 0 0 1/2-v 0
                   0 0 0 0 0 1/2-v];

for i=1:tetraNum(2)*tetraNum(1)/5
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


%KMSMを求める
for i=1:tetraNum(2)*tetraNum(1)/5
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


%KCを求める
for i=1:tetraNum(2)*tetraNum(1)/5
dr1(:,i)=1/(6*V0(i,1))*cross([x(i,2)-x(i,4); y(i,2)-y(i,4); z(i,2)-z(i,4)], [x(i,3)-x(i,4); y(i,3)-y(i,4); z(i,3)-z(i,4)]); %2行1列
dr2(:,i)=1/(6*V0(i,1))*cross([x(i,3)-x(i,1); y(i,3)-y(i,1); z(i,3)-z(i,1)], [x(i,4)-x(i,1); y(i,4)-y(i,1); z(i,4)-z(i,1)]);
dr3(:,i)=1/(6*V0(i,1))*cross([x(i,4)-x(i,2); y(i,4)-y(i,2); z(i,4)-z(i,2)], [x(i,1)-x(i,2); y(i,1)-y(i,2); z(i,1)-z(i,2)]);
dr4(:,i)=1/(6*V0(i,1))*cross([x(i,1)-x(i,3); y(i,1)-y(i,3); z(i,1)-z(i,3)], [x(i,2)-x(i,3); y(i,2)-y(i,3); z(i,2)-z(i,3)]);

K11(:,:,i)=dr1(:,i)*dr1(:,i).';    K12(:,:,i)=dr1(:,i)*dr2(:,i).';    K13(:,:,i)=dr1(:,i)*dr3(:,i).';    K14(:,:,i)=dr1(:,i)*dr4(:,i).';
K21(:,:,i)=dr2(:,i)*dr1(:,i).';    K22(:,:,i)=dr2(:,i)*dr2(:,i).';    K23(:,:,i)=dr2(:,i)*dr3(:,i).';    K24(:,:,i)=dr2(:,i)*dr4(:,i).';
K31(:,:,i)=dr3(:,i)*dr1(:,i).';    K32(:,:,i)=dr3(:,i)*dr2(:,i).';    K33(:,:,i)=dr3(:,i)*dr3(:,i).';    K34(:,:,i)=dr3(:,i)*dr4(:,i).';
K41(:,:,i)=dr4(:,i)*dr1(:,i).';    K42(:,:,i)=dr4(:,i)*dr2(:,i).';    K43(:,:,i)=dr4(:,i)*dr3(:,i).';    K44(:,:,i)=dr4(:,i)*dr4(:,i).';

KC(:,:,i)=[K11(:,:,i) K12(:,:,i) K13(:,:,i) K14(:,:,i);K21(:,:,i) K22(:,:,i) K23(:,:,i) K24(:,:,i);K31(:,:,i) K32(:,:,i) K33(:,:,i) K34(:,:,i);K41(:,:,i) K42(:,:,i) K43(:,:,i) K44(:,:,i)];

end


%% 最小二乗法 参考(http://www.eli.hokkai-s-u.ac.jp/~kikuchi/ma2/chap08.html)
for i=1:tetraNum(2)*tetraNum(1)/5
    for j=1:12
        Y(12*(j-1)+1:12*j,1,i)=KCST(j,1:12,i);
        X(12*(j-1)+1:12*j,1,i)=KMSM(j,1:12,i);   %k
        X(12*(j-1)+1:12*j,2,i)=KC(j,1:12,i);   %ks
    end  
a(:,:,i)=X(:,:,i)\Y(:,:,i);
end

ave(1,1)=mean(a(1,1,:)) %k
ave(2,1)=mean(a(2,1,:)) %ks
file_name_a = ['parameter\',muscle_name, '_a.csv'];
% file_name_a = strcat("C:\Users\bubbl\Documents\shizuya_M1\DefMuscle_for_TUS_2018\matlab\surface2grid_0608\parameter\", muscle_name, "_a.csv")
csvwrite(file_name_a, ave);
end

