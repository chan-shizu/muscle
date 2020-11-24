clc;
clear;
muscleNames = loadMuscleName();
prompt = '上腕二頭筋なら1, オトガイ舌骨筋なら2, 茎突舌骨筋なら3, 顎二腹筋前腹なら4, 顎二腹筋後腹なら5を押してね';
muscleNumber = inputdlg(prompt,...
             'choose muscle', [1 50])
muscleNumber = str2num(muscleNumber{1})
muscle_name = [muscleNames{muscleNumber,1}];
%se= csvread('se_arm.csv', 0, 0);
se_name = ['parameter\',muscle_name, '_se.csv'];
se= csvread(se_name, 0, 0);
seNum=size(se);
% file_name_data = ['muscle\data_',muscle_name, '.csv'];%"data_rot_h.csv"
% data= csvread(file_name_data, 0, 0);
% timeNum = size(data);
divisionMuscle = readmatrix("divisionMuscle.csv");
p_y=divisionMuscle(1);
p_t=divisionMuscle(2);
p_h=divisionMuscle(3);

% fiber_force=csvread('fiber_force.csv', 0, 0);
file_name_data0 = strcat("muscle\", muscle_name, "_min_final.csv");%"data0_rot_h.csv"
data0=csvread(file_name_data0, 0, 0);%/1000;
datai_x = data0(:,1).';
datai_y = data0(:,2).';
datai_z = data0(:,3).';
% data_x_name = ['output\', muscle_name,'_data_x.csv'];
% data_y_name = ['output\', muscle_name,'_data_y.csv'];
% data_z_name = ['output\', muscle_name,'_data_z.csv'];
% datai_x=csvread(data_x_name, 0, 0);
% datai_y=csvread(data_y_name, 0, 0);
% datai_z=csvread(data_z_name, 0, 0);

for n=1:6
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

%%横線
%伸ばす線がないnode
for j=1:p_h
    %それぞれの層の1番目と2番目の質点をつなぐ(5×5×6の場合1と2, 25と26など)
    for k=1:p_y-1
        x_yoko(k+(j-1)*(p_y-1)*p_t,1)=datai_x(1,se(k+(j-1)*(p_y-1)*p_t,2));
        x_yoko(k+(j-1)*(p_y-1)*p_t,2)=datai_x(1,se(k+(j-1)*(p_y-1)*p_t,3));
        y_yoko(k+(j-1)*(p_y-1)*p_t,1)=datai_y(1,se(k+(j-1)*(p_y-1)*p_t,2));
        y_yoko(k+(j-1)*(p_y-1)*p_t,2)=datai_y(1,se(k+(j-1)*(p_y-1)*p_t,3));
        z_yoko(k+(j-1)*(p_y-1)*p_t,1)=datai_z(1,se(k+(j-1)*(p_y-1)*p_t,2));
        z_yoko(k+(j-1)*(p_y-1)*p_t,2)=datai_z(1,se(k+(j-1)*(p_y-1)*p_t,3));
        plot3(x_yoko(k+(j-1)*(p_y-1)*p_t,:)*1000,y_yoko(k+(j-1)*(p_y-1)*p_t,:)*1000,z_yoko(k+(j-1)*(p_y-1)*p_t,:)*1000,'k');
        hold on
    end
    for k=1+(p_y-1)*(p_t-1):(p_y-1)*p_t
        x_yoko(k+(j-1)*(p_y-1)*p_t,1)=datai_x(1,se(k+(j-1)*(p_y-1)*p_t,2));
        x_yoko(k+(j-1)*(p_y-1)*p_t,2)=datai_x(1,se(k+(j-1)*(p_y-1)*p_t,3));
        y_yoko(k+(j-1)*(p_y-1)*p_t,1)=datai_y(1,se(k+(j-1)*(p_y-1)*p_t,2));
        y_yoko(k+(j-1)*(p_y-1)*p_t,2)=datai_y(1,se(k+(j-1)*(p_y-1)*p_t,3));
        z_yoko(k+(j-1)*(p_y-1)*p_t,1)=datai_z(1,se(k+(j-1)*(p_y-1)*p_t,2));
        z_yoko(k+(j-1)*(p_y-1)*p_t,2)=datai_z(1,se(k+(j-1)*(p_y-1)*p_t,3));
        plot3(x_yoko(k+(j-1)*(p_y-1)*p_t,:)*1000,y_yoko(k+(j-1)*(p_y-1)*p_t,:)*1000,z_yoko(k+(j-1)*(p_y-1)*p_t,:)*1000,'k');
        hold on
    end
end



%%縦線
%伸ばす線がないnode
%         for k=1:p_h;
%             for j=1:p_y:1+(p_t-2)*p_y;
%                 x_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),1)=datai_x(1,se(j+(k-1)*p_y*(p_t-1)+ul(1,1),2));
%                 x_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),2)=datai_x(1,se(j+(k-1)*p_y*(p_t-1)+ul(1,1),3));
%                 y_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),1)=datai_y(1,se(j+(k-1)*p_y*(p_t-1)+ul(1,1),2));
%                 y_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),2)=datai_y(1,se(j+(k-1)*p_y*(p_t-1)+ul(1,1),3));
%                 z_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),1)=datai_z(1,se(j+(k-1)*p_y*(p_t-1)+ul(1,1),2));
%                 z_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),2)=datai_z(1,se(j+(k-1)*p_y*(p_t-1)+ul(1,1),3));
%                 plot3(x_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),:)*1000,y_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),:)*1000,z_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),:)*1000,'k');
%                 hold on
%             end
%             for j=p_y:p_y:p_y+(p_t-2)*p_y;
%                 x_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),1)=datai_x(1,se(j+(k-1)*p_y*(p_t-1)+ul(1,1),2));
%                 x_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),2)=datai_x(1,se(j+(k-1)*p_y*(p_t-1)+ul(1,1),3));
%                 y_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),1)=datai_y(1,se(j+(k-1)*p_y*(p_t-1)+ul(1,1),2));
%                 y_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),2)=datai_y(1,se(j+(k-1)*p_y*(p_t-1)+ul(1,1),3));
%                 z_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),1)=datai_z(1,se(j+(k-1)*p_y*(p_t-1)+ul(1,1),2));
%                 z_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),2)=datai_z(1,se(j+(k-1)*p_y*(p_t-1)+ul(1,1),3));
%                 plot3(x_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),:)*1000,y_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),:)*1000,z_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),:)*1000,'k');
%                 hold on
%             end
%         end
for k=1:p_h
    for j=1:p_y:1+(p_t-2)*p_y
        x_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),1)=datai_x(1,se(j+(k-1)*p_y*(p_t-1)+ul(1,1),2));
        x_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),2)=datai_x(1,se(j+(k-1)*p_y*(p_t-1)+ul(1,1),3));
        y_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),1)=datai_y(1,se(j+(k-1)*p_y*(p_t-1)+ul(1,1),2));
        y_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),2)=datai_y(1,se(j+(k-1)*p_y*(p_t-1)+ul(1,1),3));
        z_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),1)=datai_z(1,se(j+(k-1)*p_y*(p_t-1)+ul(1,1),2));
        z_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),2)=datai_z(1,se(j+(k-1)*p_y*(p_t-1)+ul(1,1),3));
        plot3(x_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),:)*1000,y_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),:)*1000,z_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),:)*1000,'k');
        hold on
    end
    for j=p_y:p_y:p_y+(p_t-2)*p_y
        x_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),1)=datai_x(1,se(j+(k-1)*p_y*(p_t-1)+ul(1,1),2));
        x_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),2)=datai_x(1,se(j+(k-1)*p_y*(p_t-1)+ul(1,1),3));
        y_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),1)=datai_y(1,se(j+(k-1)*p_y*(p_t-1)+ul(1,1),2));
        y_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),2)=datai_y(1,se(j+(k-1)*p_y*(p_t-1)+ul(1,1),3));
        z_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),1)=datai_z(1,se(j+(k-1)*p_y*(p_t-1)+ul(1,1),2));
        z_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),2)=datai_z(1,se(j+(k-1)*p_y*(p_t-1)+ul(1,1),3));
        plot3(x_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),:)*1000,y_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),:)*1000,z_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),:)*1000,'k');
        hold on
    end
end


%%高さ線
for k=0:p_y*p_t:(p_h-2)*p_y*p_t
    for j=1:p_y
        x_high(j+k+ul(1,2),1)=datai_x(1,se(j+k+ul(1,2),2));
        x_high(j+k+ul(1,2),2)=datai_x(1,se(j+k+ul(1,2),3));
        y_high(j+k+ul(1,2),1)=datai_y(1,se(j+k+ul(1,2),2));
        y_high(j+k+ul(1,2),2)=datai_y(1,se(j+k+ul(1,2),3));
        z_high(j+k+ul(1,2),1)=datai_z(1,se(j+k+ul(1,2),2));
        z_high(j+k+ul(1,2),2)=datai_z(1,se(j+k+ul(1,2),3));
        % if max(fiber_force(1,:))==0
        %     color(j+k,1)=0;
        % else
        % color(j+k,1)=fiber_force(1,j+k)/max(fiber_force(1,:))*255/255;
        % end
        plot3(x_high(j+k+ul(1,2),:)*1000,y_high(j+k+ul(1,2),:)*1000,z_high(j+k+ul(1,2),:)*1000,'r');
        hold on
    end
    for j=p_y*(p_t-1)+1:p_y*p_t
        x_high(j+k+ul(1,2),1)=datai_x(1,se(j+k+ul(1,2),2));
        x_high(j+k+ul(1,2),2)=datai_x(1,se(j+k+ul(1,2),3));
        y_high(j+k+ul(1,2),1)=datai_y(1,se(j+k+ul(1,2),2));
        y_high(j+k+ul(1,2),2)=datai_y(1,se(j+k+ul(1,2),3));
        z_high(j+k+ul(1,2),1)=datai_z(1,se(j+k+ul(1,2),2));
        z_high(j+k+ul(1,2),2)=datai_z(1,se(j+k+ul(1,2),3));
        % if max(fiber_force(1,:))==0
        %     color(j+k,1)=0;
        % else
        % color(j+k,1)=fiber_force(1,j+k)/max(fiber_force(1,:))*255/255;
        % end
        plot3(x_high(j+k+ul(1,2),:)*1000,y_high(j+k+ul(1,2),:)*1000,z_high(j+k+ul(1,2),:)*1000,'r');
        hold on
    end
    for j=p_y+1:p_y:1+p_y*(p_t-2)
        x_high(j+k+ul(1,2),1)=datai_x(1,se(j+k+ul(1,2),2));
        x_high(j+k+ul(1,2),2)=datai_x(1,se(j+k+ul(1,2),3));
        y_high(j+k+ul(1,2),1)=datai_y(1,se(j+k+ul(1,2),2));
        y_high(j+k+ul(1,2),2)=datai_y(1,se(j+k+ul(1,2),3));
        z_high(j+k+ul(1,2),1)=datai_z(1,se(j+k+ul(1,2),2));
        z_high(j+k+ul(1,2),2)=datai_z(1,se(j+k+ul(1,2),3));
        % if max(fiber_force(1,:))==0
        %     color(j+k,1)=0;
        % else
        % color(j+k,1)=fiber_force(1,j+k)/max(fiber_force(1,:))*255/255;
        % end
        plot3(x_high(j+k+ul(1,2),:)*1000,y_high(j+k+ul(1,2),:)*1000,z_high(j+k+ul(1,2),:)*1000,'r');
        hold on
    end
    for j=p_y:p_y:p_y*(p_t-1)
        x_high(j+k+ul(1,2),1)=datai_x(1,se(j+k+ul(1,2),2));
        x_high(j+k+ul(1,2),2)=datai_x(1,se(j+k+ul(1,2),3));
        y_high(j+k+ul(1,2),1)=datai_y(1,se(j+k+ul(1,2),2));
        y_high(j+k+ul(1,2),2)=datai_y(1,se(j+k+ul(1,2),3));
        z_high(j+k+ul(1,2),1)=datai_z(1,se(j+k+ul(1,2),2));
        z_high(j+k+ul(1,2),2)=datai_z(1,se(j+k+ul(1,2),3));
        % if max(fiber_force(1,:))==0
        %     color(j+k,1)=0;
        % else
        % color(j+k,1)=fiber_force(1,j+k)/max(fiber_force(1,:))*255/255;
        % end
        plot3(x_high(j+k+ul(1,2),:)*1000,y_high(j+k+ul(1,2),:)*1000,z_high(j+k+ul(1,2),:)*1000,'r');
        hold on
    end
end
xlabel('x')
ylabel('y')
zlabel('z')
axis equal
%         %         xlim([-0.1*(10^8),0.25*(10^8)]);
%            xlim([-150000,350000]);
%          xlim([-0.2*(10^8),0.25*(10^8)]);
%         % %         ylim([14.15*(10^8),14.50*(10^8)]);
%            ylim([-150000,350000]);
%          ylim([-0.2*(10^8),0.25*(10^8)]);
%         % %         zlim([1.20*(10^8),1.55*(10^8)]);
%            zlim([-200000,1270]);
% zlim([1060, 1260]);
% Frame(l) = getframe(1);
hold off