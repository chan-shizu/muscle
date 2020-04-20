clc;
clear;
muscleNames = loadMuscleName();
muscle_name = [muscleNames{1}]
%se= csvread('se_arm.csv', 0, 0);
se_name = ['C:\Users\kou_0\OneDrive\ドキュメント\研究matlab\parameter\',muscle_name, '_se.csv']
se= csvread(se_name, 0, 0);
seNum=size(se);
file_name_data = ['C:\Users\kou_0\OneDrive\ドキュメント\研究matlab\muscle\data_',muscle_name, '.csv']%"data_rot_h.csv"
data= csvread(file_name_data, 0, 0);
timeNum = size(data);
p_y=5;
p_t=5;
p_h=6;

% fiber_force=csvread('fiber_force.csv', 0, 0);
data_x_name = ['output\', muscle_name,'_data_x.csv']
data_y_name = ['output\', muscle_name,'_data_y.csv']
data_z_name = ['output\', muscle_name,'_data_z.csv']
datai_x=csvread(data_x_name, 0, 0);
datai_y=csvread(data_y_name, 0, 0);
datai_z=csvread(data_z_name, 0, 0);

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

for i=1:timeNum(1)
    if i==1||mod(i,10)==0
        if i==1
            l=1;
        else
            l=i/10;
        end
        
        
        %%横線
        %伸ばす線がないnode
        for j=1:p_h
            for k=1:p_y-1;
                x_yoko(k+(j-1)*(p_y-1)*p_t,1)=datai_x(i,se(k+(j-1)*(p_y-1)*p_t,2));
                x_yoko(k+(j-1)*(p_y-1)*p_t,2)=datai_x(i,se(k+(j-1)*(p_y-1)*p_t,3));
                y_yoko(k+(j-1)*(p_y-1)*p_t,1)=datai_y(i,se(k+(j-1)*(p_y-1)*p_t,2));
                y_yoko(k+(j-1)*(p_y-1)*p_t,2)=datai_y(i,se(k+(j-1)*(p_y-1)*p_t,3));
                z_yoko(k+(j-1)*(p_y-1)*p_t,1)=datai_z(i,se(k+(j-1)*(p_y-1)*p_t,2));
                z_yoko(k+(j-1)*(p_y-1)*p_t,2)=datai_z(i,se(k+(j-1)*(p_y-1)*p_t,3));
                plot3(x_yoko(k+(j-1)*(p_y-1)*p_t,:)*1000,y_yoko(k+(j-1)*(p_y-1)*p_t,:)*1000,z_yoko(k+(j-1)*(p_y-1)*p_t,:)*1000,'k');
                hold on
            end
            for k=1+(p_y-1)*(p_t-1):(p_y-1)*p_t;
                x_yoko(k+(j-1)*(p_y-1)*p_t,1)=datai_x(i,se(k+(j-1)*(p_y-1)*p_t,2));
                x_yoko(k+(j-1)*(p_y-1)*p_t,2)=datai_x(i,se(k+(j-1)*(p_y-1)*p_t,3));
                y_yoko(k+(j-1)*(p_y-1)*p_t,1)=datai_y(i,se(k+(j-1)*(p_y-1)*p_t,2));
                y_yoko(k+(j-1)*(p_y-1)*p_t,2)=datai_y(i,se(k+(j-1)*(p_y-1)*p_t,3));
                z_yoko(k+(j-1)*(p_y-1)*p_t,1)=datai_z(i,se(k+(j-1)*(p_y-1)*p_t,2));
                z_yoko(k+(j-1)*(p_y-1)*p_t,2)=datai_z(i,se(k+(j-1)*(p_y-1)*p_t,3));
                plot3(x_yoko(k+(j-1)*(p_y-1)*p_t,:)*1000,y_yoko(k+(j-1)*(p_y-1)*p_t,:)*1000,z_yoko(k+(j-1)*(p_y-1)*p_t,:)*1000,'k');
                hold on
            end
        end
        
        
        
        %%縦線
        %伸ばす線がないnode
        for k=1:p_h;
            for j=1:p_y:1+(p_t-2)*p_y;
                x_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),1)=datai_x(i,se(j+(k-1)*p_y*(p_t-1)+ul(1,1),2));
                x_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),2)=datai_x(i,se(j+(k-1)*p_y*(p_t-1)+ul(1,1),3));
                y_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),1)=datai_y(i,se(j+(k-1)*p_y*(p_t-1)+ul(1,1),2));
                y_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),2)=datai_y(i,se(j+(k-1)*p_y*(p_t-1)+ul(1,1),3));
                z_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),1)=datai_z(i,se(j+(k-1)*p_y*(p_t-1)+ul(1,1),2));
                z_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),2)=datai_z(i,se(j+(k-1)*p_y*(p_t-1)+ul(1,1),3));
                plot3(x_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),:)*1000,y_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),:)*1000,z_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),:)*1000,'k');
                hold on
            end
            for j=p_y:p_y:p_y+(p_t-2)*p_y;
                x_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),1)=datai_x(i,se(j+(k-1)*p_y*(p_t-1)+ul(1,1),2));
                x_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),2)=datai_x(i,se(j+(k-1)*p_y*(p_t-1)+ul(1,1),3));
                y_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),1)=datai_y(i,se(j+(k-1)*p_y*(p_t-1)+ul(1,1),2));
                y_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),2)=datai_y(i,se(j+(k-1)*p_y*(p_t-1)+ul(1,1),3));
                z_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),1)=datai_z(i,se(j+(k-1)*p_y*(p_t-1)+ul(1,1),2));
                z_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),2)=datai_z(i,se(j+(k-1)*p_y*(p_t-1)+ul(1,1),3));
                plot3(x_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),:)*1000,y_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),:)*1000,z_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),:)*1000,'k');
                hold on
            end
        end
        
        
        %%高さ線
        for k=0:p_y*p_t:(p_h-2)*p_y*p_t;
            for j=1:p_y;
                x_high(j+k+ul(1,2),1)=datai_x(i,se(j+k+ul(1,2),2));
                x_high(j+k+ul(1,2),2)=datai_x(i,se(j+k+ul(1,2),3));
                y_high(j+k+ul(1,2),1)=datai_y(i,se(j+k+ul(1,2),2));
                y_high(j+k+ul(1,2),2)=datai_y(i,se(j+k+ul(1,2),3));
                z_high(j+k+ul(1,2),1)=datai_z(i,se(j+k+ul(1,2),2));
                z_high(j+k+ul(1,2),2)=datai_z(i,se(j+k+ul(1,2),3));
                % if max(fiber_force(i,:))==0
                %     color(j+k,1)=0;
                % else
                % color(j+k,1)=fiber_force(i,j+k)/max(fiber_force(i,:))*255/255;
                % end
                plot3(x_high(j+k+ul(1,2),:)*1000,y_high(j+k+ul(1,2),:)*1000,z_high(j+k+ul(1,2),:)*1000,'r');
                hold on
            end
            for j=p_y*(p_t-1)+1:p_y*p_t;
                x_high(j+k+ul(1,2),1)=datai_x(i,se(j+k+ul(1,2),2));
                x_high(j+k+ul(1,2),2)=datai_x(i,se(j+k+ul(1,2),3));
                y_high(j+k+ul(1,2),1)=datai_y(i,se(j+k+ul(1,2),2));
                y_high(j+k+ul(1,2),2)=datai_y(i,se(j+k+ul(1,2),3));
                z_high(j+k+ul(1,2),1)=datai_z(i,se(j+k+ul(1,2),2));
                z_high(j+k+ul(1,2),2)=datai_z(i,se(j+k+ul(1,2),3));
                % if max(fiber_force(i,:))==0
                %     color(j+k,1)=0;
                % else
                % color(j+k,1)=fiber_force(i,j+k)/max(fiber_force(i,:))*255/255;
                % end
                plot3(x_high(j+k+ul(1,2),:)*1000,y_high(j+k+ul(1,2),:)*1000,z_high(j+k+ul(1,2),:)*1000,'r');
                hold on
            end
            for j=p_y+1:p_y:1+p_y*(p_t-2);
                x_high(j+k+ul(1,2),1)=datai_x(i,se(j+k+ul(1,2),2));
                x_high(j+k+ul(1,2),2)=datai_x(i,se(j+k+ul(1,2),3));
                y_high(j+k+ul(1,2),1)=datai_y(i,se(j+k+ul(1,2),2));
                y_high(j+k+ul(1,2),2)=datai_y(i,se(j+k+ul(1,2),3));
                z_high(j+k+ul(1,2),1)=datai_z(i,se(j+k+ul(1,2),2));
                z_high(j+k+ul(1,2),2)=datai_z(i,se(j+k+ul(1,2),3));
                % if max(fiber_force(i,:))==0
                %     color(j+k,1)=0;
                % else
                % color(j+k,1)=fiber_force(i,j+k)/max(fiber_force(i,:))*255/255;
                % end
                plot3(x_high(j+k+ul(1,2),:)*1000,y_high(j+k+ul(1,2),:)*1000,z_high(j+k+ul(1,2),:)*1000,'r');
                hold on
            end
            for j=p_y:p_y:p_y*(p_t-1);
                x_high(j+k+ul(1,2),1)=datai_x(i,se(j+k+ul(1,2),2));
                x_high(j+k+ul(1,2),2)=datai_x(i,se(j+k+ul(1,2),3));
                y_high(j+k+ul(1,2),1)=datai_y(i,se(j+k+ul(1,2),2));
                y_high(j+k+ul(1,2),2)=datai_y(i,se(j+k+ul(1,2),3));
                z_high(j+k+ul(1,2),1)=datai_z(i,se(j+k+ul(1,2),2));
                z_high(j+k+ul(1,2),2)=datai_z(i,se(j+k+ul(1,2),3));
                % if max(fiber_force(i,:))==0
                %     color(j+k,1)=0;
                % else
                % color(j+k,1)=fiber_force(i,j+k)/max(fiber_force(i,:))*255/255;
                % end
                plot3(x_high(j+k+ul(1,2),:)*1000,y_high(j+k+ul(1,2),:)*1000,z_high(j+k+ul(1,2),:)*1000,'r');
                hold on
            end
        end
        xlabel('x')
        ylabel('y')
        zlabel('z')
        %       xlim([82.5,282.5]);
        xlim([-1.58*(10^20),-1.53*(10^20)]);
        %         % % ylim([-192.5,7.5]);
        ylim([[2.45*(10^20),2.50*(10^20)]]);
        %         % %zlim([1055,1255]);
        zlim([0.23*(10^20),0.28*(10^20)]);
        Frame(l) = getframe(1);
        hold off
    end
end

% v = VideoWriter('speedup');
videoName = ['\video\',muscle_name, '_upsideDown']
v = VideoWriter(videoName);
v.Quality=50;

open(v);
writeVideo(v,Frame);
close(v);
