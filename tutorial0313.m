% clc;
% clear;
warning('on','all');
warning;

muscleNames = loadMuscleName();
muscle_name = [muscleNames{1}]

% Fn_x(i,:)=FsSum_x(i,:)+Fv_x(i,:)+Mr_x(i,:)+con_x(i,:)+vis_x(i,:);
% Fn_y(i,:)=FsSum_y(i,:)+Fv_y(i,:)+Mr_y(i,:)+con_y(i,:)+vis_y(i,:);
% Fn_z(i,:)=FsSum_z(i,:)+Fv_z(i,:)+Mr_z(i,:)+con_z(i,:)+vis_z(i,:);
% v = VideoWriter('test');
% videoName = ['video\',muscle_name, '_volume_activr=0.1_dt=10^-4_c=600_kv300',];
% v = VideoWriter(videoName);
% v.Quality=50;

% open(v);
% writeVideo(v,Frame);
% close(v);

FsSumTotal = abs(FsSum_x) +  abs(FsSum_y) + abs(FsSum_z);
FvTotal = abs(Fv_x) + abs(Fv_y) + abs(Fv_z);
MrTotal = abs(Mr_x) + abs(Mr_y) + abs(Mr_z);
conTotal = abs(con_x) + abs(con_y) + abs(con_z);
visTotal = abs(vis_x) + abs(vis_y) + abs(vis_z);
% FgTotal = abs(Fg_x) + abs(Fg_y) + abs(Fg_z);
HillActive = abs(HillActive);
HillPassive = abs(HillPassive);

FsSumAve = mean(FsSumTotal');
FvAve = mean(FvTotal');
MrAve = mean(MrTotal');
conAve = mean(conTotal');
visAve = mean(visTotal');
% FgTotal = mean(FgTotal');
HillActiveAve = mean(HillActive')*1000;
HillPassiveAve = mean(HillPassive')*1000;

plot(FsSumAve,'r');
hold on
plot(FvAve,'k');
plot(MrAve,'g');
plot(conAve,'b');
plot(visAve,'c');
% plot(FgTotal,'c-');
plot(HillActiveAve,'m');
plot(HillPassiveAve,'y');
legend('FsSum','Fv','Mr','con','vis','HillActive','HillPassive');
xlabel('Time step');
ylabel('Force');
hold off
% springkArray = 1000*ones(1,957);
% fileNameSprigk = ['C:\Users\kou_0\OneDrive\ドキュメント\研究matlab\parameter\',muscle_name, '_springk.csv'];
% csvwrite(fileNameSprigk, springkArray);
%% 座標表示
% muscleNames = loadMuscleName();
% muscle_name = [muscleNames{1}]
% file_name_readcsv = ['C:\Users\kou_0\OneDrive\ドキュメント\研究matlab\muscle\', muscle_name, '_min_final.csv']
% % file_name_readcsv = strcat('C:\Users\bubbl\Documents\shizuya_M1\DefMuscle_for_TUS_2018\matlab\surface2grid_0608\muscle\', muscle_name, '_min_final.csv')
% file = csvread(file_name_readcsv);
% plot3(file(:,1), file(:,2), file(:,3), 'b.' );
% xlabel('x')
% ylabel('y')
% zlabel('z')
% axis equal
% 
% divisionMuscle = readmatrix("C:\Users\kou_0\OneDrive\ドキュメント\研究matlab\divisionMuscle.csv")
% p_y=divisionMuscle(1);
% p_t=divisionMuscle(2);
% p_h=divisionMuscle(3);
% 
% se_name = ['C:\Users\kou_0\OneDrive\ドキュメント\研究matlab\parameter\',muscle_name, '_se.csv']
% se= csvread(se_name, 0, 0);
% seNum=size(se);
% 
% for n=1:6;
%     for j=1:seNum(1);
%         if n==6;
%             ul(1,n)=seNum(1);
%         end
%         if se(j,1)>n-1;
%             ul(1,n)=j-1;  %upper limit
%             break
%         end
%     end
% end
% 
% file_name_readcsv = ['C:\Users\kou_0\OneDrive\ドキュメント\研究matlab\muscle\', muscle_name, '_min_final.csv']
% file = csvread(file_name_readcsv);
% plot3(file(:,1)*1000, file(:,2)*1000, file(:,3)*1000, 'k.' );
% hold on
% 
% for j=1:p_h
%     %それぞれの層の1番目と2番目の質点をつなぐ(5×5×6の場合1と2, 25と26など)
%     for k=1:p_y-1
%         x_yoko(k+(j-1)*(p_y-1)*p_t,1)=file(se(k+(j-1)*(p_y-1)*p_t,2),1);
%         x_yoko(k+(j-1)*(p_y-1)*p_t,2)=file(se(k+(j-1)*(p_y-1)*p_t,3),1);
%         y_yoko(k+(j-1)*(p_y-1)*p_t,1)=file(se(k+(j-1)*(p_y-1)*p_t,2),2);
%         y_yoko(k+(j-1)*(p_y-1)*p_t,2)=file(se(k+(j-1)*(p_y-1)*p_t,3),2);
%         z_yoko(k+(j-1)*(p_y-1)*p_t,1)=file(se(k+(j-1)*(p_y-1)*p_t,2),3);
%         z_yoko(k+(j-1)*(p_y-1)*p_t,2)=file(se(k+(j-1)*(p_y-1)*p_t,3),3);
%         plot3(x_yoko(k+(j-1)*(p_y-1)*p_t,:)*1000,y_yoko(k+(j-1)*(p_y-1)*p_t,:)*1000,z_yoko(k+(j-1)*(p_y-1)*p_t,:)*1000,'k');
%         hold on
%     end
%     for k=1+(p_y-1)*(p_t-1):(p_y-1)*p_t
%         x_yoko(k+(j-1)*(p_y-1)*p_t,1)=file(se(k+(j-1)*(p_y-1)*p_t,2),1);
%         x_yoko(k+(j-1)*(p_y-1)*p_t,2)=file(se(k+(j-1)*(p_y-1)*p_t,3),1);
%         y_yoko(k+(j-1)*(p_y-1)*p_t,1)=file(se(k+(j-1)*(p_y-1)*p_t,2),2);
%         y_yoko(k+(j-1)*(p_y-1)*p_t,2)=file(se(k+(j-1)*(p_y-1)*p_t,3),2);
%         z_yoko(k+(j-1)*(p_y-1)*p_t,1)=file(se(k+(j-1)*(p_y-1)*p_t,2),3);
%         z_yoko(k+(j-1)*(p_y-1)*p_t,2)=file(se(k+(j-1)*(p_y-1)*p_t,3),3);
%         plot3(x_yoko(k+(j-1)*(p_y-1)*p_t,:)*1000,y_yoko(k+(j-1)*(p_y-1)*p_t,:)*1000,z_yoko(k+(j-1)*(p_y-1)*p_t,:)*1000,'k');
%         hold on
%     end
% end
% 
% for k=1:p_h;
%     for j=1:p_y:1+(p_t-2)*p_y;
%         x_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),1)=file(se(j+(k-1)*p_y*(p_t-1)+ul(1,1),2),1);
%         x_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),2)=file(se(j+(k-1)*p_y*(p_t-1)+ul(1,1),3),1);
%         y_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),1)=file(se(j+(k-1)*p_y*(p_t-1)+ul(1,1),2),2);
%         y_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),2)=file(se(j+(k-1)*p_y*(p_t-1)+ul(1,1),3),2);
%         z_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),1)=file(se(j+(k-1)*p_y*(p_t-1)+ul(1,1),2),3);
%         z_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),2)=file(se(j+(k-1)*p_y*(p_t-1)+ul(1,1),3),3);
%         plot3(x_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),:)*1000,y_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),:)*1000,z_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),:)*1000,'b');
%         hold on
%     end
%     for j=p_y:p_y:p_y+(p_t-2)*p_y;
%         x_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),1)=file(se(j+(k-1)*p_y*(p_t-1)+ul(1,1),2),1);
%         x_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),2)=file(se(j+(k-1)*p_y*(p_t-1)+ul(1,1),3),1);
%         y_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),1)=file(se(j+(k-1)*p_y*(p_t-1)+ul(1,1),2),2);
%         y_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),2)=file(se(j+(k-1)*p_y*(p_t-1)+ul(1,1),3),2);
%         z_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),1)=file(se(j+(k-1)*p_y*(p_t-1)+ul(1,1),2),3);
%         z_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),2)=file(se(j+(k-1)*p_y*(p_t-1)+ul(1,1),3),3);
%         plot3(x_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),:)*1000,y_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),:)*1000,z_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),:)*1000,'b');
%         hold on
%     end
% end
% 
% %高さ線
% for k=0:p_y*p_t:(p_h-2)*p_y*p_t;
%     for j=1:p_y;
%         x_high(j+k+ul(1,2),1)=file(se(j+k+ul(1,2),2),1);
%         x_high(j+k+ul(1,2),2)=file(se(j+k+ul(1,2),3),1);
%         y_high(j+k+ul(1,2),1)=file(se(j+k+ul(1,2),2),2);
%         y_high(j+k+ul(1,2),2)=file(se(j+k+ul(1,2),3),2);
%         z_high(j+k+ul(1,2),1)=file(se(j+k+ul(1,2),2),3);
%         z_high(j+k+ul(1,2),2)=file(se(j+k+ul(1,2),3),3);
%         plot3(x_high(j+k+ul(1,2),:)*1000,y_high(j+k+ul(1,2),:)*1000,z_high(j+k+ul(1,2),:)*1000,'r');
%         hold on
%     end
%     for j=p_y*(p_t-1)+1:p_y*p_t;
%         x_high(j+k+ul(1,2),1)=file(se(j+k+ul(1,2),2),1);
%         x_high(j+k+ul(1,2),2)=file(se(j+k+ul(1,2),3),1);
%         y_high(j+k+ul(1,2),1)=file(se(j+k+ul(1,2),2),2);
%         y_high(j+k+ul(1,2),2)=file(se(j+k+ul(1,2),3),2);
%         z_high(j+k+ul(1,2),1)=file(se(j+k+ul(1,2),2),3);
%         z_high(j+k+ul(1,2),2)=file(se(j+k+ul(1,2),3),3);
%         plot3(x_high(j+k+ul(1,2),:)*1000,y_high(j+k+ul(1,2),:)*1000,z_high(j+k+ul(1,2),:)*1000,'r');
%         hold on
%     end
%     for j=p_y+1:p_y:1+p_y*(p_t-2);
%         x_high(j+k+ul(1,2),1)=file(se(j+k+ul(1,2),2),1);
%         x_high(j+k+ul(1,2),2)=file(se(j+k+ul(1,2),3),1);
%         y_high(j+k+ul(1,2),1)=file(se(j+k+ul(1,2),2),2);
%         y_high(j+k+ul(1,2),2)=file(se(j+k+ul(1,2),3),2);
%         z_high(j+k+ul(1,2),1)=file(se(j+k+ul(1,2),2),3);
%         z_high(j+k+ul(1,2),2)=file(se(j+k+ul(1,2),3),3);
%         plot3(x_high(j+k+ul(1,2),:)*1000,y_high(j+k+ul(1,2),:)*1000,z_high(j+k+ul(1,2),:)*1000,'r');
%         hold on
%     end
%     for j=p_y:p_y:p_y*(p_t-1);
%         x_high(j+k+ul(1,2),1)=file(se(j+k+ul(1,2),2),1);
%         x_high(j+k+ul(1,2),2)=file(se(j+k+ul(1,2),3),1);
%         y_high(j+k+ul(1,2),1)=file(se(j+k+ul(1,2),2),2);
%         y_high(j+k+ul(1,2),2)=file(se(j+k+ul(1,2),3),2);
%         z_high(j+k+ul(1,2),1)=file(se(j+k+ul(1,2),2),3);
%         z_high(j+k+ul(1,2),2)=file(se(j+k+ul(1,2),3),3);
%         plot3(x_high(j+k+ul(1,2),:)*1000,y_high(j+k+ul(1,2),:)*1000,z_high(j+k+ul(1,2),:)*1000,'r');
%         hold on
%     end
% end
% xlabel('x')
% ylabel('y')
% zlabel('z')
% % %         xlim([-0.1*(10^8),0.25*(10^8)]);
% % %         xlim([-1.58*(10^20),-1.53*(10^20)]);
% % xlim([-0.2*(10^8),0.25*(10^8)]);
% % % %         ylim([14.15*(10^8),14.50*(10^8)]);
% % %         ylim([[2.45*(10^20),2.50*(10^20)]]);
% % ylim([-0.2*(10^8),0.25*(10^8)]);
% % % %         zlim([1.20*(10^8),1.55*(10^8)]);
% % %         zlim([0.23*(10^20),0.28*(10^20)]);
% % zlim([0.95*(10^8),1.40*(10^8)]);
% % Frame(l) = getframe(1);
% axis equal
% hold off