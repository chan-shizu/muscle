clc;
clear;
muscleNames = loadMuscleName();
prompt = 'è„òrìÒì™ãÿÇ»ÇÁ1, ÉIÉgÉKÉCê„çúãÿÇ»ÇÁ2, åsìÀê„çúãÿÇ»ÇÁ3, ä{ìÒï†ãÿëOï†Ç»ÇÁ4, ä{ìÒï†ãÿå„ï†Ç»ÇÁ5ÇâüÇµÇƒÇÀ';
muscleNumber = inputdlg(prompt,...
             'choose muscle', [1 50])
muscleNumber = str2num(muscleNumber{1})
muscle_name = [muscleNames{muscleNumber,1}];
%se= csvread('se_arm.csv', 0, 0);
se_name = ['parameter\',muscle_name, '_se.csv'];
se= csvread(se_name, 0, 0);
seNum=size(se);
file_name_data = ['muscle\data_',muscle_name, '.csv'];%"data_rot_h.csv"
data= csvread(file_name_data, 0, 0);
timeNum = size(data);
divisionMuscle = readmatrix("divisionMuscle.csv");
p_y=divisionMuscle(1);
p_t=divisionMuscle(2);
p_h=divisionMuscle(3);

% fiber_force=csvread('fiber_force.csv', 0, 0);
data_x_name = ['output\', muscle_name,'_data_x.csv'];
data_y_name = ['output\', muscle_name,'_data_y.csv'];
data_z_name = ['output\', muscle_name,'_data_z.csv'];
datai_x=csvread(data_x_name, 0, 0);
datai_y=csvread(data_y_name, 0, 0);
datai_z=csvread(data_z_name, 0, 0);

checkRotation = false;%true or forse

if checkRotation == 1
    %% zé≤Ç‹ÇÌÇËÇ…âÒì]
    alpha = pi()/6;
    beta = pi()/3;
    convertData = zeros(p_t*p_y,3*p_h);
    
    average_lower_x = mean(datai_x(1,1:p_t*p_y));
    average_lower_y = mean(datai_y(1,1:p_t*p_y));
    average_lower_z = mean(datai_z(1,1:p_t*p_y));
    center_of_rotation = [1.0*average_lower_x,...
        1.0*average_lower_y,...
        1.0*average_lower_z];
    %
    Coordinate_transformed_z = [cos(alpha) sin(alpha) 0 ; -sin(alpha) cos(alpha) 0; 0 0 1];%zé≤Ç‹ÇÌÇËÇÃâÒì]
    Coordinate_transformed_y = [cos(pi/2-beta) 0 -sin(pi/2-beta);0 1 0;sin(pi/2-beta) 0 cos(pi/2-beta)];%yé≤Ç‹ÇÌÇËÇÃâÒì]
    Coordinate_transformed_file = zeros(p_y*p_t,p_h*3);
    
    for i = 1:timeNum(1)
        i
        for n=1:p_h
            convertData(1:p_t*p_y,3*(n-1)+1) = datai_x(i,(n-1)*(p_t*p_y)+1:n*(p_t*p_y));
            convertData(1:(p_t*p_y),3*(n-1)+2) = datai_y(i,(n-1)*(p_t*p_y)+1:n*(p_t*p_y));
            convertData(1:(p_t*p_y),3*(n-1)+3) = datai_z(i,(n-1)*(p_t*p_y)+1:n*(p_t*p_y));
        end
        
        
        
        for n=1:p_t*p_y
            for j=1:p_h
                x_y_z = convertData(n,1+3*(j-1):3*j) - center_of_rotation;
                Coordinate_transformed_file(n,1+3*(j-1):3*j) = Coordinate_transformed_z * x_y_z.' + center_of_rotation.';
            end
        end
        
        %         xx = zeros(1, p_t*p_y*p_h);
        %         yy= zeros(1, p_t*p_y*p_h);
        %         zz = zeros(1, p_t*p_y*p_h);
        %         xx_2 = zeros(1, p_t*p_y*p_h);
        %         yy_2 = zeros(1, p_t*p_y*p_h);
        %         zz_2 = zeros(1, p_t*p_y*p_h);
        %
        %         for i=1:p_t*p_y
        %             for j=1:p_h
        %                 xx(p_h*(i-1)+j) = convertData(i,3*(j-1)+1);
        %             end
        %         end
        %         for i=1:p_t*p_y
        %             for j=1:p_h
        %                 yy(p_h*(i-1)+j) = convertData(i,3*(j-1)+2);
        %             end
        %         end
        %         for i=1:p_t*p_y
        %             for j=1:p_h
        %                 zz(p_h*(i-1)+j) = convertData(i,3*(j-1)+3);
        %             end
        %         end
        %
        %         for i=1:p_t*p_y
        %             for j=1:p_h
        %                 xx_2(p_h*(i-1)+j) = Coordinate_transformed_file(i,3*(j-1)+1);
        %             end
        %         end
        %         for i=1:p_t*p_y
        %             for j=1:p_h
        %                 yy_2(p_h*(i-1)+j) = Coordinate_transformed_file(i,3*(j-1)+2);
        %             end
        %         end
        %         for i=1:p_t*p_y
        %             for j=1:p_h
        %                 zz_2(p_h*(i-1)+j) = Coordinate_transformed_file(i,3*(j-1)+3);
        %             end
        %         end
        %
        %         f1 = figure;
        %         plot3(xx(1,:), yy(1,:), zz(1,:), 'b.' );
        %         xlabel('x')
        %         ylabel('y')
        %         zlabel('z')
        %
        %         f2 = figure;
        %         plot3(xx_2(1,:), yy_2(1,:), zz_2(1,:), 'b.' );
        %         xlabel('x')
        %         ylabel('y')
        %         zlabel('z')
        %         axis equal
        
        
        %% yé≤Ç‹ÇÌÇËÇ…âÒì]
        for n=1:p_t*p_y
            for j=1:p_h
                x_y_z = Coordinate_transformed_file(n,1+3*(j-1):3*j) - center_of_rotation;
                Coordinate_transformed_file(n,1+3*(j-1):3*j) = Coordinate_transformed_y * x_y_z.' + center_of_rotation.';
            end
        end
        
        %         for i=1:p_t*p_y
        %             for j=1:p_h
        %                 xx3(p_h*(i-1)+j) = Coordinate_transformed_file(i,3*(j-1)+1);
        %             end
        %         end
        %         for i=1:p_t*p_y
        %             for j=1:p_h
        %                 yy3(p_h*(i-1)+j) = Coordinate_transformed_file(i,3*(j-1)+2);
        %             end
        %         end
        %         for i=1:p_t*p_y
        %             for j=1:p_h
        %                 zz3(p_h*(i-1)+j) = Coordinate_transformed_file(i,3*(j-1)+3);
        %             end
        %         end
        %
        %         f3 = figure;
        %         plot3(xx3(1,:), yy3(1,:), zz3(1,:), 'b.' );
        %         xlabel('x')
        %         ylabel('y')
        %         zlabel('z')
        %         axis equal
        
        for n = 1:p_h
            datai_x(i,(n-1)*p_y*p_t+1:n*p_y*p_t) = Coordinate_transformed_file(1:p_y*p_t,(n-1)*3+1);
            datai_y(i,(n-1)*p_y*p_t+1:n*p_y*p_t) = Coordinate_transformed_file(1:p_y*p_t,(n-1)*3+2);
            datai_z(i,(n-1)*p_y*p_t+1:n*p_y*p_t) = Coordinate_transformed_file(1:p_y*p_t,(n-1)*3+3);
        end
    end
    %     CoordinateTransformedFileTable = array2table(Coordinate_transformed_file);
    %     fileText = fileText(:, 1:3);
    %     CoordinateTransformedFileTable = [fileText CoordinateTransformedFileTable];
    %     writetable(CoordinateTransformedFileTable, file_name_readcsv, 'WriteVariableNames',false);
end

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

for i=1:timeNum(1)
    if i==1||mod(i,2)==0%10
        if i==1
            l=1;
        else
            l=i/2;%10
        end
        
        
        %%â°ê¸
        %êLÇŒÇ∑ê¸Ç™Ç»Ç¢node
        for j=1:p_h
            %ÇªÇÍÇºÇÍÇÃëwÇÃ1î‘ñ⁄Ç∆2î‘ñ⁄ÇÃéøì_ÇÇ¬Ç»ÇÆ(5Å~5Å~6ÇÃèÍçá1Ç∆2, 25Ç∆26Ç»Ç«)
            for k=1:p_y-1
                x_yoko(k+(j-1)*(p_y-1)*p_t,1)=datai_x(i,se(k+(j-1)*(p_y-1)*p_t,2));
                x_yoko(k+(j-1)*(p_y-1)*p_t,2)=datai_x(i,se(k+(j-1)*(p_y-1)*p_t,3));
                y_yoko(k+(j-1)*(p_y-1)*p_t,1)=datai_y(i,se(k+(j-1)*(p_y-1)*p_t,2));
                y_yoko(k+(j-1)*(p_y-1)*p_t,2)=datai_y(i,se(k+(j-1)*(p_y-1)*p_t,3));
                z_yoko(k+(j-1)*(p_y-1)*p_t,1)=datai_z(i,se(k+(j-1)*(p_y-1)*p_t,2));
                z_yoko(k+(j-1)*(p_y-1)*p_t,2)=datai_z(i,se(k+(j-1)*(p_y-1)*p_t,3));
                plot3(x_yoko(k+(j-1)*(p_y-1)*p_t,:)*1000,y_yoko(k+(j-1)*(p_y-1)*p_t,:)*1000,z_yoko(k+(j-1)*(p_y-1)*p_t,:)*1000,'k');
                hold on
            end
            for k=1+(p_y-1)*(p_t-1):(p_y-1)*p_t
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
        
        
        
        %%ècê¸
        %êLÇŒÇ∑ê¸Ç™Ç»Ç¢node
        %         for k=1:p_h;
        %             for j=1:p_y:1+(p_t-2)*p_y;
        %                 x_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),1)=datai_x(i,se(j+(k-1)*p_y*(p_t-1)+ul(1,1),2));
        %                 x_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),2)=datai_x(i,se(j+(k-1)*p_y*(p_t-1)+ul(1,1),3));
        %                 y_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),1)=datai_y(i,se(j+(k-1)*p_y*(p_t-1)+ul(1,1),2));
        %                 y_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),2)=datai_y(i,se(j+(k-1)*p_y*(p_t-1)+ul(1,1),3));
        %                 z_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),1)=datai_z(i,se(j+(k-1)*p_y*(p_t-1)+ul(1,1),2));
        %                 z_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),2)=datai_z(i,se(j+(k-1)*p_y*(p_t-1)+ul(1,1),3));
        %                 plot3(x_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),:)*1000,y_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),:)*1000,z_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),:)*1000,'k');
        %                 hold on
        %             end
        %             for j=p_y:p_y:p_y+(p_t-2)*p_y;
        %                 x_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),1)=datai_x(i,se(j+(k-1)*p_y*(p_t-1)+ul(1,1),2));
        %                 x_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),2)=datai_x(i,se(j+(k-1)*p_y*(p_t-1)+ul(1,1),3));
        %                 y_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),1)=datai_y(i,se(j+(k-1)*p_y*(p_t-1)+ul(1,1),2));
        %                 y_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),2)=datai_y(i,se(j+(k-1)*p_y*(p_t-1)+ul(1,1),3));
        %                 z_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),1)=datai_z(i,se(j+(k-1)*p_y*(p_t-1)+ul(1,1),2));
        %                 z_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),2)=datai_z(i,se(j+(k-1)*p_y*(p_t-1)+ul(1,1),3));
        %                 plot3(x_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),:)*1000,y_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),:)*1000,z_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),:)*1000,'k');
        %                 hold on
        %             end
        %         end
        for k=1:p_h
            for j=1:p_y:1+(p_t-2)*p_y
                x_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),1)=datai_x(i,se(j+(k-1)*p_y*(p_t-1)+ul(1,1),2));
                x_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),2)=datai_x(i,se(j+(k-1)*p_y*(p_t-1)+ul(1,1),3));
                y_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),1)=datai_y(i,se(j+(k-1)*p_y*(p_t-1)+ul(1,1),2));
                y_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),2)=datai_y(i,se(j+(k-1)*p_y*(p_t-1)+ul(1,1),3));
                z_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),1)=datai_z(i,se(j+(k-1)*p_y*(p_t-1)+ul(1,1),2));
                z_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),2)=datai_z(i,se(j+(k-1)*p_y*(p_t-1)+ul(1,1),3));
                plot3(x_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),:)*1000,y_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),:)*1000,z_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),:)*1000,'k');
                hold on
            end
            for j=p_y:p_y:p_y+(p_t-2)*p_y
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
        
        
        %%çÇÇ≥ê¸
        for k=0:p_y*p_t:(p_h-2)*p_y*p_t
            for j=1:p_y
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
            for j=p_y*(p_t-1)+1:p_y*p_t
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
            for j=p_y+1:p_y:1+p_y*(p_t-2)
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
            for j=p_y:p_y:p_y*(p_t-1)
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
        axis equal
        %         %         xlim([-0.1*(10^8),0.25*(10^8)]);
        %            xlim([-150000,350000]);
        %          xlim([-0.2*(10^8),0.25*(10^8)]);
        %         % %         ylim([14.15*(10^8),14.50*(10^8)]);
        %            ylim([-150000,350000]);
        %          ylim([-0.2*(10^8),0.25*(10^8)]);
        %         % %         zlim([1.20*(10^8),1.55*(10^8)]);
        %            zlim([-200000,1000000]);
        %                    zlim([70, 350]);
        Frame(l) = getframe(1);
        hold off
    end
end

v = VideoWriter('speedup');
videoName = ['video\',muscle_name, '_MA=1.0_new_èáìÆóÕäw'];
v = VideoWriter(videoName);
v.Quality=50;

open(v);
writeVideo(v,Frame);
close(v);
