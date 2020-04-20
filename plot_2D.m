se= csvread('se_arm.csv', 0, 0); 
seNum=size(se);
data= csvread('data_slide_t.csv', 0, 0); 

timeNum = size(data); 
p_y=5;
p_t=5;
p_h=6;

% fiber_force=csvread('fiber_force.csv', 0, 0); 
datai_y=csvread('output\data_y.csv', 0, 0);
datai_z=csvread('output\data_z.csv', 0, 0);

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


%%â°ê¸
%êLÇŒÇ∑ê¸Ç™Ç»Ç¢node
for j=1:p_h
    for k=1:p_y-1;
y_yoko(k+(j-1)*(p_y-1)*p_t,1)=datai_y(i,se(k+(j-1)*(p_y-1)*p_t,2));
y_yoko(k+(j-1)*(p_y-1)*p_t,2)=datai_y(i,se(k+(j-1)*(p_y-1)*p_t,3));
z_yoko(k+(j-1)*(p_y-1)*p_t,1)=datai_z(i,se(k+(j-1)*(p_y-1)*p_t,2));
z_yoko(k+(j-1)*(p_y-1)*p_t,2)=datai_z(i,se(k+(j-1)*(p_y-1)*p_t,3));
plot(y_yoko(k+(j-1)*(p_y-1)*p_t,:)*1000,z_yoko(k+(j-1)*(p_y-1)*p_t,:)*1000,'k');
hold on
    end
    for k=1+(p_y-1)*(p_t-1):(p_y-1)*p_t;
y_yoko(k+(j-1)*(p_y-1)*p_t,1)=datai_y(i,se(k+(j-1)*(p_y-1)*p_t,2));
y_yoko(k+(j-1)*(p_y-1)*p_t,2)=datai_y(i,se(k+(j-1)*(p_y-1)*p_t,3));
z_yoko(k+(j-1)*(p_y-1)*p_t,1)=datai_z(i,se(k+(j-1)*(p_y-1)*p_t,2));
z_yoko(k+(j-1)*(p_y-1)*p_t,2)=datai_z(i,se(k+(j-1)*(p_y-1)*p_t,3));
plot(y_yoko(k+(j-1)*(p_y-1)*p_t,:)*1000,z_yoko(k+(j-1)*(p_y-1)*p_t,:)*1000,'k');
hold on
    end
end



%%ècê¸
%êLÇŒÇ∑ê¸Ç™Ç»Ç¢node
for k=1:p_h;
for j=1:p_y:1+(p_t-2)*p_y;
y_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),1)=datai_y(i,se(j+(k-1)*p_y*(p_t-1)+ul(1,1),2));
y_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),2)=datai_y(i,se(j+(k-1)*p_y*(p_t-1)+ul(1,1),3));
z_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),1)=datai_z(i,se(j+(k-1)*p_y*(p_t-1)+ul(1,1),2));
z_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),2)=datai_z(i,se(j+(k-1)*p_y*(p_t-1)+ul(1,1),3));
plot(y_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),:)*1000,z_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),:)*1000,'k');
hold on
end
for j=p_y:p_y:p_y+(p_t-2)*p_y;
y_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),1)=datai_y(i,se(j+(k-1)*p_y*(p_t-1)+ul(1,1),2));
y_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),2)=datai_y(i,se(j+(k-1)*p_y*(p_t-1)+ul(1,1),3));
z_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),1)=datai_z(i,se(j+(k-1)*p_y*(p_t-1)+ul(1,1),2));
z_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),2)=datai_z(i,se(j+(k-1)*p_y*(p_t-1)+ul(1,1),3));
plot(y_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),:)*1000,z_tate(j+(k-1)*p_y*(p_t-1)+ul(1,1),:)*1000,'k');
hold on
end
end


%%çÇÇ≥ê¸
for k=0:p_y*p_t:(p_h-2)*p_y*p_t;
for j=1:p_y;
y_high(j+k+ul(1,2),1)=datai_y(i,se(j+k+ul(1,2),2));
y_high(j+k+ul(1,2),2)=datai_y(i,se(j+k+ul(1,2),3));
z_high(j+k+ul(1,2),1)=datai_z(i,se(j+k+ul(1,2),2));
z_high(j+k+ul(1,2),2)=datai_z(i,se(j+k+ul(1,2),3));
% if max(fiber_force(i,:))==0
%     color(j+k,1)=0;
% else
% color(j+k,1)=fiber_force(i,j+k)/max(fiber_force(i,:))*255/255;
% end
plot(y_high(j+k+ul(1,2),:)*1000,z_high(j+k+ul(1,2),:)*1000,'r');
hold on
end
for j=p_y*(p_t-1)+1:p_y*p_t;
y_high(j+k+ul(1,2),1)=datai_y(i,se(j+k+ul(1,2),2));
y_high(j+k+ul(1,2),2)=datai_y(i,se(j+k+ul(1,2),3));
z_high(j+k+ul(1,2),1)=datai_z(i,se(j+k+ul(1,2),2));
z_high(j+k+ul(1,2),2)=datai_z(i,se(j+k+ul(1,2),3));
% if max(fiber_force(i,:))==0
%     color(j+k,1)=0;
% else
% color(j+k,1)=fiber_force(i,j+k)/max(fiber_force(i,:))*255/255;
% end
plot(y_high(j+k+ul(1,2),:)*1000,z_high(j+k+ul(1,2),:)*1000,'r');
hold on
end
for j=p_y+1:p_y:1+p_y*(p_t-2);
y_high(j+k+ul(1,2),1)=datai_y(i,se(j+k+ul(1,2),2));
y_high(j+k+ul(1,2),2)=datai_y(i,se(j+k+ul(1,2),3));
z_high(j+k+ul(1,2),1)=datai_z(i,se(j+k+ul(1,2),2));
z_high(j+k+ul(1,2),2)=datai_z(i,se(j+k+ul(1,2),3));
% if max(fiber_force(i,:))==0
%     color(j+k,1)=0;
% else
% color(j+k,1)=fiber_force(i,j+k)/max(fiber_force(i,:))*255/255;
% end
plot(y_high(j+k+ul(1,2),:)*1000,z_high(j+k+ul(1,2),:)*1000,'r');
hold on
end
for j=p_y:p_y:p_y*(p_t-1);
y_high(j+k+ul(1,2),1)=datai_y(i,se(j+k+ul(1,2),2));
y_high(j+k+ul(1,2),2)=datai_y(i,se(j+k+ul(1,2),3));
z_high(j+k+ul(1,2),1)=datai_z(i,se(j+k+ul(1,2),2));
z_high(j+k+ul(1,2),2)=datai_z(i,se(j+k+ul(1,2),3));
% if max(fiber_force(i,:))==0
%     color(j+k,1)=0;
% else
% color(j+k,1)=fiber_force(i,j+k)/max(fiber_force(i,:))*255/255;
% end
plot(y_high(j+k+ul(1,2),:)*1000,z_high(j+k+ul(1,2),:)*1000,'r');
hold on
end
end

xlim([-192.5,7.5]);
ylim([1055,1255]);

Frame(l) = getframe(1);
hold off
end
end

v = VideoWriter('model_2_90.avi');
v.Quality=90;

open(v);
writeVideo(v,Frame);
close(v);
