Fn_x(i,:)=FsSum_x(i,:)+Fv_x(i,:)+Mr_x(i,:)+con_x(i,:)+vis_x(i,:);%Fg_x(i,:)+Fv_x(i,:)+FsSum_x(i,:)+
Fn_y(i,:)=FsSum_y(i,:)+Fv_y(i,:)+Mr_y(i,:)+con_y(i,:)+vis_y(i,:);%Fg_y(i,:)+Fv_y(i,:)+FsSum_y(i,:)
Fn_z(i,:)=FsSum_z(i,:)+Fv_z(i,:)+Mr_z(i,:)+con_z(i,:)+vis_z(i,:);%+Fg_z(i,:)+Fv_z(i,:)+

vn_x(i,1:pointNum(1))=vn_x(i,1:pointNum(1))+Fn_x(i,1:pointNum(1))/2*dt;
vn_y(i,1:pointNum(1))=vn_y(i,1:pointNum(1))+Fn_y(i,1:pointNum(1))/2*dt;
vn_z(i,1:pointNum(1))=vn_z(i,1:pointNum(1))+Fn_z(i,1:pointNum(1))/2*dt;