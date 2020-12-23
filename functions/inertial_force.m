% if i==1
    Mr_x(i,1:pointNum(1))=0;
    Mr_y(i,1:pointNum(1))=0;
    Mr_z(i,1:pointNum(1))=0;
% else
%     Mr_x(i,1:pointNum(1))=mass*(datai_x(i,1:pointNum(1))-datai_x(i-1,1:pointNum(1)))/dt/dt;
%     Mr_y(i,1:pointNum(1))=mass*(datai_y(i,1:pointNum(1))-datai_y(i-1,1:pointNum(1)))/dt/dt;
%     Mr_z(i,1:pointNum(1))=mass*(datai_z(i,1:pointNum(1))-datai_z(i-1,1:pointNum(1)))/dt/dt;
% % end