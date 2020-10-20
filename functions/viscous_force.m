if i==1
    vis_x(i,1:pointNum(1))=0;
    vis_y(i,1:pointNum(1))=0;
    vis_z(i,1:pointNum(1))=0;
else
    vis_x(i,1:pointNum(1))=-1*conc*(datai_x(i,1:pointNum(1))-datai_x(i-1,1:pointNum(1)))/dt;
    vis_y(i,1:pointNum(1))=-1*conc*(datai_y(i,1:pointNum(1))-datai_y(i-1,1:pointNum(1)))/dt;
    vis_z(i,1:pointNum(1))=-1*conc*(datai_z(i,1:pointNum(1))-datai_z(i-1,1:pointNum(1)))/dt;
end