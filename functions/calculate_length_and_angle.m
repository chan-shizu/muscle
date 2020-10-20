for s=1:seNum(1)
    length_s(i,s)=sqrt((datai_x(i,se(s,3))-datai_x(i,se(s,2)))^2+(datai_y(i,se(s,3))-datai_y(i,se(s,2)))^2+(datai_z(i,se(s,3))-datai_z(i,se(s,2)))^2);
    lengthen_s(i,s)=length_s(i,s)-length_s(1,s);%自然長からの長さ
    
    %角度計算
    if datai_y(i,se(s,3))-datai_y(i,se(s,2))==0&&datai_x(i,se(s,3))-datai_x(i,se(s,2))>=0
        se(s,5)=0;
    elseif datai_y(i,se(s,3))-datai_y(i,se(s,2))==0&&datai_x(i,se(s,3))-datai_x(i,se(s,2))<=0
        se(s,5)=pi;
    elseif datai_x(i,se(s,3))-datai_x(i,se(s,2))==0&&datai_y(i,se(s,3))-datai_y(i,se(s,2))>=0
        se(s,5)=pi/2;
    elseif datai_x(i,se(s,3))-datai_x(i,se(s,2))==0&&datai_y(i,se(s,3))-datai_y(i,se(s,2))<=0
        se(s,5)=-pi/2;
    elseif datai_x(i,se(s,3))-datai_x(i,se(s,2))<=0&&datai_y(i,se(s,3))-datai_y(i,se(s,2))>=0
        se(s,5)=atan((datai_y(i,se(s,3))-datai_y(i,se(s,2)))/(datai_x(i,se(s,3))-datai_x(i,se(s,2))))+pi;   %角度β(xy平面)
    elseif datai_x(i,se(s,3))-datai_x(i,se(s,2))<=0&&datai_y(i,se(s,3))-datai_y(i,se(s,2))<=0
        se(s,5)=atan((datai_y(i,se(s,3))-datai_y(i,se(s,2)))/(datai_x(i,se(s,3))-datai_x(i,se(s,2))))+pi;
    else
        se(s,5)=atan((datai_y(i,se(s,3))-datai_y(i,se(s,2)))/(datai_x(i,se(s,3))-datai_x(i,se(s,2))));
    end
    
    se(s,4)=atan((datai_z(i,se(s,3))-datai_z(i,se(s,2)))/sqrt((datai_x(i,se(s,3))-datai_x(i,se(s,2)))^2+(datai_y(i,se(s,3))-datai_y(i,se(s,2)))^2));   %角度α(xz平面
end