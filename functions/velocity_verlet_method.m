if i==1
    vn_x(i,1:pointNum(1))=0;   %1列目:node4 2列目:node5
    vn_y(i,1:pointNum(1))=0;
    vn_z(i,1:pointNum(1))=0;
    
elseif i <3
    %固定面と強制変位面以外を平等に変位させる(速度ベルレ法)
    vn_x(i,(y*t+1):(pointNum(1)-y*t))=vn_x(i-1,(y*t+1):(pointNum(1)-y*t))+Fn_x(i-1,(y*t+1):(pointNum(1)-y*t))*dt/(2*mass);
    vn_y(i,(y*t+1):(pointNum(1)-y*t))=vn_y(i-1,(y*t+1):(pointNum(1)-y*t))+Fn_y(i-1,(y*t+1):(pointNum(1)-y*t))*dt/(2*mass);
    vn_z(i,(y*t+1):(pointNum(1)-y*t))=vn_z(i-1,(y*t+1):(pointNum(1)-y*t))+Fn_z(i-1,(y*t+1):(pointNum(1)-y*t))*dt/(2*mass);
    
    datai_x(i,(y*t+1):(pointNum(1)-y*t))=datai_x(i-1,(y*t+1):(pointNum(1)-y*t))+dt*vn_x(i,(y*t+1):(pointNum(1)-y*t));
    datai_y(i,(y*t+1):(pointNum(1)-y*t))=datai_y(i-1,(y*t+1):(pointNum(1)-y*t))+dt*vn_y(i,(y*t+1):(pointNum(1)-y*t));
%         datai_z(i,(y*t+1):(pointNum(1)-y*t))=datai_z(i-1,(y*t+1):(pointNum(1)-y*t))+dt*vn_z(i,(y*t+1):(pointNum(1)-y*t));
    
    %     vn_x(i,1:(pointNum(1)-y*t))=vn_x(i-1,1:(pointNum(1)-y*t))+Fn_x(i-1,1:(pointNum(1)-y*t))*dt/(2*mass);
    %     vn_y(i,1:(pointNum(1)-y*t))=vn_y(i-1,1:(pointNum(1)-y*t))+Fn_y(i-1,1:(pointNum(1)-y*t))*dt/(2*mass);
    %     vn_z(i,1:(pointNum(1)-y*t))=vn_z(i-1,1:(pointNum(1)-y*t))+Fn_z(i-1,1:(pointNum(1)-y*t))*dt/(2*mass);
    %
    %     datai_x(i,1:(pointNum(1)-y*t))=datai_x(i-1,1:(pointNum(1)-y*t))+dt*vn_x(i,1:(pointNum(1)-y*t));
    %     datai_y(i,1:(pointNum(1)-y*t))=datai_y(i-1,1:(pointNum(1)-y*t))+dt*vn_y(i,1:(pointNum(1)-y*t));
    %     datai_z(i,1:(pointNum(1)-y*t))=datai_z(i-1,1:(pointNum(1)-y*t))+dt*vn_z(i,1:(pointNum(1)-y*t));
    for n=1:(h-1)
        if data0(1,3)<data0(y*t+1,3)
            datai_z(i,n*y*t+1:(n+1)*y*t) = data0(1,3) + n/(h-1)*(data(i,3)-data0(1,3));
        else
            datai_z(i,n*y*t+1:(n+1)*y*t) = data0(1,3) - n/(h-1)*(data0(1,3)-data(i,3));
        end
    end
    
    
elseif i==3
    vn_x(i,(y*t+1):(pointNum(1)-y*t))=vn_x(i-1,(y*t+1):(pointNum(1)-y*t))+Fn_x(i-1,(y*t+1):(pointNum(1)-y*t))*dt/(2*mass);
    vn_y(i,(y*t+1):(pointNum(1)-y*t))=vn_y(i-1,(y*t+1):(pointNum(1)-y*t))+Fn_y(i-1,(y*t+1):(pointNum(1)-y*t))*dt/(2*mass);
    vn_z(i,(y*t+1):(pointNum(1)-y*t))=0;
    
    datai_x(i,(y*t+1):(pointNum(1)-y*t))=datai_x(i-1,(y*t+1):(pointNum(1)-y*t))+dt*vn_x(i,(y*t+1):(pointNum(1)-y*t));
    datai_y(i,(y*t+1):(pointNum(1)-y*t))=datai_y(i-1,(y*t+1):(pointNum(1)-y*t))+dt*vn_y(i,(y*t+1):(pointNum(1)-y*t));
    datai_z(i,(y*t+1):(pointNum(1)-y*t))=datai_z(i-1,(y*t+1):(pointNum(1)-y*t))+dt*vn_z(i,(y*t+1):(pointNum(1)-y*t));
    
else
    vn_x(i,(y*t+1):(pointNum(1)-y*t))=vn_x(i-1,(y*t+1):(pointNum(1)-y*t))+Fn_x(i-1,(y*t+1):(pointNum(1)-y*t))*dt/(2*mass);
    vn_y(i,(y*t+1):(pointNum(1)-y*t))=vn_y(i-1,(y*t+1):(pointNum(1)-y*t))+Fn_y(i-1,(y*t+1):(pointNum(1)-y*t))*dt/(2*mass);
    vn_z(i,(y*t+1):(pointNum(1)-y*t))=vn_z(i-1,(y*t+1):(pointNum(1)-y*t))+Fn_z(i-1,(y*t+1):(pointNum(1)-y*t))*dt/(2*mass);
    
    datai_x(i,(y*t+1):(pointNum(1)-y*t))=datai_x(i-1,(y*t+1):(pointNum(1)-y*t))+dt*vn_x(i,(y*t+1):(pointNum(1)-y*t));
    datai_y(i,(y*t+1):(pointNum(1)-y*t))=datai_y(i-1,(y*t+1):(pointNum(1)-y*t))+dt*vn_y(i,(y*t+1):(pointNum(1)-y*t));
    datai_z(i,(y*t+1):(pointNum(1)-y*t))=datai_z(i-1,(y*t+1):(pointNum(1)-y*t))+dt*vn_z(i,(y*t+1):(pointNum(1)-y*t));
    
%     for layerNumber =1:h
%         datai_z(i,y*t*(layerNumber-1)+1:y*t*layerNumber)=mean(datai_z(i,y*t*(layerNumber-1)+1:y*t*layerNumber));
%     end
    
%     for layer=1:h
%         datai_x(i,y*t*(layer-1)+1:y*t*(layer-1)+5)=data0(y*t*(layer-1)+1:y*t*(layer-1)+5,1);
%         datai_y(i,y*t*(layer-1)+1:y*t*(layer-1)+5)=data0(y*t*(layer-1)+1:y*t*(layer-1)+5,2);
% %         datai_z(i,y*t*(layer-1)+1:y*t*(layer-1)+5)=data0(y*t*(layer-1)+1:y*t*(layer-1)+5,3);
%     end
end