if i==1
    vn_x(i,1:pointNum(1))=0;   %1列目:node4 2列目:node5
    vn_y(i,1:pointNum(1))=0;
    vn_z(i,1:pointNum(1))=0;
    
else
    %固定面と強制変位面以外を平等に変位させる(速度ベルレ法)
    vn_x(i,(y*t+1):(pointNum(1)-y*t))=vn_x(i-1,(y*t+1):(pointNum(1)-y*t))+Fn_x(i-1,(y*t+1):(pointNum(1)-y*t))*dt/2;
    vn_y(i,(y*t+1):(pointNum(1)-y*t))=vn_y(i-1,(y*t+1):(pointNum(1)-y*t))+Fn_y(i-1,(y*t+1):(pointNum(1)-y*t))*dt/2;
    vn_z(i,(y*t+1):(pointNum(1)-y*t))=vn_z(i-1,(y*t+1):(pointNum(1)-y*t))+Fn_z(i-1,(y*t+1):(pointNum(1)-y*t))*dt/2;
    
    datai_x(i,(y*t+1):(pointNum(1)-y*t))=datai_x(i-1,(y*t+1):(pointNum(1)-y*t))+dt*vn_x(i,(y*t+1):(pointNum(1)-y*t));
    datai_y(i,(y*t+1):(pointNum(1)-y*t))=datai_y(i-1,(y*t+1):(pointNum(1)-y*t))+dt*vn_y(i,(y*t+1):(pointNum(1)-y*t));
    %         datai_z(i,(y*t+1):(pointNum(1)-y*t))=datai_z(i-1,(y*t+1):(pointNum(1)-y*t))+dt*vn_z(i,(y*t+1):(pointNum(1)-y*t));
    for n=1:(h-1)
        if data0(1,3)<data0(y*t+1,3)
            datai_z(i,n*y*t+1:(n+1)*y*t) = dataBottom(i,3) + n/(h-1)*(dataBottom(i,3)-data(i,3));
        else
            datai_z(i,n*y*t+1:(n+1)*y*t) = dataBottom(i,3) - n/(h-1)*(dataBottom(i,3)-data(i,3));
        end
    end
end