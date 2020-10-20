%固定面と強制変位面以外を平等に変位させる(速度ベルレ法)
vn_x(i,(y*t+1):(pointNum(1)))=vn_x(i-1,(y*t+1):(pointNum(1)))+Fn_x(i-1,(y*t+1):(pointNum(1)))*dt/2;
vn_y(i,(y*t+1):(pointNum(1)))=vn_y(i-1,(y*t+1):(pointNum(1)))+Fn_y(i-1,(y*t+1):(pointNum(1)))*dt/2;
vn_z(i,(y*t+1):(pointNum(1)))=vn_z(i-1,(y*t+1):(pointNum(1)))+Fn_z(i-1,(y*t+1):(pointNum(1)))*dt/2;

%固定面
datai_x(i,1:y*t)=data0(1:y*t,1);
datai_y(i,1:y*t)=data0(1:y*t,2);
datai_z(i,1:y*t)=data0(1:y*t,3);

datai_x(i,(y*t+1):(pointNum(1)))=datai_x(i-1,(y*t+1):(pointNum(1)))+dt*vn_x(i,(y*t+1):(pointNum(1)));
datai_y(i,(y*t+1):(pointNum(1)))=datai_y(i-1,(y*t+1):(pointNum(1)))+dt*vn_y(i,(y*t+1):(pointNum(1)));
datai_z(i,(y*t+1):(pointNum(1)))=datai_z(i-1,(y*t+1):(pointNum(1)))+dt*vn_z(i,(y*t+1):(pointNum(1)));

for k=1:h
    datai_z(i,y*t*(k-1)+1:y*t*k) =  mean(datai_z(i,y*t*(k-1)+1:y*t*k));
end

for n=1:(h-1)
    if data0(1,3)<data0(y*t+1,3)
        datai_z(i,n*y*t+1:(n+1)*y*t) = data0(1,3) + n/(h-1)*(datai_z(i,end)-data0(1,3));
    else
        datai_z(i,n*y*t+1:(n+1)*y*t) = data0(1,3) - n/(h-1)*(data0(1,3)-datai_z(i,end));
    end
end