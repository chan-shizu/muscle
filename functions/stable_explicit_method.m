%     vn_x(i,(y*t+1):(pointNum(1)-y*t))=vn_x(i-1,(y*t+1):(pointNum(1)-y*t))+Fn_x(i-1,(y*t+1):(pointNum(1)-y*t))*dt/2;
%     vn_y(i,(y*t+1):(pointNum(1)-y*t))=vn_y(i-1,(y*t+1):(pointNum(1)-y*t))+Fn_y(i-1,(y*t+1):(pointNum(1)-y*t))*dt/2;
%     vn_z(i,(y*t+1):(pointNum(1)-y*t))=vn_z(i-1,(y*t+1):(pointNum(1)-y*t))+Fn_z(i-1,(y*t+1):(pointNum(1)-y*t))*dt/2;
    
    beta1 = ((1+dt*c/(mass*2)+springkVPk(1,1)*(dt^2)/(mass*4))^-1)*(1+dt*c/(mass*2))
    beta2 = ((1+dt*c/(mass*2)+springkVPk(1,1)*(dt^2)/(mass*4))^-1)/2
%     K1z = Fn_z
%     
%     k2x = 
%     
%     datai_x(i,(y*t+1):(pointNum(1)-y*t))=datai_x(i-1,(y*t+1):(pointNum(1)-y*t))+dt*vn_x(i,(y*t+1):(pointNum(1)-y*t));
%     datai_y(i,(y*t+1):(pointNum(1)-y*t))=datai_y(i-1,(y*t+1):(pointNum(1)-y*t))+dt*vn_y(i,(y*t+1):(pointNum(1)-y*t));
%     datai_z(i,(y*t+1):(pointNum(1)-y*t))=datai_z(i-1,(y*t+1):(pointNum(1)-y*t))+dt*vn_z(i,(y*t+1):(pointNum(1)-y*t));