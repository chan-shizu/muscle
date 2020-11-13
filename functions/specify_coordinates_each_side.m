if i==1
    datai_x(i,1:pointNum)=data0(:,1).';  %i時における各点の座標
    datai_y(i,1:pointNum)=data0(:,2).';
    datai_z(i,1:pointNum)=data0(:,3).';
    
else
    %固定面
    datai_x(i,1:y*t)=data0(1:y*t,1);
    datai_y(i,1:y*t)=data0(1:y*t,2);
    datai_z(i,1:y*t)=data0(1:y*t,3);
    
    %強制変位を与える面
    for k=1:y*t
        % 上面の座標を指定
        datai_x(i,(pointNum(1)-k+1))=data(i,3*(y*t-k)+1);  %x座標
        datai_y(i,(pointNum(1)-k+1))=data(i,3*(y*t-k)+2);
        datai_z(i,(pointNum(1)-k+1))=data(i,3*(y*t-k)+3);
        
        % 下面の座標を指定
        datai_x(i,k)=dataBottom(i,3*(k-1)+1);  %x座標
        datai_y(i,k)=dataBottom(i,3*(k-1)+2);
        datai_z(i,k)=dataBottom(i,3*(k-1)+3);
        %             vn_x(i,(pointNum(1)-k+1))=vn_x(i-1,(pointNum(1)-k+1))+Fn_x(i-1,(pointNum(1)-k+1))*dt/2;
        %             vn_y(i,(pointNum(1)-k+1))=vn_y(i-1,(pointNum(1)-k+1))+Fn_y(i-1,(pointNum(1)-k+1))*dt/2;
        %             vn_z(i,(pointNum(1)-k+1))=vn_z(i-1,(pointNum(1)-k+1))+Fn_z(i-1,(pointNum(1)-k+1))*dt/2;
        %             datai_x(i,(pointNum(1)-k+1))=datai_x(i-1,(pointNum(1)-k+1))+dt*vn_x(i,(pointNum(1)-k+1));
        %             datai_y(i,(pointNum(1)-k+1))=datai_y(i-1,(pointNum(1)-k+1))+dt*vn_y(i,(pointNum(1)-k+1));
        %             datai_z(i,(pointNum(1)-k+1))=datai_z(i-1,(pointNum(1)-k+1))+dt*vn_z(i,(pointNum(1)-k+1));
    end
end