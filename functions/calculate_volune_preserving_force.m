coefficient = 1;
% if data(i,3)-data(1,3)>0    %こっちが標準,それぞれの四面体ごとに体積保存係数を設定
    fvA(:,k+(j-1)*tetraNum(2))=1/6*springkVPk(2,k+(j-1)*tetraNum(2))*coefficient*(V(k+(j-1)*tetraNum(2),1)-V0(k+(j-1)*tetraNum(2),1))/V0(k+(j-1)*tetraNum(2),1)^2*cross(Vec2(:,1+4*(k-1)+(j-1)*tetraNum(2)*4,i),Vec1(:,1+4*(k-1)+(j-1)*tetraNum(2)*4,i));
    fvB(:,k+(j-1)*tetraNum(2))=1/6*springkVPk(2,k+(j-1)*tetraNum(2))*coefficient*(V(k+(j-1)*tetraNum(2),1)-V0(k+(j-1)*tetraNum(2),1))/V0(k+(j-1)*tetraNum(2),1)^2*cross(Vec1(:,2+4*(k-1)+(j-1)*tetraNum(2)*4,i),Vec2(:,2+4*(k-1)+(j-1)*tetraNum(2)*4,i));
    fvC(:,k+(j-1)*tetraNum(2))=1/6*springkVPk(2,k+(j-1)*tetraNum(2))*coefficient*(V(k+(j-1)*tetraNum(2),1)-V0(k+(j-1)*tetraNum(2),1))/V0(k+(j-1)*tetraNum(2),1)^2*cross(Vec2(:,3+4*(k-1)+(j-1)*tetraNum(2)*4,i),Vec1(:,3+4*(k-1)+(j-1)*tetraNum(2)*4,i));
    fvD(:,k+(j-1)*tetraNum(2))=1/6*springkVPk(2,k+(j-1)*tetraNum(2))*coefficient*(V(k+(j-1)*tetraNum(2),1)-V0(k+(j-1)*tetraNum(2),1))/V0(k+(j-1)*tetraNum(2),1)^2*cross(Vec1(:,4+4*(k-1)+(j-1)*tetraNum(2)*4,i),Vec2(:,4+4*(k-1)+(j-1)*tetraNum(2)*4,i));
% else
%     fvA(:,k+(j-1)*tetraNum(2))=1/6*springkVPk(2,k+(j-1)*tetraNum(2))*coefficient*(V(k+(j-1)*tetraNum(2),1)-V0(k+(j-1)*tetraNum(2),1))/V0(k+(j-1)*tetraNum(2),1)^2*cross(Vec1(:,1+4*(k-1)+(j-1)*tetraNum(2)*4,i),Vec2(:,1+4*(k-1)+(j-1)*tetraNum(2)*4,i));
%     fvB(:,k+(j-1)*tetraNum(2))=1/6*springkVPk(2,k+(j-1)*tetraNum(2))*coefficient*(V(k+(j-1)*tetraNum(2),1)-V0(k+(j-1)*tetraNum(2),1))/V0(k+(j-1)*tetraNum(2),1)^2*cross(Vec2(:,2+4*(k-1)+(j-1)*tetraNum(2)*4,i),Vec1(:,2+4*(k-1)+(j-1)*tetraNum(2)*4,i));
%     fvC(:,k+(j-1)*tetraNum(2))=1/6*springkVPk(2,k+(j-1)*tetraNum(2))*coefficient*(V(k+(j-1)*tetraNum(2),1)-V0(k+(j-1)*tetraNum(2),1))/V0(k+(j-1)*tetraNum(2),1)^2*cross(Vec1(:,3+4*(k-1)+(j-1)*tetraNum(2)*4,i),Vec2(:,3+4*(k-1)+(j-1)*tetraNum(2)*4,i));
%     fvD(:,k+(j-1)*tetraNum(2))=1/6*springkVPk(2,k+(j-1)*tetraNum(2))*coefficient*(V(k+(j-1)*tetraNum(2),1)-V0(k+(j-1)*tetraNum(2),1))/V0(k+(j-1)*tetraNum(2),1)^2*cross(Vec2(:,4+4*(k-1)+(j-1)*tetraNum(2)*4,i),Vec1(:,4+4*(k-1)+(j-1)*tetraNum(2)*4,i));
% end