%fvAを求めるのに必要なﾍﾞｸﾄﾙ
Vec1(:,1+4*(k-1)+(j-1)*tetraNum(2)*4,i)=[datai_x(i,tetra(3+5*(j-1),k));datai_y(i,tetra(3+5*(j-1),k));datai_z(i,tetra(3+5*(j-1),k))]-[datai_x(i,tetra(5+5*(j-1),k));datai_y(i,tetra(5+5*(j-1),k));datai_z(i,tetra(5+5*(j-1),k))];   %B-D
Vec2(:,1+4*(k-1)+(j-1)*tetraNum(2)*4,i)=[datai_x(i,tetra(4+5*(j-1),k));datai_y(i,tetra(4+5*(j-1),k));datai_z(i,tetra(4+5*(j-1),k))]-[datai_x(i,tetra(5+5*(j-1),k));datai_y(i,tetra(5+5*(j-1),k));datai_z(i,tetra(5+5*(j-1),k))];   %C-D

%fvB
Vec1(:,2+4*(k-1)+(j-1)*tetraNum(2)*4,i)=[datai_x(i,tetra(4+5*(j-1),k));datai_y(i,tetra(4+5*(j-1),k));datai_z(i,tetra(4+5*(j-1),k))]-[datai_x(i,tetra(2+5*(j-1),k));datai_y(i,tetra(2+5*(j-1),k));datai_z(i,tetra(2+5*(j-1),k))];   %C-A
Vec2(:,2+4*(k-1)+(j-1)*tetraNum(2)*4,i)=[datai_x(i,tetra(5+5*(j-1),k));datai_y(i,tetra(5+5*(j-1),k));datai_z(i,tetra(5+5*(j-1),k))]-[datai_x(i,tetra(2+5*(j-1),k));datai_y(i,tetra(2+5*(j-1),k));datai_z(i,tetra(2+5*(j-1),k))];   %D-A

%fvC
Vec1(:,3+4*(k-1)+(j-1)*tetraNum(2)*4,i)=[datai_x(i,tetra(5+5*(j-1),k));datai_y(i,tetra(5+5*(j-1),k));datai_z(i,tetra(5+5*(j-1),k))]-[datai_x(i,tetra(3+5*(j-1),k));datai_y(i,tetra(3+5*(j-1),k));datai_z(i,tetra(3+5*(j-1),k))];   %D-B
Vec2(:,3+4*(k-1)+(j-1)*tetraNum(2)*4,i)=[datai_x(i,tetra(2+5*(j-1),k));datai_y(i,tetra(2+5*(j-1),k));datai_z(i,tetra(2+5*(j-1),k))]-[datai_x(i,tetra(3+5*(j-1),k));datai_y(i,tetra(3+5*(j-1),k));datai_z(i,tetra(3+5*(j-1),k))];   %A-B

%fvD
Vec1(:,4+4*(k-1)+(j-1)*tetraNum(2)*4,i)=[datai_x(i,tetra(2+5*(j-1),k));datai_y(i,tetra(2+5*(j-1),k));datai_z(i,tetra(2+5*(j-1),k))]-[datai_x(i,tetra(4+5*(j-1),k));datai_y(i,tetra(4+5*(j-1),k));datai_z(i,tetra(4+5*(j-1),k))];   %A-C
Vec2(:,4+4*(k-1)+(j-1)*tetraNum(2)*4,i)=[datai_x(i,tetra(3+5*(j-1),k));datai_y(i,tetra(3+5*(j-1),k));datai_z(i,tetra(3+5*(j-1),k))]-[datai_x(i,tetra(4+5*(j-1),k));datai_y(i,tetra(4+5*(j-1),k));datai_z(i,tetra(4+5*(j-1),k))];   %B-C
