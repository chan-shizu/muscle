y=5;
t=5;
h=6;

for i=1:y*t*h
    datai_x(1,i) = i^2; 
end

i = 1;
for j=1:h-2
    for  n =1:y-2
        distancex(1+n+j*y*t,1) = abs(datai_x(i,1+n+j*y*t) - datai_x(i,1+n+(j-1)*y*t));
        distancex(1+n+j*y*t,2) = abs(datai_x(i,1+n+j*y*t) - datai_x(i,1+n+(j+1)*y*t)) ;
        distancex(1+n+j*y*t,3) = abs(datai_x(i,1+n+j*y*t) - datai_x(i,n+j*y*t));
        distancex(1+n+j*y*t,4) = abs(datai_x(i,1+n+j*y*t) - datai_x(i,2+n+j*y*t));
        distancex(1+n+j*y*t,5) = sum(distancex(1+n+j*y*t,1:4));
%         fivePointAverageDatay(1+n+j*y*t) = (datai_y(i,1+n+j*y*t) +  datai_y(i,1+n+(j-1)*y*t) +  datai_y(i,1+n+(j+1)*y*t) ...
%             +  datai_y(i,n+j*y*t)  +  datai_y(i,2+n+j*y*t))/5;
%         fivePointAverageDataz(1+n+j*y*t) = (datai_z(i,1+n+j*y*t) +  datai_z(i,1+n+(j-1)*y*t) +  datai_z(i,1+n+(j+1)*y*t) ...
%             +  datai_z(i,n+j*y*t)  +  datai_z(i,2+n+j*y*t))/5;
    end
end

for j=1:h-2
    for  n =1:y-2
        distancex(1+n+j*y*t+(t-1)*y,1) = abs(datai_x(i,1+n+j*y*t+(t-1)*y) - datai_x(i,1+n+(j-1)*y*t)+(t-1)*y);
        distancex(1+n+j*y*t+(t-1)*y,2) = abs(datai_x(i,1+n+j*y*t+(t-1)*y) - datai_x(i,1+n+(j+1)*y*t)+(t-1)*y) ;
        distancex(1+n+j*y*t+(t-1)*y,3) = abs(datai_x(i,1+n+j*y*t+(t-1)*y) - datai_x(i,n+j*y*t+(t-1)*y));
        distancex(1+n+j*y*t+(t-1)*y,4) = abs(datai_x(i,1+n+j*y*t+(t-1)*y) - datai_x(i,2+n+j*y*t+(t-1)*y));
%         fivePointAverageDatay(1+n+(t-1)*y+j*y*t) = (datai_y(i,1+n+(t-1)*y+j*y*t) +  datai_y(i,1+n+(t-1)*y+(j-1)*y*t) +  datai_y(i,1+n+(t-1)*y+(j+1)*y*t) ...
%             +  datai_y(i,n+(t-1)*y+j*y*t)  +  datai_y(i,2+n+(t-1)*y+j*y*t))/5;
%         fivePointAverageDataz(1+n+(t-1)*y+j*y*t) = (datai_z(i,1+n+(t-1)*y+j*y*t) +  datai_z(i,1+n+(t-1)*y+(j-1)*y*t) +  datai_z(i,1+n+(t-1)*y+(j+1)*y*t) ...
%             +  datai_z(i,n+(t-1)*y+j*y*t)  +  datai_z(i,2+n+(t-1)*y+j*y*t))/5;
    end
end

for n=1:h-2
    for j =1:t-2
        distancex(1+j*y+n*y*t,1) = abs(datai_x(i,1+j*y+n*y*t) - datai_x(i,1+j*y+(n-1)*y*t));
        distancex(1+j*y+n*y*t,2) = abs(datai_x(i,1+j*y+n*y*t) - datai_x(i,1+j*y+(n+1)*y*t)) ;
        distancex(1+j*y+n*y*t,3) = abs(datai_x(i,1+j*y+n*y*t) - datai_x(i,1+(j-1)*y+n*y*t));
        distancex(1+j*y+n*y*t,4) = abs(datai_x(i,1+j*y+n*y*t) - datai_x(i,1+(j+1)*y+n*y*t));
%         fivePointAverageDatax(1+j*y+n*y*t) = (datai_x(i,1+j*y+n*y*t) +  datai_x(i,1+j*y+(n-1)*y*t) +  datai_x(i,1+j*y+(n+1)*y*t) ...
%             +  datai_x(i,1+(j-1)*y+n*y*t)  +  datai_x(i,1+(j+1)*y+n*y*t))/5;
%         fivePointAverageDatay(1+j*y+n*y*t) = (datai_y(i,1+j*y+n*y*t) +  datai_y(i,1+j*y+(n-1)*y*t) +  datai_y(i,1+j*y+(n+1)*y*t) ...
%             +  datai_y(i,1+(j-1)*y+n*y*t)  +  datai_y(i,1+(j+1)*y+n*y*t))/5;
%         fivePointAverageDataz(1+j*y+n*y*t) = (datai_z(i,1+j*y+n*y*t) +  datai_z(i,1+j*y+(n-1)*y*t) +  datai_z(i,1+j*y+(n+1)*y*t) ...
%             +  datai_z(i,1+(j-1)*y+n*y*t)  +  datai_z(i,1+(j+1)*y+n*y*t))/5;
    end
end

for n=1:h-2
    for j =1:t-2
        distancex(y+j*y+n*y*t,1) = abs(datai_x(i,y+j*y+n*y*t) - datai_x(i,y+j*y+(n-1)*y*t));
        distancex(y+j*y+n*y*t,2) = abs(datai_x(i,y+j*y+n*y*t) - datai_x(i,y+j*y+(n+1)*y*t)) ;
        distancex(y+j*y+n*y*t,3) = abs(datai_x(i,y+j*y+n*y*t) - datai_x(i,y+(j-1)*y+n*y*t));
        distancex(y+j*y+n*y*t,4) = abs(datai_x(i,y+j*y+n*y*t) - datai_x(i,y+(j+1)*y+n*y*t));
%         fivePointAverageDatax(y+j*y+n*y*t) = (datai_x(i,y+j*y+n*y*t) +  datai_x(i,y+j*y+(n-1)*y*t) +  datai_x(i,y+j*y+(n+1)*y*t) ...
%             +  datai_x(i,y+(j-1)*y+n*y*t)  +  datai_x(i,y+(j+1)*y+n*y*t))/5;
%         fivePointAverageDatay(y+j*y+n*y*t) = (datai_y(i,y+j*y+n*y*t) +  datai_y(i,y+j*y+(n-1)*y*t) +  datai_y(i,y+j*y+(n+1)*y*t) ...
%             +  datai_y(i,y+(j-1)*y+n*y*t)  +  datai_y(i,y+(j+1)*y+n*y*t))/5;
%         fivePointAverageDataz(y+j*y+n*y*t) = (datai_z(i,y+j*y+n*y*t) +  datai_z(i,y+j*y+(n-1)*y*t) +  datai_z(i,y+j*y+(n+1)*y*t) ...
%             +  datai_z(i,y+(j-1)*y+n*y*t)  +  datai_z(i,y+(j+1)*y+n*y*t))/5;
    end
end