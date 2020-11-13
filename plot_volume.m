volumeData(:,1) = datai_x(i,:)';
volumeData(:,2) = datai_y(i,:)';
volumeData(:,3) = datai_z(i,:)';

for n=1:h-1
    for j=1:t-1
        for k=1:y-1
            tetraV1(k+(j-1)*(y-1)+(n-1)*(y-1)*(t-1)) = abs(dot((volumeData(n*t*y+k+(j-1)*y,:)-volumeData((n-1)*t*y+k+(j-1)*y,:)),cross((volumeData((n-1)*t*y+1+k+(j-1)*y,:)-volumeData((n-1)*t*y+k+(j-1)*y,:)),(volumeData((n-1)*t*y+k+j*y,:)-volumeData((n-1)*t*y+k+(j-1)*y,:)))))/6;
            tetraV2(k+(j-1)*(y-1)+(n-1)*(y-1)*(t-1)) = abs(dot((volumeData(n*t*y+k+(j-1)*y,:)-volumeData(n*t*y+k+(j-1)*y+1,:)),cross((volumeData((n-1)*t*y+k+j*y,:)-volumeData(n*t*y+k+(j-1)*y+1,:)),(volumeData((n-1)*t*y+1+k+(j-1)*y,:)-volumeData(n*t*y+k+(j-1)*y+1,:)))))/6;
            tetraV3(k+(j-1)*(y-1)+(n-1)*(y-1)*(t-1)) = abs(dot((volumeData(n*t*y+k+(j-1)*y,:)-volumeData(n*t*y+k+(j-1)*y+1,:)),cross((volumeData((n-1)*t*y+k+j*y,:)-volumeData(n*t*y+k+(j-1)*y+1,:)),(volumeData(n*t*y+k+j*y,:)-volumeData(n*t*y+k+(j-1)*y+1,:)))))/6;
            tetraV4(k+(j-1)*(y-1)+(n-1)*(y-1)*(t-1)) = abs(dot((volumeData(n*t*y+k+j*y+1,:)-volumeData((n-1)*t*y+k+j*y+1,:)),cross((volumeData((n-1)*t*y+k+j*y,:)-volumeData((n-1)*t*y+k+j*y+1,:)),(volumeData((n-1)*t*y+1+k+(j-1)*y,:)-volumeData((n-1)*t*y+k+j*y+1,:)))))/6;
            tetraV5(k+(j-1)*(y-1)+(n-1)*(y-1)*(t-1)) = abs(dot((volumeData(n*t*y+k+j*y+1,:)-volumeData((n-1)*t*y+1+k+(j-1)*y,:)),cross((volumeData(n*t*y+k+(j-1)*y+1,:)-volumeData((n-1)*t*y+1+k+(j-1)*y,:)),(volumeData((n-1)*t*y+k+j*y,:)-volumeData((n-1)*t*y+1+k+(j-1)*y,:)))))/6;
            tetraV6(k+(j-1)*(y-1)+(n-1)*(y-1)*(t-1)) = abs(dot((volumeData(n*t*y+k+j*y+1,:)-volumeData((n-1)*t*y+k+j*y,:)),cross((volumeData(n*t*y+k+(j-1)*y+1,:)-volumeData((n-1)*t*y+k+j*y,:)),(volumeData((n-1)*t*y+k+j*y,:)-volumeData(n*t*y+k+(j-1)*y,:)))))/6;
        end
    end
end
sumTetraV1 = sum(tetraV1);
sumTetraV2 = sum(tetraV2);
sumTetraV3 = sum(tetraV3);
sumTetraV4 = sum(tetraV4);
sumTetraV5 = sum(tetraV5);
sumTetraV6 = sum(tetraV6);
volume(i) = sumTetraV1 + sumTetraV2 + sumTetraV3 + sumTetraV4 + sumTetraV5 + sumTetraV6;
volumeRatio(i) = volume(i)/volume(1);
plot([1:i],volumeRatio);
Frame(i) = getframe(1);