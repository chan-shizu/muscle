volumeData(:,1) = data0(:,1)
volumeData(:,2) = data0(:,2)
volumeData(:,3) = data0(:,3)

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
volumeInitial = sumTetraV1 + sumTetraV2 + sumTetraV3 + sumTetraV4 + sumTetraV5 + sumTetraV6;

boxVolume = tetraV1 + tetraV2 + tetraV3 + tetraV4 + tetraV5 + tetraV6;

boxSize = size(boxVolume);
for boxNumber=1:boxSize(2)
    initialBoxIncludeTetraVolume(1+6*(boxNumber-1):1+6*boxNumber) = boxVolume(boxNumber);
end