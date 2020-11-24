for muscleFiberNumber=1:y*t %筋線維の数
    for massPointsNumber=1:h-1 %高さ方向の質点の間隔の数
        muscleFiberLength(i,muscleFiberNumber) = muscleFiberLength(i,muscleFiberNumber) + length_s(i,ul(2)+muscleFiberNumber+(massPointsNumber-1)*y*t);
    end
end

muscleFiberLengthen(i,:) = muscleFiberLength(i,:) - muscleFiberLength(1,:);
muscleFiberLengthPer(i,:) = muscleFiberLengthen(i,:) ./ muscleFiberLength(1,:);