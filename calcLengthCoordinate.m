function length = calcLengthCoordinate(fiberCoordinate)

[nump,numq] = size(fiberCoordinate)
if numq < 3
    error('calcFiberLength: xyz data is necessary!')
elseif numq > 3
    error('calcFiberLength: too many dataset!')
else
    
  
    for p=1:nump
        if p==1
            length = 0;
        else
        pos1 = fiberCoordinate(p-1,:);
        pos2 = fiberCoordinate(p,:);
        vec = pos1 - pos2;
        tmp = norm(vec);%sqrt(vec(1)^2 +vec(2)^2 +vec(3)^2 )
        %length = length + tmp;
        length = [length; length(p-1) + tmp];%配列　水平連結
        end
        
    end
end

% 入力
% fiberCoordinate =
%     0.1360    1.3221    0.1013
%     0.1397    1.3138    0.1011
%     0.1428    1.3076    0.1001
%     0.1446    1.3036    0.0993
%     0.1470    1.2984    0.0980
%     0.1507    1.2900    0.0957
%     0.1542    1.2820    0.0934
%     0.1574    1.2735    0.0896
%     0.1614    1.2639    0.0863
%     0.1642    1.2561    0.0845
%     0.1665    1.2467    0.0817
%     0.1699    1.2368    0.0798
%     0.1729    1.2300    0.0784
%     0.1757    1.2239    0.0771
%     0.1784    1.2180    0.0758
%     0.1818    1.2113    0.0740
%
% 出力
% lengthCoordinate =
%          0
%     0.0091
%     0.0161
%     0.0206
%     0.0264
%     0.0359
%     0.0449
%     0.0548
%     0.0656
%     0.0742
%     0.0842
%     0.0948
%     0.1024
%     0.1092
%     0.1158
%     0.1236


