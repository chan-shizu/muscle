Lcesh(1,muscleFiberNumber) = 0.5 * muscleFiberLength(1,muscleFiberNumber);
f_Lce(i,muscleFiberNumber) = exp(-( muscleFiberLengthen(i,muscleFiberNumber)/Lcesh(1,muscleFiberNumber)).^2);%占部さんの修論の式(3.31)
if i == 1
    Vce(i,muscleFiberNumber) = 0;
else
    %                 Vce(i,muscleFiberNumber = -(length_s(i,s) - length_s(i-1,s))/0.1;%0.1ではなくdtでは?
    Vce(i,muscleFiberNumber) = -(muscleFiberLength(i,muscleFiberNumber) - muscleFiberLength(i-1,muscleFiberNumber))/dt;
end
Vvm = 6 * muscleFiberLength(1,muscleFiberNumber);%V_vmは最大等尺性収縮中の最大速度でV_vm=6?l_ce
Vmax(1,muscleFiberNumber) = Vvm .* ( 1 - Ver .* (1 - ActivityLevel(i) * f_Lce(i,muscleFiberNumber)));

%占部さんの修論式(3-32)
if Vce(i,muscleFiberNumber) <= -Vmax(1,muscleFiberNumber)
    f_Vce(i,muscleFiberNumber) = 0;
elseif (Vce(i,muscleFiberNumber) > -Vmax(1,muscleFiberNumber)) &&  (Vce(i,muscleFiberNumber) < 0)
    f_Vce(i,muscleFiberNumber) = (Vsh * Vmax(1,muscleFiberNumber) + Vsh * Vce(i,muscleFiberNumber)) / (Vsh *Vmax(1,muscleFiberNumber) - Vce(i,muscleFiberNumber));
elseif Vce(i,muscleFiberNumber) >= 0
    f_Vce(i,muscleFiberNumber) = (Vsh * Vshl * Vmax(1,muscleFiberNumber) + Vml * Vce(i,muscleFiberNumber)) / (Vsh * Vshl * Vmax(1,muscleFiberNumber) + Vce(i,muscleFiberNumber));%vmax抜けてない?
    %             else
    %                 f_Vce(i,muscleFiberNumber = 0
end

length_per(i,muscleFiberNumber) = muscleFiberLengthen(i,muscleFiberNumber)/ muscleFiberLength(1,muscleFiberNumber);
if muscleFiberLengthen(i,muscleFiberNumber)<=0
    HillPassive(i,muscleFiberNumber) =0;
else
    HillPassive(i,muscleFiberNumber) = Fmax/exp(PEsh)* (exp(muscleFiberLengthen(i,muscleFiberNumber)/ muscleFiberLength(1,muscleFiberNumber)* PEsh/PExm) -1); %
end
HillActive(i,muscleFiberNumber) = ActivityLevel(i) * f_Lce(i,muscleFiberNumber) * f_Vce(i,muscleFiberNumber) * Fmax'; %占部さんの修論，式(3-29)
% if i==1
%     springF(i,s)=0;
% else
%     springF(i,s)=0;%HillActive(i,muscleFiberNumber+HillPassive(i,muscleFiberNumber;%受動的な力と能動的な力の和です
%     muscleFiberF(i,s)=HillActive(i,muscleFiberNumber+HillPassive(i,muscleFiberNumber;%受動的な力と能動的な力の和
%     %                         springF(i,s)=-(HillActive(i,muscleFiberNumber+HillPassive(i,muscleFiberNumber);%*scale;%受動的な力と能動的な力の和
% end

for massPointsNumber=1:h-1 %高さ方向の質点の間隔の数
    muscleFiberF(i,ul(2)+muscleFiberNumber+(massPointsNumber-1)*y*t)=(HillActive(i,muscleFiberNumber)+HillPassive(i,muscleFiberNumber));%受動的な力と能動的な力の和(i,muscleFiberNumber) = muscleFiberLength(muscleFiberNumber) + length_s(i,ul(2)+muscleFiberNumber+(massPointsNumber-1)*y*t);
end

