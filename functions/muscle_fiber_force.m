Lcesh(1,s-ul(1,2)) = 0.5 * length_s(1,s);
f_Lce(i,s-ul(1,2)) = exp(-( lengthen_s(i,s)/Lcesh(1,s-ul(1,2))).^2);%占部さんの修論の式(3.31)
if i == 1
    Vce(i,s-ul(1,2)) = 0;
else
    %                 Vce(i,s-ul(1,2)) = -(length_s(i,s) - length_s(i-1,s))/0.1;%0.1ではなくdtでは?
    Vce(i,s-ul(1,2)) = -(length_s(i,s) - length_s(i-1,s))/dt;
end
Vvm = 6 * length_s(1,s);%V_vmは最大等尺性収縮中の最大速度でV_vm=6?l_ce
Vmax(1,s-ul(1,2)) = Vvm .* ( 1 - Ver .* (1 - ActivityLevel(i) * f_Lce(i,s-ul(1,2))));

%占部さんの修論式(3-32)
if Vce(i,s-ul(1,2)) <= -Vmax(1,s-ul(1,2))
    f_Vce(i,s-ul(1,2)) = 0;
elseif (Vce(i,s-ul(1,2)) > -Vmax(1,s-ul(1,2)))  &&  (Vce(i,s-ul(1,2)) < 0)
    f_Vce(i,s-ul(1,2)) = (Vsh * Vmax(1,s-ul(1,2)) + Vsh * Vce(i,s-ul(1,2))) / (Vsh *Vmax(1,s-ul(1,2)) - Vce(i,s-ul(1,2)));
elseif Vce(i,s-ul(1,2)) >= 0
    f_Vce(i,s-ul(1,2)) = (Vsh * Vshl * Vmax(1,s-ul(1,2)) + Vml * Vce(i,s-ul(1,2))) / (Vsh * Vshl * Vmax(1,s-ul(1,2)) + Vce(i,s-ul(1,2)));%vmax抜けてない?
    %             else
    %                 f_Vce(i,s-ul(1,2)) = 0
end

length_per(i,s-ul(1,2)) =lengthen_s(i,s)/ length_s(1,s);
if lengthen_s(i,s)<=0
    HillPassive(i,s-ul(1,2)) =0;
else
    HillPassive(i,s-ul(1,2)) = Fmax/exp(PEsh)* (exp(lengthen_s(i,s)/ length_s(1,s)* PEsh/PExm) -1); %
end
HillActive(i,s-ul(1,2)) = ActivityLevel(i) .* f_Lce(i,s-ul(1,2)) .* f_Vce(i,s-ul(1,2)) .* Fmax'; %占部さんの修論，式(3-29)
if i==1
    springF(i,s)=0;
else
    springF(i,s)=0;%HillActive(i,s-ul(1,2))+HillPassive(i,s-ul(1,2));%受動的な力と能動的な力の和
    muscleFiberF(i,s)=HillActive(i,s-ul(1,2))+HillPassive(i,s-ul(1,2));%受動的な力と能動的な力の和
    %                         springF(i,s)=-(HillActive(i,s-ul(1,2))+HillPassive(i,s-ul(1,2)));%*scale;%受動的な力と能動的な力の和
end