for k=ul(2)+1:ul(3)
    %                 Fs_x(i,se(k,2),7) = muscleFiberF(i,se(k,2));
    %                 muscleFiberF_x(i,se(k,2)) = muscleFiberF(i,k)*cos(se(k,4))*cos(se(k,5));
    %                 muscleFiberF_y(i,se(k,2)) = muscleFiberF(i,k)*cos(se(k,4))*sin(se(k,5));
    %                 muscleFiberF_z(i,se(k,2)) = muscleFiberF(i,k)*sin(se(k,4));
    muscleFiberF_x(i,se(k,3)) = -muscleFiberF(i,k)*cos(se(k,4))*cos(se(k,5));
    muscleFiberF_y(i,se(k,3)) = -muscleFiberF(i,k)*cos(se(k,4))*sin(se(k,5));
    muscleFiberF_z(i,se(k,3)) = -muscleFiberF(i,k)*sin(se(k,4));
    
%     muscleFiberF_x(i,se_step1(k,3)) = -muscleFiberF(i,k)*cos(se_step1(k,4))*cos(se_step1(k,5));
%     muscleFiberF_y(i,se_step1(k,3)) = -muscleFiberF(i,k)*cos(se_step1(k,4))*sin(se_step1(k,5));
%     muscleFiberF_z(i,se_step1(k,3)) = -muscleFiberF(i,k)*sin(se_step1(k,4));
    
%     muscleFiberF_x((muscleFiberF_x(i,:) < 0.001) & (muscleFiberF_x(i,:) > 0)) = 0.001;
%     muscleFiberF_y((muscleFiberF_y(i,:) < 0.001) & (muscleFiberF_y(i,:) > 0)) = 0.001;
%     muscleFiberF_z((muscleFiberF_z(i,:) < 0.001) & (muscleFiberF_z(i,:) > 0)) = 0.001;
    
%     muscleFiberF_x((muscleFiberF_x > 1)) = 0;
%     muscleFiberF_y((muscleFiberF_y > 1)) = 0;
%     muscleFiberF_z((muscleFiberF_z > 5)) = 0;
%     
%     muscleFiberF_x((muscleFiberF_x < -1)) = 0;
%     muscleFiberF_y((muscleFiberF_y < -1)) = 0;
%     muscleFiberF_z((muscleFiberF_z < -5)) = 0;
%     muscleFiberF_x(i,se(k,3)) = muscleFiberF(i,k)*cos(se(k,4))*cos(se(k,5));
%     muscleFiberF_y(i,se(k,3)) = muscleFiberF(i,k)*cos(se(k,4))*sin(se(k,5));
%     muscleFiberF_z(i,se(k,3)) = muscleFiberF(i,k)*sin(se(k,4));
    
%     muscleFiberF_x(i,se(k,2)) = -muscleFiberF(i,k)*cos(se(k,4))*cos(se(k,5));
%     muscleFiberF_y(i,se(k,2)) = -muscleFiberF(i,k)*cos(se(k,4))*sin(se(k,5));
%     muscleFiberF_z(i,se(k,2)) = -muscleFiberF(i,k)*sin(se(k,4));
end
