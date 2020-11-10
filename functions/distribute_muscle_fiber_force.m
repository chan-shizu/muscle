for k=ul(2)+1:ul(3)
    %                 Fs_x(i,se(k,2),7) = muscleFiberF(i,se(k,2));
    %                 muscleFiberF_x(i,se(k,2)) = muscleFiberF(i,k)*cos(se(k,4))*cos(se(k,5));
    %                 muscleFiberF_y(i,se(k,2)) = muscleFiberF(i,k)*cos(se(k,4))*sin(se(k,5));
    %                 muscleFiberF_z(i,se(k,2)) = muscleFiberF(i,k)*sin(se(k,4));
    muscleFiberF_x(i,se(k,3)) = -muscleFiberF(i,k)*cos(se(k,4))*cos(se(k,5));
    muscleFiberF_y(i,se(k,3)) = -muscleFiberF(i,k)*cos(se(k,4))*sin(se(k,5));
    muscleFiberF_z(i,se(k,3)) = -muscleFiberF(i,k)*sin(se(k,4));
end
% for k=ul(2)+1:ul(2)+t*y*(h/2-1)
%     muscleFiberF_x(i,se(k,3)) = muscleFiberF(i,k)*cos(se(k,4))*cos(se(k,5));
%     muscleFiberF_y(i,se(k,3)) = muscleFiberF(i,k)*cos(se(k,4))*sin(se(k,5));
%     muscleFiberF_z(i,se(k,3)) = muscleFiberF(i,k)*sin(se(k,4));
% end
% for k=ul(2)+t*y*(h/2-1)+1:ul(3)
%     muscleFiberF_x(i,se(k,3)) = -muscleFiberF(i,k)*cos(se(k,4))*cos(se(k,5));
%     muscleFiberF_y(i,se(k,3)) = -muscleFiberF(i,k)*cos(se(k,4))*sin(se(k,5));
%     muscleFiberF_z(i,se(k,3)) = -muscleFiberF(i,k)*sin(se(k,4));
% end