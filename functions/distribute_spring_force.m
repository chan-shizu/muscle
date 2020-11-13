Fs_x(i,1:pointNum(1),n)=0;
Fs_y(i,1:pointNum(1),n)=0;
Fs_z(i,1:pointNum(1),n)=0;
senum=(1:seNum);

if n==1
    Fsl_x(i,1:pointNum(1))=0;
    Fsl_y(i,1:pointNum(1))=0;
    Fsl_z(i,1:pointNum(1))=0;
    Fsr_x(i,1:pointNum(1))=0;
    Fsr_y(i,1:pointNum(1))=0;
    Fsr_z(i,1:pointNum(1))=0;
    
    %左だけ検索
    for k=1:ul(1,n)
        for j=1:ul(1,n)
            %             if n ==1
            %                 ul_number = 1;
            %             else
            %                 ul_number = n-1;
            %             end
            %             for k=ul_number:ul(1,n)
            %                 for j=ul_number:ul(1,n)
            match_l(j,1)=se(j,2)==se(k,2);%(false(0)とtrue(1)の配列)
        end
        match_lNum=senum(match_l);   %同じ質点を有する三角Noを抽出
        sN=size(match_lNum);
        if sN(2)>1
            matchse_l=[];
            for s=1:sN(2)   %matchnumに入っている三角Noのxy方向力をmatchf_Aに表示
                matchse_l(1,s)=springF(i,match_lNum(1,s))*cos(se(match_lNum(1,s),4))*cos(se(match_lNum(1,s),5));  %x
                matchse_l(2,s)=springF(i,match_lNum(1,s))*cos(se(match_lNum(1,s),4))*sin(se(match_lNum(1,s),5));   %y
                matchse_l(3,s)=springF(i,match_lNum(1,s))*sin(se(match_lNum(1,s),4));   %z
            end
            Fsl_x(i,se(k,2))=sum(matchse_l(1,:));
            Fsl_y(i,se(k,2))=sum(matchse_l(2,:));
            Fsl_z(i,se(k,2))=sum(matchse_l(3,:));
            
        else
            Fsl_x(i,se(k,2))=springF(i,k)*cos(se(k,4))*cos(se(k,5));
            Fsl_y(i,se(k,2))=springF(i,k)*cos(se(k,4))*sin(se(k,5));
            Fsl_z(i,se(k,2))=springF(i,k)*sin(se(k,4));
        end
        match_l(:,1)=0;
    end
    
    %右だけ検索
    for k=1:ul(1,n)
        for j=1:ul(1,n)
            %             if n ==1
            %                 ul_number = 1;
            %             else
            %                 ul_number = n-1;
            %             end
            %             for k=ul_number:ul(1,n)
            %                 for j=ul_number:ul(1,n)
            match_r(j,1)=se(j,3)==se(k,3);
        end
        match_rNum=senum(match_r);   %同じ質点を有する三角Noを抽出
        sN=size(match_rNum);
        if sN(2)>1
            matchse_r=[];
            for s=1:sN(2)   %matchnumに入っている三角Noのxy方向力をmatchf_Aに表示
                matchse_r(1,s)=-springF(i,match_rNum(1,s))*cos(se(match_rNum(1,s),4))*cos(se(match_rNum(1,s),5));  %x
                matchse_r(2,s)=-springF(i,match_rNum(1,s))*cos(se(match_rNum(1,s),4))*sin(se(match_rNum(1,s),5));   %y
                matchse_r(3,s)=-springF(i,match_rNum(1,s))*sin(se(match_rNum(1,s),4));   %z
                
            end
            Fsr_x(i,se(k,3))=sum(matchse_r(1,:));
            Fsr_y(i,se(k,3))=sum(matchse_r(2,:));
            Fsr_z(i,se(k,3))=sum(matchse_r(3,:));
            
        else
            Fsr_x(i,se(k,3))=-springF(i,k)*cos(se(k,4))*cos(se(k,5));
            Fsr_y(i,se(k,3))=-springF(i,k)*cos(se(k,4))*sin(se(k,5));
            Fsr_z(i,se(k,3))=-springF(i,k)*sin(se(k,4));
        end
        match_r(:,1)=0;
    end
    Fs_x(i,:,1)=Fsl_x(i,:)+Fsr_x(i,:);
    Fs_y(i,:,1)=Fsl_y(i,:)+Fsr_y(i,:);
    Fs_z(i,:,1)=Fsl_z(i,:)+Fsr_z(i,:);
    
else
    Fsl_x(i,1:pointNum(1))=0;
    Fsl_y(i,1:pointNum(1))=0;
    Fsl_z(i,1:pointNum(1))=0;
    Fsr_x(i,1:pointNum(1))=0;
    Fsr_y(i,1:pointNum(1))=0;
    Fsr_z(i,1:pointNum(1))=0;
    
    
    %左だけ検索
    for k=ul(1,n-1)+1:ul(1,n)
        for j=ul(1,n-1)+1:ul(1,n)
            %             if n ==1
            %                 ul_number = 1;
            %             else
            %                 ul_number = n-1;
            %             end
            %             for k=ul_number:ul(1,n)
            %                 for j=ul_number:ul(1,n)
            match_l(j,1)=se(j,2)==se(k,2);
        end
        match_lNum=senum(match_l);   %同じ質点を有する三角Noを抽出
        sN=size(match_lNum);
        if sN(2)>1
            matchse_l=[];
            for s=1:sN(2)  %matchnumに入っている三角Noのxy方向力をmatchf_Aに表示
                %                                 if se(match_lNum(1,s),1) == 2 %筋線維ではないとき
                %                                     matchse_r(1,s)=0  %x
                %                                     matchse_r(2,s)=0  %y
                %                                     matchse_r(3,s)=0  %z
                %                                 else
                matchse_l(1,s)=springF(i,match_lNum(1,s))*cos(se(match_lNum(1,s),4))*cos(se(match_lNum(1,s),5));  %x
                matchse_l(2,s)=springF(i,match_lNum(1,s))*cos(se(match_lNum(1,s),4))*sin(se(match_lNum(1,s),5));   %y
                matchse_l(3,s)=springF(i,match_lNum(1,s))*sin(se(match_lNum(1,s),4));   %z
                %                                 end
            end
            Fsl_x(i,se(k,2))=sum(matchse_l(1,:));
            Fsl_y(i,se(k,2))=sum(matchse_l(2,:));
            Fsl_z(i,se(k,2))=sum(matchse_l(3,:));
            
        else
            Fsl_x(i,se(k,2))=springF(i,k)*cos(se(k,4))*cos(se(k,5));
            Fsl_y(i,se(k,2))=springF(i,k)*cos(se(k,4))*sin(se(k,5));
            Fsl_z(i,se(k,2))=springF(i,k)*sin(se(k,4));
        end
        match_l(:,1)=0;
    end
    
    %右だけ検索
    for k=ul(1,n-1)+1:ul(1,n)
        for j=ul(1,n-1)+1:ul(1,n)
            %             if n ==1
            %                 ul_number = 1;
            %             else
            %                 ul_number = n-1;
            %             end
            %             for k=ul_number:ul(1,n)
            %                 for j=ul(1,n-1):ul(1,n)
            match_r(j,1)=se(j,3)==se(k,3);
        end
        match_rNum=senum(match_r);   %同じ質点を有する三角Noを抽出
        sN=size(match_rNum);
        if sN(2)>1
            matchse_r=[];
            for s=1:sN(2)   %matchnumに入っている三角Noのxy方向力をmatchf_Aに表示
                %                                 if se(match_rNum(1,s),1) == 2 %筋線維ではないとき
                %                                     matchse_r(1,s)=0  %x
                %                                     matchse_r(2,s)=0  %y
                %                                     matchse_r(3,s)=0  %z
                %                                 else
                matchse_r(1,s)=-springF(i,match_rNum(1,s))*cos(se(match_rNum(1,s),4))*cos(se(match_rNum(1,s),5));  %x
                matchse_r(2,s)=-springF(i,match_rNum(1,s))*cos(se(match_rNum(1,s),4))*sin(se(match_rNum(1,s),5));   %y
                matchse_r(3,s)=-springF(i,match_rNum(1,s))*sin(se(match_rNum(1,s),4));   %z
                %                                 end
            end
            Fsr_x(i,se(k,3))=sum(matchse_r(1,:));
            Fsr_y(i,se(k,3))=sum(matchse_r(2,:));
            Fsr_z(i,se(k,3))=sum(matchse_r(3,:));
            
        else
            Fsr_x(i,se(k,3))=-springF(i,k)*cos(se(k,4))*cos(se(k,5));
            Fsr_y(i,se(k,3))=-springF(i,k)*cos(se(k,4))*sin(se(k,5));
            Fsr_z(i,se(k,3))=-springF(i,k)*sin(se(k,4));
        end
        match_r(:,1)=0;
    end
    Fs_x(i,:,n)=Fsl_x(i,:)+Fsr_x(i,:);
    Fs_y(i,:,n)=Fsl_y(i,:)+Fsr_y(i,:);
    Fs_z(i,:,n)=Fsl_z(i,:)+Fsr_z(i,:);
end