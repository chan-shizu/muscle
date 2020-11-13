FvD(1:pointNum(1),1:3)=0;
%a行b列目のnode番号を有するtetraNumに1
for m=1:tetraNum(1)/5
    for k=1:tetraNum(2)
        for j=1:tetraNum(1)/5
            matchD(1+(j-1)*tetraNum(2):tetraNum(2)*j,1)=tetra(5+5*(j-1),:).'==tetra(5+5*(m-1),k);    %同じ質点のある三角Noが　ある場合：1 ない：0
        end
        matchDnum=tetranum(matchD);   %同じ質点を有する三角Noを抽出
        mN=size(matchDnum);      %1つ以上なら同じ質点を有する三角Noがほかにあるということ
        if mN(2)>1
            matchf_D=[];
            %高速化のための事前割り当て
            matcn_D = zeros(3, mN(2));
            for n=1:mN(2)   %matchnumに入っている三角Noのxy方向力をmatchf_Dに表示
                matchf_D(1,n)=fvD(1,matchDnum(1,n));   %x
                matchf_D(2,n)=fvD(2,matchDnum(1,n));   %y
                matchf_D(3,n)=fvD(3,matchDnum(1,n));   %z
            end
            
            FvD(tetra(5+5*(m-1),k),1)=sum(matchf_D(1,:));
            FvD(tetra(5+5*(m-1),k),2)=sum(matchf_D(2,:));
            FvD(tetra(5+5*(m-1),k),3)=sum(matchf_D(3,:));
            
        else
            FvD(tetra(5+5*(m-1),k),1)=fvD(1,matchDnum);
            FvD(tetra(5+5*(m-1),k),2)=fvD(2,matchDnum);
            FvD(tetra(5+5*(m-1),k),3)=fvD(3,matchDnum);
        end
    end
end