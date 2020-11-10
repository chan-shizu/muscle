FvC(1:pointNum(1),1:3)=0;
%a行b列目のnode番号を有するtetraNumに1
for m=1:tetraNum(1)/5
    for k=1:tetraNum(2)
        for j=1:tetraNum(1)/5
            matchC(1+(j-1)*tetraNum(2):tetraNum(2)*j,1)=tetra(4+5*(j-1),:).'==tetra(4+5*(m-1),k);    %同じ質点のある三角Noが　ある場合：1 ない：0
        end
        matchCnum=tetranum(matchC);   %同じ質点を有する三角Noを抽出
        mN=size(matchCnum);      %1つ以上なら同じ質点を有する三角Noがほかにあるということ
        if mN(2)>1
            matchf_C=[];
            %プログラム高速化のための事前割り当て
            matcn_C = zeros(3, mN(2));
            for n=1:mN(2)   %matchnumに入っている三角Noのxy方向力をmatchf_Cに表示
                matchf_C(1,n)=fvC(1,matchCnum(1,n));   %x
                matchf_C(2,n)=fvC(2,matchCnum(1,n));   %y
                matchf_C(3,n)=fvC(3,matchCnum(1,n));   %z
            end
            
            FvC(tetra(4+5*(m-1),k),1)=sum(matchf_C(1,:));
            FvC(tetra(4+5*(m-1),k),2)=sum(matchf_C(2,:));
            FvC(tetra(4+5*(m-1),k),3)=sum(matchf_C(3,:));
            
        else
            FvC(tetra(4+5*(m-1),k),1)=fvC(1,matchCnum);
            FvC(tetra(4+5*(m-1),k),2)=fvC(2,matchCnum);
            FvC(tetra(4+5*(m-1),k),3)=fvC(3,matchCnum);
        end
    end
end