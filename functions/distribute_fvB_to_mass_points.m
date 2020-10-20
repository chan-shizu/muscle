FvB(1:pointNum(1),1:3)=0;
%a行b列目のnode番号を有するtetraNumに1
for m=1:tetraNum(1)/5
    for k=1:tetraNum(2)
        for j=1:tetraNum(1)/5
            matchB(1+(j-1)*tetraNum(2):tetraNum(2)*j,1)=tetra(3+5*(j-1),:).'==tetra(3+5*(m-1),k);    %同じ質点のある三角Noが　ある場合：1 ない：0
        end
        matchBnum=tetranum(matchB);   %同じ質点を有する三角Noを抽出
        mN=size(matchBnum);      %1つ以上なら同じ質点を有する三角Noがほかにあるということ
        if mN(2)>1
            matchf_B=[];
            %プログラム高速化のための事前割り当て
            matcn_B = zeros(3, mN(2));
            for n=1:mN(2)  %matchnumに入っている三角Noのxy方向力をmatchf_Bに表示
                matchf_B(1,n)=fvB(1,matchBnum(1,n));   %x
                matchf_B(2,n)=fvB(2,matchBnum(1,n));   %y
                matchf_B(3,n)=fvB(3,matchBnum(1,n));   %z
            end
            
            FvB(tetra(3+5*(m-1),k),1)=sum(matchf_B(1,:));
            FvB(tetra(3+5*(m-1),k),2)=sum(matchf_B(2,:));
            FvB(tetra(3+5*(m-1),k),3)=sum(matchf_B(3,:));
            
        else
            FvB(tetra(3+5*(m-1),k),1)=fvB(1,matchBnum);
            FvB(tetra(3+5*(m-1),k),2)=fvB(2,matchBnum);
            FvB(tetra(3+5*(m-1),k),3)=fvB(3,matchBnum);
        end
    end
end