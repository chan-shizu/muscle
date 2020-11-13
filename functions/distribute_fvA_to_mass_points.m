FvA(1:pointNum(1),1:3)=0;
tetranum=(1:tetraNum(1)*tetraNum(2)/5);
%a行b列目のnode番号を有するtetraNumに1
for m=1:tetraNum(1)/5
    for k=1:tetraNum(2)
        for j=1:tetraNum(1)/5
            matchA(1+(j-1)*tetraNum(2):tetraNum(2)*j,1)=tetra(2+5*(j-1),:).'==tetra(2+5*(m-1),k);    %同じ質点のある三角Noが　ある場合：1 ない：0
        end
        matchAnum=tetranum(matchA);   %同じ質点を有する三角Noを抽出
        mN=size(matchAnum);      %1つ以上なら同じ質点を有する三角Noがほかにあるということ
        if mN(2)>1
            matchf_A=[];
            %プログラム高速化のための事前割り当て
            matcn_A = zeros(3, mN(2));
            for n=1:mN(2)   %matchnumに入っている三角Noのxy方向力をmatchf_Aに表示
                matchf_A(1,n)=fvA(1,matchAnum(1,n));   %x
                matchf_A(2,n)=fvA(2,matchAnum(1,n));   %y
                matchf_A(3,n)=fvA(3,matchAnum(1,n));   %z
            end
            
            FvA(tetra(2+5*(m-1),k),1)=sum(matchf_A(1,:));
            FvA(tetra(2+5*(m-1),k),2)=sum(matchf_A(2,:));
            FvA(tetra(2+5*(m-1),k),3)=sum(matchf_A(3,:));
            
        else
            FvA(tetra(2+5*(m-1),k),1)=fvA(1,matchAnum);
            FvA(tetra(2+5*(m-1),k),2)=fvA(2,matchAnum);
            FvA(tetra(2+5*(m-1),k),3)=fvA(3,matchAnum);
        end
    end
end