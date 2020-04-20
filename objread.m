function obj = objread(fname)
v = []; vt = []; vn = []; f.v = []; f.vt = []; f.vn = [];
fid = fopen(fname);
while 1    
    tline = fgetl(fid);
    if ~ischar(tline),   break,   end
     ln = sscanf(tline,'%s',1); % line type 
    switch ln
        case 'v'   % mesh vertexs
            v = [v; sscanf(tline(2:end),'%f')'];
        case 'vt'  % texture coordinate
            vt = [vt; sscanf(tline(3:end),'%f')'];
        case 'vn'  % normal coordinate
            vn = [vn; sscanf(tline(3:end),'%f')'];
        case 'f'   % face definition
            fv = []; fvt = []; fvn = [];
            str = textscan(tline(2:end),'%s'); str = str{1};
       
           nf = length(findstr(str{1},'/')); % nf=0:v, nf=1:v,vt, nf=2:v,vt or v,vn
           
           [tok str] = strtok(str,'/');     % vertex only
            for k = 1:length(tok) fv = [fv str2num(tok{k})]; end
           
            if (nf == 1) 
                [tok str] = strtok(str,'/');   % add texture coordinates
                for k = 1:length(tok) fvt = [fvt str2num(tok{k})]; end
            else    
                nff = length(findstr(str{1},'//')); % nff=0:v,vt,vn nff=1:v,vn
                if (nff == 0) 
                    [tok str] = strtok(str,'/');   % add texture coordinates
                    for k = 1:length(tok) fvt = [fvt str2num(tok{k})]; end
                    [tok str] = strtok(str,'/');   % add texture coordinates
                    for k = 1:length(tok) fvn = [fvn str2num(tok{k})]; end
                else
                    [tok str] = strtok(str,'//');   % add normal coordinates
                    for k = 1:length(tok) fvn = [fvn str2num(tok{k})]; end
                end
            end
            f.v = [f.v; fv]; f.vt = [f.vt; fvt]; f.vn = [f.vn; fvn];
    end
end
fclose(fid);
% set up matlab object 
obj.v = v; obj.vt = vt; obj.vn = vn; obj.f = f;
end