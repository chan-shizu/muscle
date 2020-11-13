clear
%%
%obj�t�@�C���̓����Ă���t�H���_�̈ʒu���w�肷��
dataPath = 'obj\';

%�؂̖��O���w�肷��
objList = [
    {
    'biceps'
    };
];
%dataPath�̒��ɂ���AobjList�ɂ��閼�O��obj�t�@�C����ǉ����Ă���

%�o�͐��ύX����
global outputPath
% outputpath_name = strcat('outputVrml\', objList(1),'.wrl')
% outputPath = outputpath_name
outputPath = 'outputVrml\outputfile_L_Digastric_Anterior_notepc.wrl';


%%
% %�t�@�C�����̑Sobj�𓝍�����ꍇ
% cd obj
%     list = dir;
% cd ..
% p = 0;
% for q = 1:size(list, 1)
%     if contains((list(q).name), '.obj')
%         p = p+1;
%         tmp =  char(list(q).name);
%         objList{p, 1} = tmp(1:end-4);
%     end
% end

%%
%loadVRML
inputVrmlFile = 'inputVrml\Develop_Bone_only.wrl';%input�ƂȂ鍜�̈ʒu���w��
fileID = fopen(inputVrmlFile, 'r');
nLine = 0;
while 1
    tmp = fgetl(fileID);
    if tmp==-1
        break;
    else
        nLine = nLine+1;
        inputVrml{nLine,1} = cellstr(tmp);
    end
end
fclose(fileID);
disp('load vrml')

%%
nFiles = size(objList, 1);
muscleVrml = [];
for q = 1:nFiles
    muscleName = objList(q);
    objPath = cell2mat(strcat(dataPath, muscleName, '.obj'));
    obj = objread(objPath);

    %Header
    vrmlHeader = returnHeader(muscleName, obj);
    
    vertices = obj.v;
    vrmlVertices = [];
    for p = 1:size(vertices, 1)%�Ō�̍s�����R���}��t���Ȃ�(�ȉ�����)
        if p ~= size(vertices, 1)
            vrmlVertices = [vrmlVertices; cellstr(cell2mat(strcat([{'               '}, num2str(vertices(p, 1)), {' '}, num2str(vertices(p, 2)), {' '}, num2str(vertices(p,3)), ','])))];
        else
            vrmlVertices = [vrmlVertices; cellstr(cell2mat(strcat([{'               '}, num2str(vertices(p, 1)), {' '}, num2str(vertices(p, 2)), {' '}, num2str(vertices(p,3))])))];
        end
    end
    
    %header, footer�ȊO�̋L�q����
    vrmlMiddle1 = [
        {cellstr('            ]')};
        {cellstr('          }')};
        {cellstr(strcat({'          # '}, num2str(size(obj.f.v, 1)), {' polygons'}))};
        {cellstr('          coordIndex [')};
        ];
        
    fv = obj.f.v-1;%�C���f�b�N�X�ނ�0�X�^�[�g
    vrmlCoordIdx = [];
    for p = 1:size(fv, 1)
        if p ~= size(vertices, 1)
            vrmlCoordIdx = [vrmlCoordIdx; cellstr(cell2mat(strcat([{'            '}, num2str(fv(p, 1)), {' '}, num2str(fv(p, 2)), {' '}, num2str(fv(p,3)), ', -1,'])))];
        else
            vrmlCoordIdx = [vrmlCoordIdx; cellstr(cell2mat(strcat([{'            '}, num2str(fv(p, 1)), {' '}, num2str(fv(p, 2)), {' '}, num2str(fv(p,3)), ', -1'])))];
        end 
    end
    
    vrmlMiddle2 = [
        {cellstr('          ]')};
        {cellstr('          normal Normal {')};
        {cellstr(strcat({'            vector [   #'}, num2str(size(obj.vn, 1)), {' normals'}))};
        ];
        
    normals = obj.vn;
    vrmlNormals = [];
    for p = 1:size(vertices, 1)
        if p ~= size(vertices, 1)
            vrmlNormals = [vrmlNormals; cellstr(cell2mat(strcat([{'               '}, num2str(normals(p, 1)), {' '}, num2str(normals(p, 2)), {' '}, num2str(normals(p,3)), ','])))];
        else
            vrmlNormals = [vrmlNormals; cellstr(cell2mat(strcat([{'               '}, num2str(normals(p, 1)), {' '}, num2str(normals(p, 2)), {' '}, num2str(normals(p,3))])))];
        end
    end
    
    vrmlMiddle3 = [
        {cellstr('            ]')};
        {cellstr('          }')};
        {cellstr('          normalIndex [')};
        ];
        
    fvn = obj.f.vn-1;%�C���f�b�N�X��0�X�^�[�g
    vrmlNormalIdx = [];
    for p = 1:size(fv, 1)
        if p ~= size(vertices, 1)
            vrmlNormalIdx = [vrmlNormalIdx; cellstr(cell2mat(strcat([{'            '}, num2str(fvn(p, 1)), {' '}, num2str(fvn(p, 2)), {' '}, num2str(fvn(p,3)), ', -1,'])))];
        else
            vrmlNormalIdx = [vrmlNormalIdx; cellstr(cell2mat(strcat([{'            '}, num2str(fvn(p, 1)), {' '}, num2str(fvn(p, 2)), {' '}, num2str(fvn(p,3)), ', -1'])))];
        end
    end
    
    %Footer
    vrmlFooter = returnFooter(muscleName);

    %����
    muscleVrml = vertcat(muscleVrml, vrmlHeader, vrmlVertices, vrmlMiddle1, vrmlCoordIdx, vrmlMiddle2, vrmlNormals, vrmlMiddle3, vrmlNormalIdx, vrmlFooter);
    disp(strcat(muscleName, ' : converted'))
end
%%
%input�̍�vrml�Ƌ؂�vrml�ɕϊ��������̂�����������
vrml = vertcat(inputVrml, muscleVrml);

%%
%vrml�̏o��
writevrml(vrml);

%%
function VrmlHeader = returnHeader(muscleName, obj)
    VrmlHeader = [
    {cellstr('Transform {')};
    {cellstr('  children [')};
    {cellstr(strcat({'  DEF '}, muscleName, ' Transform {'))};
    {cellstr('   translation 0 0 0')};
    {cellstr('   rotation -1 0 0 1.5708')};%���W�A����90�xx�����}�C�i�X�����ɉ�]
    {cellstr('   scale 0.001 0.001 0.001')};
    {cellstr('   children [')};
    {cellstr('    Transform { # pivot node')};
    {cellstr('    translation 0 0 0')};
    {cellstr('    children [')};
    {cellstr(strcat({'    DEF '}, muscleName, '_geometry Group {'))};
    {cellstr('     children [')};
    {cellstr('      Shape {')};
    {cellstr('        appearance Appearance {')};
    {cellstr('          material DEF mtl1 Material {')};
    {cellstr('            diffuseColor 0.4 0.4 0.4')};
    {cellstr('            ambientIntensity 0.3')};
    {cellstr('            specularColor 0.7 0.7 0.7')};
    {cellstr('            emissiveColor 0 0 0')};
    {cellstr('            shininess 0.25')};
    {cellstr('            transparency 0')};
    {cellstr('          }')};
    {cellstr('        }')};
    {cellstr('        geometry IndexedFaceSet {')};
    {cellstr('          coord Coordinate {')};
    {cellstr(strcat({'            point [   #'}, num2str(size(obj.v, 1)), ' vertices'))};
    ];
end

function VrmlFooter = returnFooter(muscleName)
    VrmlFooter = [
        {cellstr('          ]')};
        {cellstr('          normalPerVertex TRUE')};
        {cellstr('          solid           FALSE')};
        {cellstr('        }  # End of IndexedFaceSet')};
        {cellstr('      }  # End of Shape node')};
        {cellstr('     ]    # End of geometry children nodes')};
        {cellstr(strcat({'    }    # End of '''}, muscleName, ''' Group node'))};
        {cellstr('    ]    # End of children nodes')};
        {cellstr(strcat({'    }    # End of pivot node for'}, muscleName))};
        {cellstr('   ]    # End of children nodes')};
        {cellstr(strcat({'  }    # End of '''}, muscleName, ''' Transform node'))};
        {cellstr('  ]    # End of children nodes')};
        {cellstr('}    # End of outer Transform node')};
        ];
end

function vrmlData = obj2vrml(data)
    m = size(data,2);
    if m~=3
        error('obj���b�V�����O�p�`�ł͂���܂���B')
    end
    vrmlData = [];
    for p = 1:size(data, 1)
        if p ~= size(data, 1)
            vrmlData = [vrmlData; cellstr(cell2mat(strcat([num2str(data(p, 1)), {' '}, num2str(data(p, 2)), {' '}, num2str(data(p,3)), ','])))];
        else
            vrmlData = [vrmlData; cellstr(cell2mat(strcat([num2str(data(p, 1)), {' '}, num2str(data(p, 2)), {' '}, num2str(data(p,3))])))];
        end
    end
end

function writevrml(vrml)
%�o��
    global outputPath
    nOutput = size(vrml, 1);
    outputID = fopen(outputPath, 'w');
    for n = 1:nOutput
        if iscell(vrml{n,1})
            fprintf(outputID, cell2mat(vrml{n,1}));
        else
            fprintf(outputID, (vrml{n,1}));
        end
    fprintf(outputID, '\n');
    end
    fclose(outputID);

    disp("vrml output")
end