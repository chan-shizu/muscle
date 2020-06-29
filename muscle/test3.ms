mybox1 = Box heightsegs:5 lengthsegs:4 widthsegs:4 
convertTo mybox1 TriMeshGeometry
mybox1.verts.count
--myVert1 = mybox1.verts[1].pos
/*for i=28 to i=35 do
(
	print mybox1.verts[i].pos 
)*/

--move $Box001.verts[#{2, 4, 6}] [0,-1.48116,0]

--$Box001.selectedVerts[1].pos
--move $Box001.selectedVerts [0,-13.3767,0]
adata = (dotnetClass "System.IO.File").ReadAllLines "C:\Users\shizu\OneDrive\ÉhÉLÉÖÉÅÉìÉg\å§ãÜmatlab\muscle\L_Digastric_Posterior_min_n_3dsmax.csv"
print adata.count
for i =  1 to adata.count do 
(
	pInfo = (filterString adata[i] ",")
	mybox1.verts[i].pos.x  = (pInfo[1] as float)*1000
	--print scale_x
	mybox1.verts[i].pos.y = (pInfo[2] as float)*1000
	--print scale_y
	mybox1.verts[i].pos.z = (pInfo[3] as float)*1000
	--print scale_z
)