/*mybox1 = Box heightsegs:5 lengthsegs:4 widthsegs:4 
convertTo mybox1 TriMeshGeometry
mybox1.verts.count
--myVert1 = mybox1.verts[1].pos*/
/*for i=28 to i=35 do
(
	print mybox1.verts[i].pos 
)*/

--move $Box001.verts[#{2, 4, 6}] [0,-1.48116,0]

--$Box001.selectedVerts[1].pos
--move $Box001.selectedVerts [0,-13.3767,0]
/*adata = (dotnetClass "System.IO.File").ReadAllLines "C:\Users\shizu\OneDrive\ドキュメント\研究matlab\muscle\L_Digastric_Posterior_min_n_3dsmax.csv"
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
)*/

mybox1 = Box heightsegs:11 lengthsegs:4 widthsegs:4 
convertTo mybox1 TriMeshGeometry
mybox1.verts.count

/*各時刻における筋肉の各質点の座標を入力*/
xdata = (dotnetClass "System.IO.File").ReadAllLines "D:\Documents\matlabまとめ\研究matlab\3ds_max\biceps_ver3_data_x.csv"
ydata = (dotnetClass "System.IO.File").ReadAllLines "D:\Documents\matlabまとめ\研究matlab\3ds_max\biceps_ver3_data_y.csv"
zdata = (dotnetClass "System.IO.File").ReadAllLines "D:\Documents\matlabまとめ\研究matlab\3ds_max\biceps_ver3_data_z.csv"
print xdata.count


		/*行ごとに読み取り*/
		xInfo = (filterString xdata[1] ",")
		yInfo = (filterString ydata[1] ",")
		zInfo = (filterString zdata[1] ",")
		
		/*if (n == 1) then
		(
			at time n
		)
		else
		(
			at time n
		)*/
			

		/*for i = 1 to xInfo.count do*/ 
for i = 1 to 300 do 
(
	if i < 26 then
		(
			print(1000)
			mybox1.verts[i].pos.x  = (xInfo[i] as float)*100
		--print scale_x
			mybox1.verts[i].pos.y = (yInfo[i] as float)*100
			--print scale_y
			mybox1.verts[i].pos.z = (zInfo[i] as float)*100
			--print scale_z
		)
	else if i < 51 then
		(
			print(2000)
		mybox1.verts[i].pos.x  = (xInfo[i+160] as float)*100
			--print scale_x
			mybox1.verts[i].pos.y = (yInfo[i+160] as float)*100
		--print scale_y
			mybox1.verts[i].pos.z = (zInfo[i+160] as float)*100
			--print scale_z
		)
	else	
		(
			print(3000)
			mybox1.verts[i].pos.x = (xInfo[i-25] as float)*100
			--print scale_x
			mybox1.verts[i].pos.y = (yInfo[i-25] as float)*100
			--print scale_y
			mybox1.verts[i].pos.z = (zInfo[i-25] as float)*100
			--print scale_z
		)
	)
)

