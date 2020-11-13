mybox1 = Box heightsegs:11 lengthsegs:4 widthsegs:4 
convertTo mybox1 TriMeshGeometry
mybox1.verts.count

/*各時刻における筋肉の各質点の座標を入力*/
xdata = (dotnetClass "System.IO.File").ReadAllLines "D:\Documents\matlabまとめ\研究matlab\3ds_max\biceps_ver3_data_x.csv"
ydata = (dotnetClass "System.IO.File").ReadAllLines "D:\Documents\matlabまとめ\研究matlab\3ds_max\biceps_ver3_data_y.csv"
zdata = (dotnetClass "System.IO.File").ReadAllLines "D:\Documents\matlabまとめ\研究matlab\3ds_max\biceps_ver3_data_z.csv"
print xdata.count

with animate on
(	
	
	for n=1 to 100  do
	(
		at time n mybox1.pos.x = n
		/*行ごとに読み取り*/
		xInfo = (filterString xdata[3*n] ",")
		yInfo = (filterString ydata[3*n] ",")
		zInfo = (filterString zdata[3*n] ",")
		
		/*if (n == 1) then
		(
			at time n
		)
		else
		(
			at time n
		)*/
				
		for i = 1 to xInfo.count do 
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
)