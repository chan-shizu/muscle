mybox1 = Box heightsegs:5 lengthsegs:4 widthsegs:4 
convertTo mybox1 TriMeshGeometry
mybox1.verts.count

xdata = (dotnetClass "System.IO.File").ReadAllLines "C:\Users\shizu\OneDrive\ドキュメント\研究matlab\output\3dsmax_test_column_x.csv"
ydata = (dotnetClass "System.IO.File").ReadAllLines "C:\Users\shizu\OneDrive\ドキュメント\研究matlab\output\3dsmax_test_column_y.csv"
zdata = (dotnetClass "System.IO.File").ReadAllLines "C:\Users\shizu\OneDrive\ドキュメント\研究matlab\output\3dsmax_test_column_z.csv"
print xdata.count
	
animate on
(	
	for n=1 to 500 by 50 do
	(
		xInfo = (filterString xdata[n] ",")
		yInfo = (filterString ydata[n] ",")
		zInfo = (filterString zdata[n] ",")
		
		/*if (n == 1) then
		(
			at time n
		)
		else
		(
			at time n
		)*/
			
		(	
			for i = 1 to xInfo.count do 
			(
				case  of
				(
					(i < 26):
					(
						print(1000)
						mybox1.verts[i].pos.x  = (xInfo[i] as float)
						--print scale_x
						mybox1.verts[i].pos.y = (yInfo[i] as float)
						--print scale_y
						mybox1.verts[i].pos.z = (zInfo[i] as float)
						--print scale_z
					)
					(i < 51):	
					(
						print(2000)
						mybox1.verts[i].pos.x  = (xInfo[i+64] as float)
						--print scale_x
						mybox1.verts[i].pos.y = (yInfo[i+64] as float)
						--print scale_y
						mybox1.verts[i].pos.z = (zInfo[i+64] as float)
						--print scale_z
					)
					default:	
					(
						print(3000)
						mybox1.verts[i].pos.x = (xInfo[i-25] as float)
						--print scale_x
						mybox1.verts[i].pos.y = (yInfo[i-25] as float)
						--print scale_y
						mybox1.verts[i].pos.z = (zInfo[i-25] as float)
						--print scale_z
					)
				)
			)
		)
	)
)