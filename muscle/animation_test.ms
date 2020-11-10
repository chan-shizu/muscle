mybox1 = Box heightsegs:5 lengthsegs:4 widthsegs:4 
convertTo mybox1 TriMeshGeometry
animate on
(
	for i=1 to 2 do
	(
		at time i (
			--mybox1.pos = [i*10, i*20, i*30]
			mybox1.verts[i].pos.x  = i
			mybox1.scale = [i, i, i]
		)
	)
)