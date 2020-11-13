clearlistener()
adata = (dotnetClass "System.IO.File").ReadAllLines "C:\Users\shizu\OneDrive\ÉhÉLÉÖÉÅÉìÉg\å§ãÜmatlab\muscle\L_Digastric_Anterior_min_n.csv"

for i = 1 to adata.count do 
(
	if i==1 then (
		pInfo = (filterString adata[i] ",")
		print pInfo
		scale_x = (pInfo[1] as float)*1000
		print scale_x
		scale_y = (pInfo[2] as float)*1000
		print scale_y
		scale_z = (pInfo[3] as float)*1000
		print scale_z
		mybox1 = box()
		mybox1.pos = [scale_x , scale_y, scale_z]
		/*s = sphere segments:(pInfo[6] as integer)
		s.radius = pInfo[2] as integer 
		s.pos = [pInfo[3] as integer,pInfo[4] as integer, pInfo[5] as integer]
		s.name = pInfo[1]*/
	)
	
	if i==53 then (
		pInfo = (filterString adata[i] ",")
		print pInfo
		scale_x = (pInfo[1] as float)*1000
		print scale_x
		scale_y = (pInfo[2] as float)*1000
		print scale_y
		scale_z = (pInfo[3] as float)*1000
		print scale_z
		mybox53 = box()
		mybox53.pos = [scale_x , scale_y, scale_z]
		/*s = sphere segments:(pInfo[6] as integer)
		s.radius = pInfo[2] as integer 
		s.pos = [pInfo[3] as integer,pInfo[4] as integer, pInfo[5] as integer]
		s.name = pInfo[1]*/
	)
)