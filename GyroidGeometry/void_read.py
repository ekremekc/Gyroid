import matplotlib.pyplot as plt

path = '/home/ekrem/Dev/Gyroid/GyroidGeometry/Voids'
num_files = 20

# for asc_file in range(num_files):
# 	coords = []
# 	file_name = path+"/void"+str(asc_file)+".asc"
# 	with open(file_name) as f:
# 		
# 		x,y = [], []
# 		z = 0
# 		for l in f:
# 		    row = l.split()
# 		    x.append(float(row[0]))
# 		    y.append(float(row[1]))
coords = []         
for asc_file in range(num_files):
	
    file_name = path+"/void"+str(asc_file)+".asc"
    with open(file_name) as f:
		
        x,y = [], []
        z = 0
        for l in f:
            row = l.split()
            x.append(float(row[0]))
            y.append(float(row[1]))
            
        plt.plot(x,y)
        
# coords.append([x,y])




plt.show()