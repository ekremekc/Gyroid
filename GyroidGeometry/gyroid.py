import matplotlib.pyplot as plt

from numpy import pi, cos, sin, linspace, meshgrid


n = 5

start_x = 0
L_x = 50
total_x = start_x+L_x

start_y = 2.5
L_y = 50
total_y = start_y+L_y


t = 0

x = linspace(start_x,total_x,250)
y = linspace(start_y,total_y,250)


def gyroid(x, y, n, L_x,L_y, t):
    tanım1 = (sin(2*pi*n*x/L_x) * cos(2*pi*n*y/L_y) + sin(2*pi*n*y/L_y) + cos(2*pi*n*x/L_x))
    return tanım1*tanım1 - t**2

XX, YY = meshgrid(x, y)
z = gyroid(XX, YY, n, L_x,L_y, t)

thickness = 0.15
contour = plt.contour(XX, YY, z,levels=[thickness])

# plt.colorbar()
plt.savefig('gyroid.pdf')
plt.show()

# dat0= contour.allsegs[0][0]
# plt.plot(dat0[:,0],dat0[:,1])
# plt.show()


# fig = plt.figure()
# ax = plt.axes(projection='3d')
# ax.contour3D(XX, YY, z, 100, cmap='jet')
# ax.set_xlabel('x')
# ax.set_ylabel('y')
# ax.set_zlabel('z')
# plt.savefig("3d.pdf")

top    = []
bottom = []
right  = []
left   = []

for voids in contour.allsegs:
    for index,void in enumerate(voids):
        with open('Voids/void'+str(index)+".asc", 'w') as fh:
            for x, y in void: 
                if x==start_x:
                    bottom.append([x,y])
                elif x==total_x:
                    right.append([x,y])
                elif y==start_y:
                    left.append([x,y])
                elif y==total_y:
                    top.append([x,y])
                fh.write('{} {}\n'.format(x, y))
                
        X, Y = void[:,0], void[:,1]
        plt.plot(X, Y)
plt.show()

# plt.plot(void[:,0], void[:,1])
# plt.show()

# def write_edge(liste, filename):
#     # Sort the coordinates first
#     if liste[0][0]==liste[1][0]:
#         liste = sorted(liste, key=lambda x: x[1])
#     else:
#         liste = sorted(liste, key=lambda x: x[0])
#     #Then write them into a file to be able to read from FreeCAD
#     with open('Edges/'+filename+".asc", 'w') as edge:
#         for x, y in liste: 
#             edge.write('{} {}\n'.format(x, y))
         


# write_edge(top, "top")
# write_edge(bottom, "bottom")
# write_edge(right, "right")
# write_edge(left, "left")











