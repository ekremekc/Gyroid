import FreeCAD as App
import Draft
import os


doc = App.newDocument()
path = '/home/ekrem/Dev/Gyroid/GyroidGeometry/Voids'
num_files = len([f for f in os.listdir(path)if os.path.isfile(os.path.join(path, f))])

L = 50
h = 50

App.activeDocument().addObject('PartDesign::Body','Body000')
Gui.activeView().setActiveObject('pdbody', App.activeDocument().Body000)



#pl = FreeCAD.Placement()
#pl.Rotation.Q = (0.0, 0.0, 0, 1.0)
#pl.Base = FreeCAD.Vector(0.0, 0.0, 0.0)
#rec = Draft.makeRectangle(length=L, height=h, placement=pl, face=True, support=None)
#doc.recompute()

splines = []
for asc_file in range(num_files):
	coords = []
	file_name = path+"/void"+str(asc_file)+".asc"
	with open(file_name) as f:
		
		x,y = [], []
		z = 0
		for l in f:
		    row = l.split()
		    x.append(float(row[0]))
		    y.append(float(row[1]))
		    vector = App.Vector(float(row[0]),float(row[1]),z)
		    coords.append(vector)

	#splines.append(Draft.make_bspline(coords, closed=True))
	splines.append(Draft.make_wire(coords, closed=True, face=None))
	if asc_file ==0:
		App.getDocument('Unnamed').getObject('Wire').adjustRelativeLinks(App.getDocument('Unnamed').getObject('Body000'))
		App.getDocument('Unnamed').getObject('Body000').ViewObject.dropObject(App.getDocument('Unnamed').getObject('Wire'),None,'',[])
	elif asc_file <10:
		App.getDocument('Unnamed').getObject('Wire00'+str(asc_file)).adjustRelativeLinks(App.getDocument('Unnamed').getObject('Body000'))
		App.getDocument('Unnamed').getObject('Body000').ViewObject.dropObject(App.getDocument('Unnamed').getObject('Wire00'+str(asc_file)),None,'',[])
	else:
		App.getDocument('Unnamed').getObject('Wire0'+str(asc_file)).adjustRelativeLinks(App.getDocument('Unnamed').getObject('Body000'))
		App.getDocument('Unnamed').getObject('Body000').ViewObject.dropObject(App.getDocument('Unnamed').getObject('Wire0'+str(asc_file)),None,'',[])



doc.recompute()
Gui.SendMsgToActiveView("ViewFit")
