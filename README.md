
# Create and draw geometry description for SiD:

```
cd geometry_scripts
root


.L SolGeomSiD.cxx
SolGeom()
SolGeom solGeomObj
solGeomObj.GeoPrint("GeoSiD.txt")
solGeomObj.Draw()

```

# Print SiD material budget

```
root

.L plot_materialSiD.c 
plot_materialSiD()
```

# Compare tracking resolution for IDEA, CLD and SiD

```
root
.L LoadAll.c
LoadAll()

.L CompGeom.c
CompGeom(50)  //50 deg is the track angle

```
