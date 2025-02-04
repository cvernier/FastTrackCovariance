
# Create and draw geometry description for SiD:

```
cd geometry_scripts
root


# SolGeomSiD.cxx (or the one for IDEA or CLD should be your starting point)
# Modify the geometry therein as needed and print the txt file with the desired name
# You will use that txt file later to compare the resolutions for different detectors in CompGeom

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

# Now call CompGeom as in this example: 


CompGeom(45, 
         "GeoIDEA_BASE", "GeoCLD", "GeoSiD",  // geometry file names (without .txt)
         "IDEA", "CLD", "SiD",                 // detector names
         2.0, 2.0, 5.0);                       // magnetic fields


```
