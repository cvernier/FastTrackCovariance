
# Create and draw geometry description for SiD:

```
cd geometry_scripts
root -l


# SolGeomSiD.cxx (or the one for IDEA or CLD should be your starting point)
# Modify the geometry therein as needed and print the txt file with the desired name
# You will use that txt file later to compare the resolutions for different detectors in CompGeom

.L SolGeomIDEA.cxx
SolGeom()
SolGeom solGeomObj
solGeomObj.GeoPrint("GeoIDEA_test.txt")
solGeomObj.Draw()

# Open the txt file (located in geometry_scripts/, inspect it and modify as necessary, e.g. if you want to remove some parts)
# !!! Then, make sure to move it to geometry_files/ ! This is where it will be read from when running the CompGeom scripts!

```

# Print SiD material budget

```
root

.L plot_materialSiD.c 
plot_materialSiD()
```

# Compare tracking resolution for IDEA, CLD and SiD

```
root -l
.L LoadAll.c
LoadAll("")

.L CompGeom.c

# Now call CompGeom as in this example: 


CompGeom(45, 
         "GeoIDEA_BASE", "GeoCLD", "GeoSiD",  // geometry file names (without .txt)
         "IDEA", "CLD", "SiD",                 // detector names for plots
         2.0, 2.0, 5.0);                       // magnetic fields

CompGeom(45,"GeoIDEA_BASE", "GeoIDEA_VTX", "GeoIDEA_GT", "old", "new", "b", 2., 2., 2.)



```
