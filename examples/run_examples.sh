# exit when any command fails
set -e

echo "BasicShapes"
cd BasicShapes
cd shapes
python shapes.py
cd ..
cd platon
python platon.py
cd ..
cd ..

echo "Lattices"
cd Lattices
cd honeycomb
python extrudedHoneycomb.py
cd ..
cd octetTruss
python testOctet.py
cd ..
cd ..

echo "TPMS"
cd TPMS
cd gyroid
python gyroid.py
cd ..
cd tpms
python tpms.py
cd ..
cd tpms_shell
python tpms_shell.py
cd ..
cd tpms_sphere
python tpms_sphere.py
cd ..
cd ..

echo "3Doperations"
cd 3Doperations
cd rasterEllipsoid
python rasterEllipsoid.py
cd ..
cd repeatGeom
python repeat.py
cd ..
cd Voronoi
python testNeper.py
cd ..
cd VoronoiGyroid
python voronoi_gyroid.py
cd ..
cd ..

echo "Mesh"
cd Mesh
cd mmg
python test_mmg.py
cd ..
cd mmg-voro
python test_mmg.py
cd ..
cd ..