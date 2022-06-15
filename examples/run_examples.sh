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
python honeycomb.py
cd ..
cd octetTruss
python octetTruss.py
cd ..
cd ..

echo "--TPMS--"
cd TPMS
echo "gyroid"
cd gyroid
python gyroid.py
cd ..
echo "tpms"
cd tpms
python tpms.py
cd ..
echo "tpmsShell"
cd tpmsShell
python tpmsShell.py
cd ..
# echo "tpmsSphere"
# cd tpmsSphere
# python tpmsSphere.py
# cd ..
cd ..

echo "3Doperations"
cd 3Doperations
cd rasterEllipsoid
python rasterEllipsoid.py
cd ..
cd repeatGeom
python repeatGeom.py
cd ..
cd voronoi
python voronoi.py
cd ..
cd voronoiGyroid
python voronoiGyroid.py
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