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

echo "3Doperations"
cd 3Doperations
cd rasterEllipsoid
python rasterEllipsoid.py
cd ..
cd repeatShape
python repeatShape.py
cd ..
cd voronoi
python voronoi.py
cd ..
cd voronoiGyroid
python voronoiGyroid.py
cd ..
cd ..

echo "Mesh"
if mmg3d_O3 -h; then
    echo "MMG3D is installed"
    
    cd Mesh
    cd mmg
    python test_mmg3d.py
    cd ..
    cd mmg-voro
    python test_mmg.py
    cd ..
    cd ..
else
    echo "MMG3D is not installed"
    exit 1
fi

echo "--TPMS--"
cd TPMS
echo "coordinate_system"
cd coordinate_system
python cylindrical.py
python cylindrical_graded.py
python spherical.py
cd ..
echo "grading"
cd grading
python cell_size.py
python cell_type.py
python density.py
cd ..
echo "gyroid"
cd gyroid
python gyroid.py
cd ..
echo "surface"
cd surface
python fischerKochS.py
cd ..
echo "tpms"
cd tpms
python tpms.py
cd ..
echo "tpmsShell"
cd tpmsShell
python tpmsShell.py
cd ..
echo "tpmsSphere"
cd tpmsSphere
python tpmsSphere.py
cd ..
cd ..
