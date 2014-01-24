g++ -c Geometry3D.cpp -o geometry3D.o
g++ -c test_geometry.cpp -o test_geometry.o
g++ -lm geometry3D.o test_geometry.o
