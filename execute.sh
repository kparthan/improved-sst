#g++ -c Geometry3D.cpp -o geometry3D.o -std=c++11
g++ -c test_geometry.cpp -o test_geometry.o -std=c++11
g++ -lm test_geometry.o #geometry3D.o
