// Gmsh project created on Fri May 12 12:05:32 2023
SetFactory("OpenCASCADE");
//+
Point(1) = {0.25, 0.25, 0, 1.0};
//+
Point(2) = {0.25, 0, 0, 1.0};
//+
Point(3) = {0, 0.25, 0, 1.0};
//+
Point(4) = {1, 0, 0, 1.0};
//+
Point(5) = {0, 1, 0, 1.0};
//+
Point(6) = {1, 1, 0, 1.0};
//+
Point(7) = {1, 0.75, 0, 1.0};
//+
Point(8) = {0.75, 0.75, 0, 1.0};
//+
Point(9) = {0.75, 1, 0, 1.0};
//+
Point(10) = {0, 0, 0, 1.0};
//+
Line(1) = {5, 9};
//+
Line(2) = {8, 9};
//+
Line(3) = {6, 7};
//+
Line(4) = {8, 7};
//+
Line(5) = {6, 8};
//+
Line(6) = {9, 6};
//+
Line(7) = {7, 4};
//+
Line(8) = {4, 2};
//+
Line(9) = {2, 1};
//+
Line(10) = {1, 3};
//+
Line(11) = {3, 10};
//+
Line(12) = {10, 2};
//+
Line(13) = {10, 1};
//+
Line(14) = {1, 8};
//+
Line(15) = {5, 3};
//+
Curve Loop(1) = {15, -10, 14, 2, -1};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {5, 2, 6};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {4, -3, 5};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {7, 8, 9, 14, 4};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {9, -13, 12};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {11, 13, 10};
//+
Plane Surface(6) = {6};
//+
Transfinite Curve {15, 1, 2, 10, 8, 9, 7, 12, 11, 6, 3} = 10 Using Progression 1;
//+
Transfinite Curve {13, 5, 5} = 20 Using Progression 1;
//+
Transfinite Curve {14, 14} = 40 Using Progression 1;
//+
Transfinite Curve {4} = 10 Using Progression 1;
