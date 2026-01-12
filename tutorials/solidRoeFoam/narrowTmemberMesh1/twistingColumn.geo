//+
//Point(1) = {0, 0, 0, 0.0};
//+
//Point(2) = {48, 44, 0, 0.0};
//+
//Point(3) = {48, 60, 0, 0.0};
//+
//Point(4) = {0, 44, 0, 0.0};
//+

Point(1) = {0, 0, 0, 1.0};
Point(2) = {48, 44, 0, 1.0};
Point(3) = {48, 60, 0, 1.0};
Point(4) = {0, 44, 0, 1.0};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(5) = {3, 4, 1, 2};
Plane Surface(6) = {5};
Extrude {0, 0, 1} {
  Surface{6}; 
}

Transfinite Line {9, 4, 11, 2} = 20  Using Progression 1; //left right
Transfinite Line {8, 3, 10, 1} = 20  Using Progression 1; //top bottom 
Transfinite Line {13, 22, 14,18 } = 1  Using Progression 1; // z-axix
//Transfinite Line {8, 9, 10, 11, 22, 13, 2, 1 , 3, 14, 4, 18} = 4  Using Progression 1;
Transfinite Surface {28};
Transfinite Surface {27};
Transfinite Surface {19};
Transfinite Surface {15};
Transfinite Surface {23};
Transfinite Surface {6};

Recombine Surface {28, 27, 15, 6, 19, 23};

Transfinite Volume{1} = {10, 14, 5, 6, 1, 2, 3, 4};

Physical Surface("sides") = {28, 6};
Physical Surface("left") = {19};
Physical Surface("bottom") = {23};
Physical Surface("right") = {27};
Physical Surface("top") = {15};
//+
Physical Volume("volume") = {1};
