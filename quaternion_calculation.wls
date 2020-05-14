#!/usr/bin/env wolframscript

<< Quaternions`

phi = (1 + Sqrt[5])/2;
vertexJson = {{1,1},{-1,1},{-1,1},{-1,1}};
vertexJsonRotated = vertexJson[[{2,3,4,1}]];
vertexQ = Apply[Quaternion, vertexJsonRotated.{1, phi}];

target = Quaternion[0, 0, 0, 1] * 8;
factor = Divide[target, vertexQ];
Print["product"];
Print[factor ** vertexQ];

Print["factor"];
Print[factor];
Print[factor//N];

Print["original"];
Print[vertexQ//N];

