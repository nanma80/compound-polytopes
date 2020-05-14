#!/usr/bin/env wolframscript

<< Quaternions`

phi = (1 + Sqrt[5])/2;
vertexJson = {{1,1},{-1,1},{-1,1},{-1,1}};
vertexQ = Apply[Quaternion, vertexJson.{1, phi}];

target = Quaternion[0, 0, 0, 1] * 8;
factor = Divide[target, vertexQ];
Print[factor];
Print[factor ** vertexQ];
Print[factor//N];

Print[vertexQ//N];

Print[Quaternion[0.618, 0.618, -0.618, 2.618] ** Quaternion[2.618, 0.618, 0.618, 0.618]]
