#!/usr/bin/env wolframscript

numericalValue = True;

epsilon = 1.0*^-5;
dataFolder = "data";

<< Quaternions`

reflect[point_, mirror_] := FullSimplify[point - 2 * (point.mirror)/(mirror.mirror)*mirror];
doubleReflect[point_, mirror1_, mirror2_] := reflect[reflect[point, mirror1], mirror2];
normalizedInnerProd[u_, v_] := u.v/Norm[u]/Norm[v];
perspectiveProjection[point_, fullProjectionMatrix_, wDistance_ ] := projectToDirection[point, fullProjectionMatrix[[1;;3]]] / (1 + point.(fullProjectionMatrix[[4]])/wDistance);
projectToDirection[vertex4d_, projectionMatrix_] := projectionMatrix.vertex4d;
getEdgesFromFace[face_] := Table[{face[[i+1]], face[[Mod[i+1, Length[face]] + 1]]}, {i, 0, Length[face] - 1}];
sameVertex[point1_, point2_] := point1 == point2;
cubeVertices[dimension_] := 2 * Map[IntegerDigits[#, 2][[2;;dimension+1]]&, Table[k, {k, 2^dimension, 2 * 2^dimension - 1}]] - 1;
allPlusMinus[v_] := Union[Map[(v #)&, cubeVertices[Length[v]]]];
(* SameMatrix[m1_, m2_] := ( m1 == m2 ) || (Max[Max[N[Abs[m1 - m2]]]] < epsilon); *)
SameMatrix[m1_, m2_] := (Max[Max[N[Abs[m1 - m2]]]] < epsilon);
sameCenter[edge1_, edge2_] := SameMatrix[Total[edge1], Total[edge2]];
SameTotal[edge1_, edge2_] := SameMatrix[Total[edge1], Total[edge2]];
getQString[q_, separator_] := If[IntegerQ[q], IntegerString[q], IntegerString[Numerator[q]] <> separator <> IntegerString[Denominator[q]]];
getFileName[symbol_] := "complex_polytope_" <> IntegerString[symbol[[1]]] <> "_" <> getQString[symbol[[2]], "over"] <> "_" <> IntegerString[symbol[[3]]];
getShapeName[symbol_] := IntegerString[symbol[[1]]] <> "{" <> getQString[symbol[[2]], "/"] <> "}" <> IntegerString[symbol[[3]]];
maxNorm[direction_, realVertices_] := Max[Map[(#.direction)^2/direction.direction&, realVertices]];
getReal[complexVector_] := Flatten[ Map[{Re[#], Im[#]}&, complexVector] ];

shape120cell = "polytope_5_3_3";
Print["shape: " <> shape120cell];
polytope120cellFileName = FileNameJoin[{dataFolder, shape120cell <> ".wl"}];

polytope120cellData = Get[polytope120cellFileName];
polytopeVertices = polytope120cellData["vertices"];

shape600cell = "polytope_3_3_5";
polytope600cellFileName = FileNameJoin[{dataFolder, shape600cell <> ".wl"}];

polytope600cellData = Get[polytope600cellFileName];
polytope600cellVertices = polytope600cellData["vertices"];
Print["Loaded 600-cell data"];

seed1 = polytopeVertices[[600]];
baseVertices = Select[polytopeVertices, Simplify[normalizedInnerProd[seed1, #] + 1/4] == 0& ];

compounds = {};

For[seed2Index = 1, seed2Index <= Length[baseVertices], seed2Index++,
    seed2 = baseVertices[[seed2Index]];
    base2Vertices = Select[baseVertices, Simplify[normalizedInnerProd[seed2, #] + 1/4] == 0& ];
    templateVertices = Join[ {seed1, seed2}, base2Vertices ];
    AppendTo[compounds, templateVertices];
];

same5Cell[array1_, array2_] := Length[ Intersection[array1, array2, SameTest -> sameVertex] ] == 5;

compounds = Union[compounds, SameTest -> same5Cell];
(* length = 7 *)

Print[Length[compounds]];

For[i = 1, i<=Length[compounds], i++,
    templateVertices = compounds[[i]];
    distanceTable = Table[ Apply[normalizedInnerProd, templateVertices[[{i, j} ]]]//Simplify, {i, 5}, {j, 5}];
    Print[distanceTable == IdentityMatrix[5]*(5/4) - 1/4];
];


icosianLists = polytope600cellVertices/2;

quaternionRotation[leftList_, rightList_] := List @@ ( (Quaternion @@ leftList) ** (Quaternion @@ rightList) );

cellsAs = Table[
    Table[
        quaternionRotation[compounds[[k, j]], icosianLists[[i]]],
        {i, 1, Length[icosianLists]},
        {j, 1, Length[compounds[[1]]]}
    ]
    , {k, 1, Length[compounds]}
];

Print[Length/@cellsAs];

For[k=1, k<=Length[compounds], k++,
    Print[ Intersection[Flatten[cellsAs[[k]], 1], polytopeVertices, SameTest -> sameVertex]//Length ];
];

cellsBs = Table[
    Table[
        quaternionRotation[icosianLists[[i]], compounds[[k, j]]],
        {i, 1, Length[icosianLists]},
        {j, 1, Length[compounds[[1]]]}
    ]
    , {k, 1, Length[compounds]}
];

mirrorVector[v_] := {v[[1]], v[[2]], v[[3]], -v[[4]]};
mirrorCellsAs = Map[mirrorVector, cellsAs, {3}];

Print[Length/@cellsBs];

For[k=1, k<=Length[compounds], k++,
    Print[ Intersection[Flatten[cellsBs[[k]], 1], polytopeVertices, SameTest -> sameVertex]//Length ];
];

overlapTable = Table[ 
    Length[Intersection[cellsAs[[i]], cellsBs[[j]], SameTest -> same5Cell]], 
    {i, Length[cellsAs]}, {j, Length[cellsBs]}];

(* 
overlapTable =
{{120, 0, 0, 0, 0, 0, 0}, {0, 20, 20, 20, 20, 20, 20}, {0, 20, 20, 20, 20, 20, 20}, {0, 20, 20, 20, 20, 20, 20},
 {0, 20, 20, 20, 20, 20, 20}, {0, 20, 20, 20, 20, 20, 20}, {0, 20, 20, 20, 20, 20, 20}}

In[95]:= overlapTable//MatrixForm                                                                                       
Out[95]//MatrixForm= 120   0     0     0     0     0     0

                     0     20    20    20    20    20    20

                     0     20    20    20    20    20    20

                     0     20    20    20    20    20    20

                     0     20    20    20    20    20    20

                     0     20    20    20    20    20    20

                     0     20    20    20    20    20    20
 *)

mirrorOverlapTable = Table[ 
    Length[Intersection[mirrorCellsAs[[i]], cellsBs[[j]], SameTest -> same5Cell]], 
    {i, Length[cellsAs]}, {j, Length[cellsBs]}];

(* 
In[70]:= mirrorOverlapTable//InputForm                                                                                  
Out[70]//InputForm=
{{120, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 120}, {0, 0, 120, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 120, 0},
 {0, 0, 0, 0, 120, 0, 0}, {0, 120, 0, 0, 0, 0, 0}, {0, 0, 0, 120, 0, 0, 0}}

In[71]:= mirrorOverlapTable//MatrixForm                                                                                 
Out[71]//MatrixForm= 120   0     0     0     0     0     0

                     0     0     0     0     0     0     120

                     0     0     120   0     0     0     0

                     0     0     0     0     0     120   0

                     0     0     0     0     120   0     0

                     0     120   0     0     0     0     0

                     0     0     0     120   0     0     0

 *)

For[k = 1, k<= Length[cellsAs], k++,
    Print[Select[Range[120], same5Cell[cellsAs[[k, #]], mirrorCellsAs[[k, #]] ]&]//Length];
];
(* 20, 0, 0, 0, 0, 0, 0 *)
(* No McMullen 5-cell is the mirror of itself *)

For[k = 1, k<= Length[cellsAs], k++,
    Print[Length[
            Intersection[cellsAs[[k]], mirrorCellsAs[[k]], SameTest -> same5Cell]
        ]
    ];
];
(* 120, 20, 20, 20, 20, 20, 20 *)
