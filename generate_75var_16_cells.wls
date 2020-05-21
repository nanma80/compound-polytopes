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
baseVertices = Select[polytopeVertices, Simplify[normalizedInnerProd[seed1, #]] == 0& ];

compounds = {};

For[seed2Index = 1, seed2Index <= Length[baseVertices], seed2Index++,
    seed2 = baseVertices[[seed2Index]];
    base2Vertices = Select[baseVertices, Simplify[normalizedInnerProd[seed2, #]] == 0& ];
    templateVertices = Join[ {seed1, -seed1, seed2, -seed2}, base2Vertices ];
    AppendTo[compounds, templateVertices];
];


same5Cell[array1_, array2_] := Length[ Intersection[array1, array2, SameTest -> sameVertex] ] == Length[array1];

compounds = Union[compounds, SameTest -> same5Cell];
(* length == 9, from a common vertex *)
(* 600 * 9 / 8 == 675 *)

compounds = compounds[[{1}]];

icosianLists = polytope600cellVertices/2;

quaternionRotation[leftList_, rightList_] := List @@ ( (Quaternion @@ leftList) ** (Quaternion @@ rightList) );

cellsA = 
    Table[
        quaternionRotation[compounds[[1, j]], icosianLists[[i]]],
        {i, 1, Length[icosianLists]},
        {j, 1, 8}
    ];

cellsA = Union[cellsA, SameTest -> same5Cell];

Print[ Intersection[Flatten[cellsA, 1], polytopeVertices, SameTest -> sameVertex]//Length ];

outputFileName = FileNameJoin[{dataFolder, "compound_75var_a_16_cells_componentwise_vertices.wl"}];

If[FileExistsQ[outputFileName], DeleteFile[outputFileName]];
Save[outputFileName, cells];
Print["Saved polytope component data to " <> outputFileName];

