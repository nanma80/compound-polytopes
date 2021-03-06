#!/usr/bin/env wolframscript

numericalValue = True;

epsilon = 1.0*^-5;
dataFolder = "data";


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

shape = "polytope_4_3_3";
Print["shape: " <> shape];
polytopeFileName = FileNameJoin[{dataFolder, shape <> ".wl"}];

polytopeData = Get[polytopeFileName];
polytopeVertices = polytopeData["vertices"];

seed1 = polytopeVertices[[1]];
baseVertices = Select[polytopeVertices, Simplify[normalizedInnerProd[seed1, #]] == 0& ];

compounds = {};

For[seed2Index = 1, seed2Index <= Length[baseVertices], seed2Index++,
    seed2 = baseVertices[[seed2Index]];
    base2Vertices = Select[baseVertices, Simplify[normalizedInnerProd[seed2, #]] == 0& ];
    templateVertices = Join[ {seed1, -seed1, seed2, -seed2}, base2Vertices ];
    AppendTo[compounds, templateVertices];
];


sameSet[array1_, array2_] := Length[ Intersection[array1, array2, SameTest -> sameVertex] ] == Length[array1];

compounds = Union[compounds, SameTest -> sameSet];
(* length == 1, from a common vertex *)

mirrors = {
            {1, 0, 0, 1 },
            {1, 0, 1, 0 },
            {1, 1, 0, 0 },
            {1, 0, 0, 0 }
        };

cells = compounds[[{1}]];
countStart = -1;
mirrorCombinations = Subsets[Range[4], {2}];

overlapping[array1_, array2_] := Length[ Intersection[array1, array2, SameTest -> sameVertex] ] > 0;

For[i = 0, i < 10 && countStart != Length[cells], i++,
    countStart = Length[cells];

    For[j = 1, j <= Length[mirrors], j++,
        newCells = Map[reflect[#, mirrors[[j]] ]&, cells, {2}];
        nonOverlappingNewCells = Complement[newCells, cells, SameTest -> sameSet];
        cells = Join[cells, nonOverlappingNewCells];
        Print[cells//Length];
    ];

    Print["Generated " <> IntegerString[Length[cells]] <> " components so far."];
];

Print[ Intersection[Flatten[cells, 1], polytopeVertices, SameTest -> sameVertex]//Length ];

outputFileName = FileNameJoin[{dataFolder, "compound_2_16_cells_in_8_cell_componentwise_vertices.wl"}];

If[FileExistsQ[outputFileName], DeleteFile[outputFileName]];
Save[outputFileName, cells];
Print["Saved polytope component data to " <> outputFileName];

