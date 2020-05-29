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


shape = "polytope_5_3_3";
Print["shape: " <> shape];
polytopeFileName = FileNameJoin[{dataFolder, shape <> ".wl"}];

polytopeData = Get[polytopeFileName];
polytopeVertices = polytopeData["vertices"];

seed1 = polytopeVertices[[600]];
baseVertices = Select[polytopeVertices, Simplify[normalizedInnerProd[seed1, #]] == 0& ];

compounds = {};

For[seed2Index = 1, seed2Index <= Length[baseVertices], seed2Index++,
    seed2 = baseVertices[[seed2Index]];
    base2Vertices = Select[baseVertices, 
        Simplify[normalizedInnerProd[seed2, #]] == 0
        && #.{1, 0.0001, 0.00000001, 0.000000000001} > 0
        & ];
    orthogonalFrame = {seed1, seed2, base2Vertices[[1]], base2Vertices[[2]]};
    otherVertices = cubeVertices[4].orthogonalFrame/2//Simplify;
    templateVertices = Join[ orthogonalFrame, -orthogonalFrame, otherVertices ];
    AppendTo[compounds, templateVertices];
];


sameSet[array1_, array2_] := Length[ Intersection[array1, array2, SameTest -> sameVertex] ] == Length[array1];

compounds = Union[compounds, SameTest -> sameSet];
(* length == 9, from a common vertex *)
(* 600 * 9 / 24 == 225 *)

phi = (1 + Sqrt[5])/2;

mirrors = {
            {0, 1, 0, 0 },
            {1, phi, -1/phi, 0  },
            {1, -phi, 1/phi, 0  },
            {1, 0, phi, -1/phi  }
        };

Print[compounds[[1]]];

vertex = compounds[[1]][[{1}]];
edge1 = compounds[[1]][[{1, 24}]];
edge2 = compounds[[1]][[{23, 24}]];
face = compounds[[1]][[{1, 23, 24}]];
cell = compounds[[1]][[{1, 2, 21, 22, 23, 24}]];

cells = {cell};
Print[Length[cells[[1]]]];
Print[];

countStart = -1;
mirrorCombinations = Subsets[Range[4], {2}];

overlapping[array1_, array2_] := Length[ Intersection[array1, array2, SameTest -> sameVertex] ] > 0;

For[i = 0, i < 100 && countStart != Length[cells], i++,
    countStart = Length[cells];

    For[j = 1, j <= Length[mirrors], j++,

        newCells = Map[reflect[#, mirrors[[j]] ]&, cells, {2}];
        nonOverlappingNewCells = Complement[newCells, cells, SameTest -> sameSet];
        cells = Join[cells, nonOverlappingNewCells];
        Print[cells//Length];
        (* If[Length[cells] == 200, Break[]]; *)
    ];

    Print["Generated " <> IntegerString[Length[cells]] <> " components so far."];
    (* If[Length[cells] == 200, Break[]]; *)
];

Print[ Intersection[Flatten[cells, 1], polytopeVertices, SameTest -> sameVertex]//Length ];
Print[ToString[Length[cells]/200]];

(* outputFileName = FileNameJoin[{dataFolder, "compound_200_24_cells_componentwise_vertices.wl"}];

If[FileExistsQ[outputFileName], DeleteFile[outputFileName]];
Save[outputFileName, cells];
Print["Saved polytope component data to " <> outputFileName];

 *)

(*  
compound_200_24_cells_vertex_indexing.json
Vertex count: 600
Edge count: 9600
Face count: 19200
Cell count: 4800
Component count: 200

vertex: got 600/600.
edge: 7200/9600 or 2400/9600. Not transitive.
face: over half 13704 > 9600
cell: over half 3868 > 2400. Skip
By h4, an edge (edge1) only generates 7200 edges. If it were fully regular, it should generate 9600. Another edge (edge2), generates the other 2400 edges.

 *)

v = {{(3 + Sqrt[5])/2, (-1 + Sqrt[5])/2, (-1 + Sqrt[5])/2, (-1 + Sqrt[5])/2}, {-1, 1, 1, Sqrt[5]}, {0, 2, -2, 0}, {(3 - Sqrt[5])/2, (-1 - Sqrt[5])/2, (-1 - Sqrt[5])/2, (1 + Sqrt[5])/2}, {(-3 - Sqrt[5])/2, (1 - Sqrt[5])/2, (1 - Sqrt[5])/2, (1 -Sqrt[5])/2}, {1, -1, -1, -Sqrt[5]}, {0, -2, 2, 0}, {(-3 + Sqrt[5])/2, (1 + Sqrt[5])/2, (1 + Sqrt[5])/2, (-1 - Sqrt[5])/2}, {-1, -1, 1, -Sqrt[5]}, {(1 - Sqrt[5])/2, (-3 - Sqrt[5])/2, (1 - Sqrt[5])/2, (1 - Sqrt[5])/2}, {-1, 1, -1, -Sqrt[5]},{(1 - Sqrt[5])/2, (1 - Sqrt[5])/2, (-3 - Sqrt[5])/2, (1 - Sqrt[5])/2}, {-2, 0, 2, 0}, {(-1 - Sqrt[5])/2, (-1 - Sqrt[5])/2, (3 - Sqrt[5])/2, (1 + Sqrt[5])/2}, {-2, 2, 0, 0}, {(-1 - Sqrt[5])/2, (3 - Sqrt[5])/2, (-1 - Sqrt[5])/2, (1 + Sqrt[5])/2}, {(1 + Sqrt[5])/2, (-3 + Sqrt[5])/2, (1 + Sqrt[5])/2, (-1 - Sqrt[5])/2}, {2, -2, 0, 0}, {(1 + Sqrt[5])/2, (1 + Sqrt[5])/2, (-3 + Sqrt[5])/2, (-1 - Sqrt[5])/2}, {2, 0, -2, 0}, {(-1 + Sqrt[5])/2, (-1 + Sqrt[5])/2, (3 + Sqrt[5])/2, (-1 + Sqrt[5])/2}, {1, -1, 1, Sqrt[5]}, {(-1 + Sqrt[5])/2, (3 + Sqrt[5])/2, (-1 + Sqrt[5])/2, (-1 + Sqrt[5])/2}, {1, 1, -1, Sqrt[5]}};
