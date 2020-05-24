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
SameMatrix[m1_, m2_] := (Max[Max[N[Abs[m1 - m2]]]] < epsilon);
sameCenter[edge1_, edge2_] := SameMatrix[Total[edge1], Total[edge2]];
SameTotal[edge1_, edge2_] := SameMatrix[Total[edge1], Total[edge2]];
getQString[q_, separator_] := If[IntegerQ[q], IntegerString[q], IntegerString[Numerator[q]] <> separator <> IntegerString[Denominator[q]]];
getFileName[symbol_] := "complex_polytope_" <> IntegerString[symbol[[1]]] <> "_" <> getQString[symbol[[2]], "over"] <> "_" <> IntegerString[symbol[[3]]];
getShapeName[symbol_] := IntegerString[symbol[[1]]] <> "{" <> getQString[symbol[[2]], "/"] <> "}" <> IntegerString[symbol[[3]]];
maxNorm[direction_, realVertices_] := Max[Map[(#.direction)^2/direction.direction&, realVertices]];
getReal[complexVector_] := Flatten[ Map[{Re[#], Im[#]}&, complexVector] ];

shape = "polytope_5_3_3";
Print["Loading data from " <> shape];
polytopeFileName = FileNameJoin[{dataFolder, shape <> ".wl"}];

polytopeData = Get[polytopeFileName];
polytopeVertices = polytopeData["vertices"];

seed1 = polytopeVertices[[600]];
baseVertices = Select[polytopeVertices, Simplify[normalizedInnerProd[seed1, #]] == 1/2 & ];

orthogonalIndices =
    Select[Subsets[Range[Length[baseVertices]], {4}], 
        Simplify[baseVertices[[#[[1]]]].baseVertices[[#[[2]]]]] == 0
        && Simplify[baseVertices[[#[[1]]]].baseVertices[[#[[3]]]]] == 0
        && Simplify[baseVertices[[#[[1]]]].baseVertices[[#[[4]]]]] == 0
        && Simplify[baseVertices[[#[[2]]]].baseVertices[[#[[3]]]]] == 0
        && Simplify[baseVertices[[#[[2]]]].baseVertices[[#[[4]]]]] == 0
        && Simplify[baseVertices[[#[[3]]]].baseVertices[[#[[4]]]]] == 0 &
    ];

compounds = {};

For[ i = 1, i <= Length[orthogonalIndices], i++,
    neighbors = baseVertices[[orthogonalIndices[[i]]]];
    vectors = Map[#-seed1 &, neighbors]//Simplify;
    templateVertices = Map[Total[Join[vectors[[#]], {seed1}]]&, Subsets[Range[Length[vectors]]] ]//Simplify;
    AppendTo[compounds, templateVertices];
];

sameSet[array1_, array2_] := Length[ Intersection[array1, array2, SameTest -> sameVertex] ] == Length[array1];

compounds = Union[compounds, SameTest -> sameSet];
(* length == 18 *)

(* outputFileName = FileNameJoin[{dataFolder, "compound_seeds_8_cells_componentwise_vertices.wl"}];

If[FileExistsQ[outputFileName], DeleteFile[outputFileName]];
Save[outputFileName, compounds];
Print["Saved polytope component data to " <> outputFileName];
 *)

(* verified until here *)

phi = (1 + Sqrt[5])/2;

mirrors = {
            {0, 1, 0, 0 },
            {1, phi, -1/phi, 0  },
            {1, -phi, 1/phi, 0  },
            {1, 0, phi, -1/phi  }
        };

(* 1: original 75 *)
cells = {compounds[[2]]};
countStart = -1;
mirrorCombinations = Subsets[Range[4], {2}];

overlapping[array1_, array2_] := Length[ Intersection[array1, array2, SameTest -> sameVertex] ] > 0;

For[i = 0, i < 100 && countStart != Length[cells], i++,
    countStart = Length[cells];

    (* For[j = 1, j <= Length[mirrorCombinations], j++,
        mirrorCombination = mirrorCombinations[[j]];
        innerCount = Length[cells];
        For[k = 1, k <= 12, k++,
            newCells = Map[doubleReflect[#, mirrors[[mirrorCombination[[1]]]], mirrors[[mirrorCombination[[2]] ]] ]&, cells, {2}];
            nonOverlappingNewCells = Complement[newCells, cells, SameTest -> overlapping];
            cells = Join[cells, nonOverlappingNewCells];

            Print[cells//Length];

            newCells = Map[doubleReflect[#, mirrors[[mirrorCombination[[2]]]], mirrors[[mirrorCombination[[1]] ]] ]&, cells, {2}];
            nonOverlappingNewCells = Complement[newCells, cells, SameTest -> overlapping];
            cells = Join[cells, nonOverlappingNewCells];

            Print[cells//Length];

            If[Length[cells] == innerCount, Break[]];
            innerCount = Length[cells];
            If[Length[cells] == 600, Break[]];
        ];
        If[Length[cells] == 600, Break[]];
    ]; *)

    For[j = 1, j <= Length[mirrors], j++,
        newCells = Map[reflect[#, mirrors[[j]] ]&, cells, {2}];
        nonOverlappingNewCells = Complement[newCells, cells, SameTest -> sameSet];
        cells = Join[cells, nonOverlappingNewCells];
        Print[cells//Length];
    ];

    Print["Generated " <> IntegerString[Length[cells]] <> " components so far."];
    If[Length[cells] == 600, Break[]];
];

Print[ Intersection[Flatten[cells, 1], polytopeVertices, SameTest -> sameVertex]//Length ];

outputFileName = FileNameJoin[{dataFolder, "compound_600_8_cells_componentwise_vertices.wl"}];

If[FileExistsQ[outputFileName], DeleteFile[outputFileName]];
Save[outputFileName, cells];
Print["Saved polytope component data to " <> outputFileName];

