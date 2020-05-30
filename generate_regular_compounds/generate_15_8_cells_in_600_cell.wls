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

shape = "polytope_3_3_5";
Print["Loading data from " <> shape];
polytopeFileName = FileNameJoin[{dataFolder, shape <> ".wl"}];

polytopeData = Get[polytopeFileName];
polytopeVertices = polytopeData["vertices"];

seed1 = polytopeVertices[[1]];
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
(* length == 10 *)
<< Quaternions`

icosianLists = polytopeVertices/2;

quaternionRotation[leftList_, rightList_] := List @@ ( (Quaternion @@ leftList) ** (Quaternion @@ rightList) );

cells = 
    Table[
        quaternionRotation[compounds[[1, j]], icosianLists[[i]]],
        {i, 1, Length[icosianLists]},
        {j, 1, Length[compounds[[1]]]}
    ];

cells = Union[cells, SameTest -> sameSet];

Print[Length[cells]];
Print[ Intersection[Flatten[cells, 1], polytopeVertices, SameTest -> sameVertex]//Length ];

outputFileName = FileNameJoin[{dataFolder, "compound_15_8_cells_in_600_cell_componentwise_vertices.wl"}];

If[FileExistsQ[outputFileName], DeleteFile[outputFileName]];
Save[outputFileName, cells];
Print["Saved polytope component data to " <> outputFileName];

