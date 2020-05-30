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

shape600 = "polytope_3_3_5";
Print["shape: " <> shape600];
polytopeFileName600 = FileNameJoin[{dataFolder, shape600 <> ".wl"}];

polytopeData600 = Get[polytopeFileName600];
polytopeVertices = polytopeData600["vertices"];


shape120 = "polytope_5_3_3";
Print["shape: " <> shape120];
polytopeFileName120 = FileNameJoin[{dataFolder, shape120 <> ".wl"}];

polytopeData120 = Get[polytopeFileName120];

(* 600 cell component vertex *)
factor = (polytopeData120["vertices"][[1]]//Norm)/(polytopeData600["cellComponentVertices1"][[1, 1]]//Norm//Simplify)//FullSimplify;
(* 6 - 2 Sqrt[5] *)

(* 
cellComponentVertices1, cellComponentEdges1,cellComponentFaces1, cellComponentCells1,
cellComponentVertices2, cellComponentEdges2, cellComponentFaces2,cellComponentCells2

cellComponentVertices1: 5 * 120. each element is a length 4 list
cellComponentEdges1: 5 * 720. each element is a length 2 list. each element has length 4
cellComponentFaces1: 5 * 1200. each face is has 3 length 4 lists.
cellComponentCells1: 5 * 600 * 4 * 3. each element is length 4

*)

verticesInAllCells = Union[Flatten[polytopeData600["cellComponentCells1"], 3]];

Print[Intersection[polytopeData120["vertices"], verticesInAllCells*factor, SameTest -> sameVertex]//Length];
(* 600 *)

outputData = polytopeData600["cellComponentCells1"] * factor;

outputFileName = FileNameJoin[{dataFolder, "compound_5_600_cells_in_120_cell_componentwise_vertices.wl"}];

If[FileExistsQ[outputFileName], DeleteFile[outputFileName]];
Save[outputFileName, outputData];
Print["Saved polytope component data to " <> outputFileName];

