#!/usr/bin/env wolframscript

dataFolder = "data";
epsilon = 0.0000001;

projectTo3d[vertex4d_] := {vertex4d[[1]], vertex4d[[2]], vertex4d[[3]]};
reflect[point_, mirror_] := FullSimplify[point - 2 * (point.mirror)/(mirror.mirror)*mirror];
doubleReflect[point_, mirror1_, mirror2_] := reflect[reflect[point, mirror1], mirror2];
getEdgesFromFace[face_] := Table[{face[[i+1]], face[[Mod[i+1, Length[face]] + 1]]}, {i, 0, Length[face] - 1}];
explodedFace[face_, explosionFactor_] := Map[(# + Mean[face] * explosionFactor)&, face];
sameVertex[point1_, point2_] := point1 == point2;
sameCenter[edge1_, edge2_] := Total[edge1] == Total[edge2];
sameCellCenter[cell1_, cell2_] := getCellCenter[cell1] == getCellCenter[cell2];
getCenter[edge_] := Total[edge];
getCellCenter[cell_] := Simplify[Total[Flatten[cell, 1]]];

toGoldenField[value_, lookupTable_] := 
	Intersection[lookupTable, {{N[value]}}, SameTest -> (Abs[(#1[[1]] - #2[[1]])] < epsilon&)][[1,2]];

getMap[list_] := Association[Table[{list[[k]] -> k}, {k, 1, Length[list]}]];

phi = (1+ Sqrt[5])/2;
numberLimit = {3, 3};
lookupTable = Flatten[Table[
	{N[i + j * GoldenRatio], {i, j}},
	{i, -numberLimit[[1]], numberLimit[[1]]},
	{j, -numberLimit[[2]], numberLimit[[2]]}
], 1];

shapeVar = "compound_120var_a_5_cells";

dataFileVarName = FileNameJoin[{dataFolder, shapeVar <> "_componentwise_vertices.wl"}];

Print["Reading data for " <> shapeVar <> " from " <> dataFileVarName];

inputDataVar = Get[dataFileVarName];

shapeOriginal = "compound_120_5_cells";

dataFileOriginalName = FileNameJoin[{dataFolder, shapeOriginal <> "_componentwise_vertices.wl"}];

Print["Reading data for " <> shapeOriginal <> " from " <> dataFileOriginalName];

inputDataOriginal = Get[dataFileOriginalName];

leftMatrix[v_] := {
    {v[[1]], -v[[2]], -v[[3]], -v[[4]]},
    {v[[2]], v[[1]], -v[[4]], v[[3]]},
    {v[[3]], v[[4]], v[[1]], -v[[2]]},
    {v[[4]], -v[[3]], v[[2]], v[[1]]}
};

rightMatrix[v_] := {
    {v[[1]], -v[[2]], -v[[3]], -v[[4]]},
    {v[[2]], v[[1]], v[[4]], -v[[3]]},
    {v[[3]], -v[[4]], v[[1]], v[[2]]},
    {v[[4]], v[[3]], -v[[2]], v[[1]]}
};

vector = {0, 1, 0, 0};
leftMatrix[vector]
Eigenvalues[leftMatrix[vector]]
Eigenvectors[leftMatrix[vector]]

(* sample value of 600 cell vertex, from 335 data*)
(* {phi, 1, 0, 1/phi} *)

(* cell center from 533, compatible *)
(* {1, phi, 1/phi, 0} *)

