#!/usr/bin/env wolframscript

dataFolder = "data";
epsilon = 0.0000001;
shape = "compound_10_600_cells_in_120_cell";

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
normalizedInnerProd[u_, v_] := u.v/Norm[u]/Norm[v];
toGoldenField[value_, lookupTable_] := 
	Intersection[lookupTable, {{N[value]}}, SameTest -> (Abs[(#1[[1]] - #2[[1]])] < epsilon&)][[1,2]];
cubeVertices[dimension_] := 2 * Map[IntegerDigits[#, 2][[2;;dimension+1]]&, Table[k, {k, 2^dimension, 2 * 2^dimension - 1}]] - 1;
getMap[list_] := Association[Table[{list[[k]] -> k}, {k, 1, Length[list]}]];

phi = (1+Sqrt[5])/2;
numberLimit = {3, 3};
lookupTable = Flatten[Table[
	{N[i + j * GoldenRatio], {i, j}},
	{i, -numberLimit[[1]], numberLimit[[1]]},
	{j, -numberLimit[[2]], numberLimit[[2]]}
], 1];

dataFileName = FileNameJoin[{dataFolder, shape <> "_componentwise_vertices.wl"}];
outputFileName = FileNameJoin[{dataFolder, "json", shape <> "_vertex_indexing.json"}];

Print["Reading data for " <> shape <> " from " <> dataFileName];

inputData = Get[dataFileName];


Print["Converting " <> IntegerString[Length[inputData]] <> " components to index based."];
possibleValues = Union[Flatten[inputData], SameTest -> (Abs[N[#1]-N[#2]]<epsilon&)];
filteredTable = Intersection[lookupTable, Table[{N[v]}, {v, possibleValues}], SameTest -> (Abs[(#1[[1]] - #2[[1]])] < epsilon&)];

Print["Mapping input data to golden field"];
cellsGoldenField=Map[toGoldenField[#, filteredTable]&, inputData, {5}];
Print["Mapped input data to golden field"];

signs = Map[#[[{2,3,4}]]&, IntegerDigits[Range[2^3]+2^3-1, 2]]*2-1;
goldenMultiply[x_, y_] := {
    x[[1]] * y[[1]] + x[[2]] * y[[2]],
    x[[1]] * y[[2]] + x[[2]] * y[[1]] + x[[2]] * y[[2]]
};
goldenScale[vector_, scale_] := Map[goldenMultiply[scale, #]&, vector];
goldenDot[vector1_, vector2_] := Total[Table[goldenMultiply[vector1[[k]], vector2[[k]]], {k, Length[vector1]}]];
toReal[goldenVector_] := goldenVector[[1]] + phi * goldenVector[[2]]//Simplify;


overlapping[array1_, array2_] := Length[ Intersection[array1, array2, SameTest -> sameVertex] ] == Length[array1];

components = cellsGoldenField;

Print["Constructed 16-cells. Component count: " <> IntegerString[Length[components]]];
vertices = Union[Flatten[components, 3]];
verticesMap = getMap[vertices];
Print["Mapped vertices. Vertex count: " <> IntegerString[Length[vertices]]];
faces = Union[Flatten[components, 2], SameTest -> overlapping];
Print["Flattened faces. Face count: " <> IntegerString[Length[faces]]];
edgesPerFace = Map[getEdgesFromFace, faces];
Print["Created edgesPerFace"];
edges = Union[Flatten[edgesPerFace, 1], SameTest -> overlapping];
Print["Flattened edges. Edge count: " <> IntegerString[Length[edges]]];
edgeIndices = Map[Sort[Map[verticesMap, #]]&, edges];
Print["Mapped edges"];

faceByVertexIndices = Map[Sort[Map[verticesMap, #]]&, faces];
facesMap = getMap[faceByVertexIndices];
Print["Mapped faces"];

cells = Union[Flatten[components, 1], SameTest -> overlapping];
Print["Flattened cells. Cell count: " <> IntegerString[Length[cells]]];

cellsByVertexIndices = Map[Sort, Map[verticesMap, cells, {3}], {2}];
cellsByFaceIndices = Map[Sort, Map[facesMap, cellsByVertexIndices, {2}]];
cellsMap = getMap[cellsByFaceIndices];
Print["Mapped cells"];

componentsByVertexIndices = Map[Sort, Map[verticesMap, components, {4}], {3}];
componentsByFaceIndices = Map[Sort, Map[facesMap, componentsByVertexIndices, {3}], {2}];
componentsByCellIndices = Map[Sort, Map[cellsMap, componentsByFaceIndices, {2}], {1}];
Print["Mapped components"];

allData = 
	<|
		"vertices" -> vertices,
		"edges" -> edgeIndices - 1,
		"faces" -> faceByVertexIndices - 1,
		"cells" -> cellsByFaceIndices - 1,
		"components" -> componentsByCellIndices - 1
	|>;

If[FileExistsQ[outputFileName], DeleteFile[outputFileName] ];
Export[outputFileName, allData, "RawJSON", Compact -> True];
Print["saved json data to " <> outputFileName];

Print["Vertex count: " <> IntegerString[Length[vertices]]];
Print["Edge count: " <> IntegerString[Length[edges]]];
Print["Face count: " <> IntegerString[Length[faces]]];
Print["Cell count: " <> IntegerString[Length[cells]]];
Print["Component count: " <> IntegerString[Length[components]]];

(* 

compound_5_600_cells_in_120_cell_vertex_indexing.json
Vertex count: 600
Edge count: 3600
Face count: 6000
Cell count: 3000
Component count: 5

compound_10_600_cells_in_120_cell_vertex_indexing.json
Vertex count: 600
Edge count: 7200
Face count: 12000
Cell count: 6000
Component count: 10

crop command:
mogrify -crop 1164x1164+493+0 *.png

 *)