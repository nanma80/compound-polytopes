#!/usr/bin/env wolframscript

dataFolder = "data";
epsilon = 0.0000001;
shape = "compound_15_8_cells_in_600_cell";
swapDimensions = False;

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

phi = (1+Sqrt[5])/2;
numberLimit = {3, 3};
lookupTable = Flatten[Table[
	{N[i + j * GoldenRatio], {i, j}},
	{i, -numberLimit[[1]], numberLimit[[1]]},
	{j, -numberLimit[[2]], numberLimit[[2]]}
], 1];

dataFileName = FileNameJoin[{dataFolder, shape <> "_componentwise_vertices.wl"}];
outputFileName = FileNameJoin[{dataFolder, "json", shape <> If[swapDimensions, "_swap", ""] <> "_vertex_indexing.json"}];

Print["Reading data for " <> shape <> " from " <> dataFileName];

inputData = Get[dataFileName];


Print["Converting " <> IntegerString[Length[inputData]] <> " components to index based."];
possibleValues = Union[Flatten[inputData], SameTest -> (Abs[N[#1]-N[#2]]<epsilon&)];
filteredTable = Intersection[lookupTable, Table[{N[v]}, {v, possibleValues}], SameTest -> (Abs[(#1[[1]] - #2[[1]])] < epsilon&)];

cellsGoldenField=Map[toGoldenField[#, filteredTable]&, inputData, {3}];
Print["Mapped input data to golden field"];
get5Cells[vertices_] := Module[
	{},
	cells = {};
	For[i=1, i<=Length[vertices], i++,
		verticesInCell = vertices[[Complement[Range[5], {i}]]];
		facesInCell = Subsets[verticesInCell, {3}];
		AppendTo[cells, facesInCell];
	];
	cells
];

goldenMultiply[x_, y_] := {
    x[[1]] * y[[1]] + x[[2]] * y[[2]],
    x[[1]] * y[[2]] + x[[2]] * y[[1]] + x[[2]] * y[[2]]
};
goldenScale[vector_, scale_] := Map[goldenMultiply[scale, #]&, vector];
goldenDot[vector1_, vector2_] := Total[Table[goldenMultiply[vector1[[k]], vector2[[k]]], {k, Length[vector1]}]];
toReal[goldenVector_] := goldenVector[[1]] + phi * goldenVector[[2]]//Simplify;

getFaceIndex[activeTwo_, others_] := Map[Join[#, others]&, {{}, {activeTwo[[1]]}, activeTwo, {activeTwo[[2]]}}];

get8Cells[vertices_] := Module[
	{positiveVertices, signIndex},
	cells = {};
	seed1 = vertices[[1]];

	offsets = Map[# - seed1&, vertices];
	minDistance = Sort[Map[toReal[goldenDot[#, #]]&, offsets]][[2]];
	neighborOffsets = Select[offsets, toReal[goldenDot[#, #]] == minDistance&];

	getFace[faceIndex_] := Map[Total[Join[neighborOffsets[[#]], {seed1}]]&, faceIndex ]//Simplify;

	cells = {};
	For[i = 1, i <= 4, i++,
		activeThree = Complement[Range[4], {i}];
		For[activeI = 0, activeI <= 1, activeI++,
			appendI = If[activeI == 1, {i}, {}];
			facesInCell = {};
			faceIndices = {
				getFaceIndex[activeThree[[{1, 2}]], Join[{}, appendI]],
				getFaceIndex[activeThree[[{1, 2}]], Join[{activeThree[[3]]}, appendI]],
				getFaceIndex[activeThree[[{1, 3}]], Join[{}, appendI]],
				getFaceIndex[activeThree[[{1, 3}]], Join[{activeThree[[2]]}, appendI]],
				getFaceIndex[activeThree[[{2, 3}]], Join[{}, appendI]],
				getFaceIndex[activeThree[[{2, 3}]], Join[{activeThree[[1]]}, appendI]]
			};
			facesInCell = Map[getFace, faceIndices];
			AppendTo[cells, facesInCell];
		];
	];
	cells
];


overlapping[array1_, array2_] := Length[ Intersection[array1, array2, SameTest -> sameVertex] ] == Length[array1];
swap[v_] := v[[{1, 2, 4, 3}]];
If[swapDimensions, cellsGoldenField =  Map[swap, cellsGoldenField, {2}]];

components = Map[get8Cells, cellsGoldenField];
Print["Constructed 8-cells. Component count: " <> IntegerString[Length[components]]];
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

MinRotation[list_] := 
	SortBy[{
			RotateLeft[list, Position[list, Min[list]][[1,1]]-1],
			RotateLeft[Reverse[list], Position[Reverse[list], Min[Reverse[list]]][[1,1]]-1]
	}, #[[2]]&][[1]];

faceByVertexIndices = Map[MinRotation[Map[verticesMap, #]]&, faces];
facesMap = getMap[faceByVertexIndices];
Print["Mapped faces"];

cells = Union[Flatten[components, 1], SameTest -> overlapping];
Print["Flattened cells. Cell count: " <> IntegerString[Length[cells]]];

cellsByVertexIndices = Map[MinRotation, Map[verticesMap, cells, {3}], {2}];
cellsByFaceIndices = Map[Sort, Map[facesMap, cellsByVertexIndices, {2}]];
cellsMap = getMap[cellsByFaceIndices];
Print["Mapped cells"];

componentsByVertexIndices = Map[MinRotation, Map[verticesMap, components, {4}], {3}];
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

quaternion to import seeds: [1+phi, -1+phi, 1-phi, -1+phi]

compound_seeds_8_cells_vertex_indexing.json
Vertex count: 120
Edge count: 496
Face count: 432
Cell count: 144
Component count: 18

compound_75_8_cells_vertex_indexing.json
Vertex count: 600
Edge count: 2400
Face count: 1800
Cell count: 600
Component count: 75

compound_600_8_cells_vertex_indexing.json
Vertex count: 600
Edge count: 9600
Face count: 14400
Cell count: 4800
Component count: 600

compound_300_8_cells_vertex_indexing.json
Vertex count: 600
Edge count: 6000
Face count: 7200
Cell count: 2400
Component count: 300

compound_75var_a_8_cells_vertex_indexing.json
Vertex count: 600
Edge count: 2400
Face count: 1800
Cell count: 600
Component count: 75

compound_675_8_cells_vertex_indexing.json
Vertex count: 600
Edge count: 9600
Face count: 16200
Cell count: 5400
Component count: 675

crop command:
mogrify -crop 1164x1164+493+0 *.png

 *)