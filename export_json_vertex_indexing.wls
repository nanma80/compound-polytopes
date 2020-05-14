#!/usr/bin/env wolframscript

dataFolder = "data";
epsilon = 0.0000001;
shape = "compound_840_5_cells";

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

phi = GoldenRatio;
numberLimit = {3, 3};
lookupTable = Flatten[Table[
	{N[i + j * GoldenRatio], {i, j}},
	{i, -numberLimit[[1]], numberLimit[[1]]},
	{j, -numberLimit[[2]], numberLimit[[2]]}
], 1];

dataFileName = FileNameJoin[{dataFolder, shape <> "_componentwise_vertices.wl"}];
outputFileName = FileNameJoin[{dataFolder, shape <> "_vertex_indexing.json"}];

Print["Reading data for " <> shape <> " from " <> dataFileName];

inputData = Get[dataFileName];


Print["Converting " <> IntegerString[Length[inputData]] <> " cells to index based."];
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

overlapping[array1_, array2_] := Length[ Intersection[array1, array2, SameTest -> sameVertex] ] == Length[array1];

components = Map[get5Cells, cellsGoldenField];
Print["Constructed 5-cells. Component count: " <> IntegerString[Length[components]]];
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
