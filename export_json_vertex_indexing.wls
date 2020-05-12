#!/usr/bin/env wolframscript

dataFolder = "data";
epsilon = 0.0000001;
shape = "compound_120_5_cells";

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

vertices = {};
edges = {};
faces = {};
cells = {};
components = {};

(* prototype *)
edge1 = Subsets[Range[5]-1, {2}];
face1 = Subsets[Range[5]-1, {3}];
cell1 = Table[ Select[Range[Length[face1]], MemberQ[face1[[#]], k] == False& ] - 1, {k, 0, 4}];
component1 = {Range[5] - 1};

For[i = 1, i <= Length[cellsGoldenField], i++,
	verticesBeforeCount = Length[vertices];
	edgesBeforeCount = Length[edges];
	facesBeforeCount = Length[faces];
	cellsBeforeCount = Length[cells];

	vertices = Join[vertices, cellsGoldenField[[i]]];
	edges = Join[edges, edge1 + verticesBeforeCount];
	faces = Join[faces, face1 + verticesBeforeCount];
	cells = Join[cells, cell1 + facesBeforeCount];
	components = Join[components, component1 + cellsBeforeCount];
];

allData = 
	<|
		"vertices" -> vertices,
		"edges" -> edges,
		"faces" -> faces,
		"cells" -> cells,
		"components" -> components
	|>;

If[FileExistsQ[outputFileName], DeleteFile[outputFileName] ];
Export[outputFileName, allData, "RawJSON", Compact -> True];
Print["saved json data to " <> outputFileName];
