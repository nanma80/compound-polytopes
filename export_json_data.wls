#!/usr/bin/env wolframscript

shapes = {
    (* "polytope_3_3_4" *)
    (* "polytope_3_4_3_alt" *)
    "polytope_3_3_5"
};

dataFolder = "data";
epsilon = 1.0*^-5;

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

getIndexBasedPolytopeData[cells_] := Module[
	{
		cellsGoldenField, vertices, verticesMap, faces, edgesPerFace, edges, 
		edgeIndices, edgeCenters, edgeCentersMap, edgeCenterPerFace, faceIndices,
		faceCenters, faceCentersMap, faceCentersPerCell, cellIndices
	},
	
	Print["Converting " <> IntegerString[Length[cells]] <> " cells to index based."];
	cellsGoldenField=Map[toGoldenField[#, lookupTable]&, cells, {4}];
	Print["Mapped values to the golden ratio field"];

	vertices = Union[Flatten[cellsGoldenField,2]];
	verticesMap = getMap[vertices];

	faces = Union[Flatten[cellsGoldenField,1], SameTest -> sameCenter];
	edgesPerFace = Map[getEdgesFromFace, faces];
	edges = Union[Flatten[edgesPerFace,1], SameTest -> sameCenter];
	edgeIndices = Map[Sort[Map[verticesMap, #]]&, edges];
	edgeCenters = Map[Total, edges];
	edgeCentersMap = getMap[edgeCenters];
	edgeCenterPerFace = Map[getCenter, edgesPerFace, {2}];

	faceIndices = Map[Sort[Map[edgeCentersMap, #]]&, edgeCenterPerFace];
	faceCenters = Map[Total, faces];
	faceCentersMap = getMap[faceCenters];

	faceCentersPerCell = Map[Total[#, {2}]&, cellsGoldenField];
	cellIndices = Map[Sort[Map[faceCentersMap, #]]&, faceCentersPerCell];

	<|
		"vertices" -> vertices,
		"edges" -> edgeIndices - 1,
		"faces" -> faceIndices - 1,
		"cells" -> cellIndices - 1
	|>
];

phi = GoldenRatio;
numberLimit = {12, 6};
lookupTable = Flatten[Table[
	{N[i + j * GoldenRatio], {i, j}},
	{i, -numberLimit[[1]], numberLimit[[1]]},
	{j, -numberLimit[[2]], numberLimit[[2]]}
], 1];

For[i=1,i<=Length[shapes],i++,
	shape = shapes[[i]];
	outputFileName = shape <> ".json";
	dataFileName = FileNameJoin[{dataFolder, shape <> ".wl"}];
	outputFileName = FileNameJoin[{dataFolder, outputFileName}];

	Print["Reading data for " <> shape <> " from " <> dataFileName];

	inputData = Get[dataFileName];
	(* Print[Keys[inputData]]; *)

(* 	{cells, faces, edges, vertices, edgeComponentVertices1, edgeComponentEdges1, edgeComponentFaces1, edgeComponentCells1, edgeComponentVertices2, edgeComponentEdges2, edgeComponentFaces2, edgeComponentCells2, faceComponentVertices1, faceComponentEdges1, faceComponentFaces1, faceComponentCells1, faceComponentVertices2, faceComponentEdges2, faceComponentFaces2, faceComponentCells2, cellComponentVertices1, cellComponentEdges1, cellComponentFaces1, cellComponentCells1, cellComponentVertices2, cellComponentEdges2, cellComponentFaces2, cellComponentCells2} *)

	scalingFactor = 6;
	allData = <|
		"original" -> getIndexBasedPolytopeData[inputData["cells"] * scalingFactor]
	|>;

	If[MemberQ[Keys[inputData], "edgeComponentCells1"], 
		Print["edgeComponentCells1"];
		allData[["edgeCenterCompound1"]] = Map[getIndexBasedPolytopeData[#]& , inputData["edgeComponentCells1"] * scalingFactor] ;
	];

	If[MemberQ[Keys[inputData], "edgeComponentCells2"], 
		Print["edgeComponentCells2"];
		allData[["edgeCenterCompound2"]] = Map[getIndexBasedPolytopeData[#]& , inputData["edgeComponentCells2"] * scalingFactor] ;
	];

	If[MemberQ[Keys[inputData], "faceComponentCells1"], 
		Print["faceComponentCells1"];
		allData[["faceCenterCompound1"]] = Map[getIndexBasedPolytopeData[#]& , inputData["faceComponentCells1"] * scalingFactor] ;
	];

	If[MemberQ[Keys[inputData], "faceComponentCells2"], 
		Print["faceComponentCells2"];
		allData[["faceCenterCompound2"]] = Map[getIndexBasedPolytopeData[#]& , inputData["faceComponentCells2"] * scalingFactor] ;
	];

(* 	If[MemberQ[Keys[inputData], "cellComponentCells1"], 
		allData[["cellCenterCompound1"]] = Map[getIndexBasedPolytopeData[#]& , inputData["cellComponentCells1"] * scalingFactor] ;
	];

	If[MemberQ[Keys[inputData], "cellComponentCells2"], 
		allData[["cellCenterCompound2"]] = Map[getIndexBasedPolytopeData[#]& , inputData["cellComponentCells2"] * scalingFactor] ;
	];
 *)
	If[FileExistsQ[outputFileName], DeleteFile[outputFileName] ];
	Export[outputFileName, allData, "RawJSON", Compact -> True];
	Print["saved json data to " <> outputFileName];
];
