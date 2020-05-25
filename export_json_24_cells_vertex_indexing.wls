#!/usr/bin/env wolframscript

dataFolder = "data";
epsilon = 0.0000001;
shape = "compound_seeds_24_cells";

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

cellsGoldenField=Map[toGoldenField[#, filteredTable]&, inputData, {3}];
Print["Mapped input data to golden field"];

signs = Map[#[[{2,3,4}]]&, IntegerDigits[Range[2^3]+2^3-1, 2]]*2-1;
goldenMultiply[x_, y_] := {
    x[[1]] * y[[1]] + x[[2]] * y[[2]],
    x[[1]] * y[[2]] + x[[2]] * y[[1]] + x[[2]] * y[[2]]
};
goldenScale[vector_, scale_] := Map[goldenMultiply[scale, #]&, vector];
goldenDot[vector1_, vector2_] := Total[Table[goldenMultiply[vector1[[k]], vector2[[k]]], {k, Length[vector1]}]];
toReal[goldenVector_] := goldenVector[[1]] + phi * goldenVector[[2]]//Simplify;

get16Cells[vertices_] := Module[
	{positiveVertices, signIndex},
	cells = {};
	positiveVertices = Select[vertices, Map[toReal, #].{1, 0.0001, 0.00000001, 0.000000000001} > 0&];

	For[signIndex = 1, signIndex <= 2^4, signIndex++,
		verticesInCell = Table[positiveVertices[[k]] * signs[[signIndex]][[k]], {k, 4}];
		facesInCell = Subsets[verticesInCell, {3}];
		AppendTo[cells, facesInCell];
	];
	cells
];

getFacesInOctahedralCell[octahedralVertices_] := Module[
	{positiveVertices, signIndex, face},
	faces = {};
	center = Mean[octahedralVertices];
	relativeVertices = Map[# - center &, octahedralVertices];
	positiveVertices = Select[relativeVertices, Map[toReal, #].{1, 0.0001, 0.00000001, 0.000000000001} > 0&];

	For[signIndex = 1, signIndex <= 2^3, signIndex++,
		face = Table[positiveVertices[[k]] * signs[[signIndex]][[k]] + center, {k, 3}];
		AppendTo[faces, face];
	];
	faces
];


get24Cells[vertices_] := Module[
	{positiveVertices, signIndex},
	cells = {};
	seed1 = vertices[[1]];
	baseVertices = Select[vertices, Simplify[goldenDot[seed1, #]] == {0, 0}& ];
	seed2 = baseVertices[[1]];
	base2Vertices = Select[baseVertices, 
	    Simplify[goldenDot[seed2, #]] == {0, 0}
	    && Map[toReal, #].{1, 0.0001, 0.00000001, 0.000000000001} > 0
	    & ];
	orthogonalFrame = {seed1, seed2, base2Vertices[[1]], base2Vertices[[2]]};

	connectingDimensions = Subsets[Range[4], {2}];
	For[cdIndex = 1, cdIndex <= Length[connectingDimensions], cdIndex++,
		connectingDimension = connectingDimensions[[cdIndex]];
		otherDimensions = Complement[Range[4], connectingDimension];
		For[cd1Sign = -1, cd1Sign <= 1, cd1Sign += 2,
			For[cd2Sign = -1, cd2Sign <= 1, cd2Sign += 2,
				coefs = {};
				coef1 = {0, 0, 0, 0};
				coef1[[connectingDimension[[1]] ]] = cd1Sign;
				AppendTo[coefs, coef1];

				coef2 = {0, 0, 0, 0};
				coef2[[connectingDimension[[2]] ]] = cd2Sign;
				AppendTo[coefs, coef2];

				For[od1Sign = -1, od1Sign <= 1, od1Sign += 2,
					For[od2Sign = -1, od2Sign <= 1, od2Sign += 2,
						coefn = {1, 1, 1, 1};
						coefn[[connectingDimension[[1]] ]] = cd1Sign;
						coefn[[connectingDimension[[2]] ]] = cd2Sign;
						coefn[[otherDimensions[[1]] ]] = od1Sign;
						coefn[[otherDimensions[[2]] ]] = od2Sign;
						coefn = coefn/2;
						AppendTo[coefs, coefn];
					];
				];
				verticesInCell = coefs.orthogonalFrame;
				facesInCell = getFacesInOctahedralCell[verticesInCell];
				AppendTo[cells, facesInCell];
			];
		];
	];
	cells
];

(* vertices = cellsGoldenField[[1]]; *)



octahedralVertices = verticesInCell;




overlapping[array1_, array2_] := Length[ Intersection[array1, array2, SameTest -> sameVertex] ] == Length[array1];

components = Map[get24Cells, cellsGoldenField];
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

seeds:
Vertex count: 120
Edge count: 744
Face count: 864
Cell count: 216
Component count: 9

The right quaternion to use is [1+phi, -1+phi, 1-phi, -1+phi]


crop command:
mogrify -crop 1164x1164+493+0 *.png

 *)