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
snap[value_, knownValues_] := Module[
    {epsilon, snapValue, knownIndex, knownValue},
    epsilon = 0.00001;
    For[knownIndex = 1, knownIndex <= Length[knownValues], knownIndex++,
        knownValue = knownValues[[knownIndex]];
        For[offset = -4, offset <= 4, offset ++,
            If[Abs[value - knownValue - offset] < epsilon, Return[knownValue + offset, Module] ];
            If[Abs[value + knownValue + offset] < epsilon, Return[-knownValue - offset, Module] ]
        ];
    ];
    value
];

nVector[r1_, p1_] := Module[{n1Matrix, n11, n12, epsilon}, 
    epsilon = 0.000001;
    n1Matrix = (r1 - IdentityMatrix[2])/(Exp[2 Pi I / p1] - 1);
    n11 = Sqrt[Re[n1Matrix[[1]][[1]]]];
    n12 = n1Matrix[[1]][[2]]/n11;
    If[Abs[Im[n12]] < epsilon, n12 = Re[n12]];
    {n11, n12}
];

getParts[edge1_, numParts_] := Module[
    {alpha, betaStep},
    alpha = ArcCos[edge1[[1]].edge1[[2]]/Norm[edge1[[1]]]^2];
    betaStep = alpha/numParts;
    Table[(Sin[alpha - k * betaStep] edge1[[1]] + Sin[k * betaStep] edge1[[2]])/Sin[alpha], {k, 1, numParts - 1}]
];

doubleRotation[theta_] := ArrayFlatten[{
    {RotationMatrix[theta], 0},
    {0, RotationMatrix[theta]}
}];

isVertexFirst[vertices_] := Module[
    {center},
    center = Select[vertices, Abs[N[#[[1]]]] + Abs[N[#[[2]]]] + Abs[N[#[[3]]]] < epsilon&];
    Length[center] > 0
];

isEdgeFaceFirst[edges_] := Module[
    {centers},
    centers = Mean/@edges;
    Length[Select[centers, Abs[N[#[[1]]]] + Abs[N[#[[2]]]] + Abs[N[#[[3]]]] < epsilon&]]>0
];

isCellFirst[cells_] := Module[
    {centers},
    centers = Map[Mean[Mean[#]]&, cells];
    Length[Select[centers, Abs[N[#[[1]]]] + Abs[N[#[[2]]]] + Abs[N[#[[3]]]] < epsilon&]]>0
];

getCenterType[vertices_, edges_, faces_, cells_] := 
    If[isVertexFirst[vertices],
        "Vertex",
        If[isEdgeFaceFirst[edges],
            "Edge",
            If[isEdgeFaceFirst[faces],
                "Face",
                If[isCellFirst[cells],
                    "Cell",
                    "None"
                ]
            ]
        ]
    ];

getInnerProds[points_] := Module[
    {fInnerProds},
    fInnerProds = {};
    For[fIndex = 2, fIndex <= Length[points], fIndex++,
        AppendTo[fInnerProds, normalizedInnerProd[points[[1]], points[[fIndex]]]];
    ];
    Union[fInnerProds, SameTest -> (Abs[N[#1] - N[#2]] < epsilon&)]
];

shape = "polytope_5_3_3";
Print["shape: " <> shape];
polytopeFileName = FileNameJoin[{dataFolder, shape <> ".wl"}];

polytopeData = Get[polytopeFileName];
polytopeVertices = polytopeData["vertices"];

seed1 = polytopeVertices[[600]];
baseVertices = Select[polytopeVertices, Simplify[normalizedInnerProd[seed1, #] + 1/4] == 0& ];

compounds = {};

For[seed2Index = 1, seed2Index <= Length[baseVertices], seed2Index++,
    seed2 = baseVertices[[seed2Index]];
    base2Vertices = Select[baseVertices, Simplify[normalizedInnerProd[seed2, #] + 1/4] == 0& ];
    templateVertices = Join[ {seed1, seed2}, base2Vertices ];
    AppendTo[compounds, templateVertices];
];


same5Cell[array1_, array2_] := Length[ Intersection[array1, array2, SameTest -> sameVertex] ] == 5;

compounds = Union[compounds, SameTest -> same5Cell];

compounds = compounds[[{2}]];

Print[Length[compounds]];

For[i = 1, i<=Length[compounds], i++,
    templateVertices = compounds[[i]];
    distanceTable = Table[ Apply[normalizedInnerProd, templateVertices[[{i, j} ]]]//Simplify, {i, 5}, {j, 5}];
    Print[distanceTable == IdentityMatrix[5]*(5/4) - 1/4];
];

edge12 = compounds[[1]][[{1, 2}]];
face123 = compounds[[1]][[{1, 2, 3}]];
cell1234 = compounds[[1]][[{1, 2, 3, 4}]];

phi = (1 + Sqrt[5])/2;

mirrors = {
            {0, 1, 0, 0 },
            {1, phi, -1/phi, 0  },
            {1, -phi, 1/phi, 0  },
            {1, 0, phi, -1/phi  }
        };

(* cells = compounds; *)
(* cells = { edge12 }; *)
cells = { face123 };
(* cells = { cell1234 }; *)

Print[cells];

countStart = -1;
(* mirrorCombinations = Subsets[Range[4], {2}]; *)

overlapping[array1_, array2_] := Length[ Intersection[array1, array2, SameTest -> sameVertex] ] == Length[array1];

For[i = 0, i < 100 && countStart != Length[cells], i++,
    countStart = Length[cells];

    For[j = 1, j <= Length[mirrors], j++,

        newCells = Map[reflect[#, mirrors[[j]] ]&, cells, {2}];
        nonOverlappingNewCells = Complement[newCells, cells, SameTest -> overlapping];
        cells = Join[cells, nonOverlappingNewCells];
        Print[cells//Length];
    ];

    Print["Generated " <> IntegerString[Length[cells]] <> " elements so far."];
];

Print[ Intersection[Flatten[cells, 1], polytopeVertices, SameTest -> sameVertex]//Length ];

outputFileName = FileNameJoin[{dataFolder, "compound_720_5_cells_faces.wl"}];

If[FileExistsQ[outputFileName], DeleteFile[outputFileName]];
Save[outputFileName, cells];
Print["Saved polytope component data to " <> outputFileName];

(* 

Compound of 720 5-cells:

600 vertices
7200 edges, one edge can generate all 7200
7200 faces, one face can generate all 7200
3600 cells, one cell can generate all 3600
720 components

 *)