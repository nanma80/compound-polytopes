#!/usr/bin/env wolframscript

numericalValue = True;

shapes = {
    (* "polytope_3_3_4" *)
    (* "polytope_3_4_3_alt" *)
    (* "polytope_3_3_5" *)

    (* "polytope_3_3_3" *)
    (* "polytope_3_4_3" *)
    (* "polytope_4_3_3" *)
    "polytope_5_3_3"
};

epsilon = 1.0*^-5;
dataFolder = "data";

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

For[shapeIndex = 1, shapeIndex <= Length[shapes], shapeIndex++,
    shape = shapes[[shapeIndex]];
    Print["shape: " <> shape];
    polytopeFileName = FileNameJoin[{dataFolder, shape <> ".wl"}];

    polytopeData = Get[polytopeFileName];
    polytopeEdges = polytopeData["edges"];
    polytopeFaces = polytopeData["faces"];
    polytopeVertices = polytopeData["vertices"];
    polytopeCells = polytopeData["cells"];

    edgeCenters = Map[Mean, polytopeEdges];
    faceCenters = Map[Mean, polytopeFaces];
    cellCenters = Map[Mean[Mean[#]]&, polytopeCells];

    vInnerProds = Sort[getInnerProds[polytopeVertices]//N];
    eInnerProds = Sort[getInnerProds[edgeCenters]//N];
    fInnerProds = Sort[getInnerProds[faceCenters]//N];
    cInnerProds = Sort[getInnerProds[cellCenters]//N];

    Print["vertex count:"];
    Print[Length[polytopeVertices]];

    Print["edge count:"];
    Print[polytopeEdges//Length];
    Print["vertex inner prods contained by edge inner prods"];
    Print[Length[vInnerProds] == Length[Intersection[vInnerProds, eInnerProds, SameTest -> (Abs[N[#1] - N[#2]] < epsilon&)]] ];

    Print["face count:"];
    Print[polytopeFaces//Length];
    Print["vertex inner prods contained by face inner prods"];
    Print[Length[vInnerProds] == Length[Intersection[vInnerProds, fInnerProds, SameTest -> (Abs[N[#1] - N[#2]] < epsilon&)]] ];

    Print["cell count:"];
    Print[polytopeCells//Length];
    Print["vertex inner prods contained by cell inner prods"];
    Print[Length[vInnerProds] == Length[Intersection[vInnerProds, cInnerProds, SameTest -> (Abs[N[#1] - N[#2]] < epsilon&)]] ];
];

(*
333: inner prod of vertices not in face centers
433: some inner prods of vertices not in edge centers
533: some inner prods of vertices not in edge centers
*)