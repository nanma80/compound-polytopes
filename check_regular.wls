#!/usr/bin/env wolframscript

numericalValue = True;

shapes = {
    "polytope_3_3_4"
    , "polytope_3_4_3_alt"
    , "polytope_3_3_5"
};

epsilon = 1.0*^-5;
dataFolder = "data";

reflect[point_, mirror_] := FullSimplify[point - 2 * (point.mirror)/(mirror.mirror)*mirror];
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

getIntersectionMatrix[points1_, points2_] := Table[
    Length[Intersection[points1[[i]], points2[[j]]]], 
    {i, Length[points1]}, {j, Length[points2]}];

For[shapeIndex = 1, shapeIndex <= Length[shapes], shapeIndex++,
    shape = shapes[[shapeIndex]];
    Print["shape: " <> shape];
    polytopeFileName = FileNameJoin[{dataFolder, shape <> ".wl"}];

    polytopeData = Get[polytopeFileName];

    edgeComponentEdges = polytopeData["edgeComponentEdges1"];
    edgeCentersComponent = Map[Mean, edgeComponentEdges, {2}];
    edgeCenters = Flatten[edgeCentersComponent, 1]//Simplify;
    (* 96 edge centers. 4 components * 24 edges/component *)
    (* All unique, {2, 1, 1, 0}/6, all perms, all signs *)
    ta = Tally[edgeCenters];
    Print["edge center counts"];
    Print[Length[edgeCenters]];
    Print["edge centers of the edge Component, tally"];
    Print[Tally[Map[#[[2]]&, ta]]];


    faceComponentFaces = polytopeData["faceComponentFaces1"];
    faceCentersComponent = Map[Mean, faceComponentFaces, {2}];
    faceCenters = Flatten[faceCentersComponent, 1]//Simplify;
    (* 128 face centers. 4 components * 32 faces/component *)
    (* 96 face centers appearing in one component + 8 face centers shared in all 4 components *)
    (* appearing once: {2, 2, 1, 0}/9, all perms, all signs *)
    (* appearing 4 times: {1, 0, 0, 0}/3, all perms, all signs *)
    (* therefore the component is not face transitive *)
    ta = Tally[faceCenters];
    Print["face center counts"];
    Print[Length[faceCenters]];
    Print["face centers of the face Component, tally"];
    Print[Tally[Map[#[[2]]&, ta]]];

    Print[];
];


(*
16-cell:
edge center based components: 3 16-cells in 24-cell. This is regular.
face center based components: not regular. Some faces shared by 4 components, others only 1 component

24-cell, 600-cell:
edge and face center based components: not regular

shape: polytope_3_3_4
edge center counts
72
edge centers of the edge Component, tally
{{3, 24}}
face center counts
128
face centers of the face Component, tally
{{1, 96}, {4, 8}}

shape: polytope_3_4_3_alt
edge center counts
384
edge centers of the edge Component, tally
{{1, 288}, {4, 24}}
face center counts
384
face centers of the face Component, tally
{{1, 288}, {4, 24}}

shape: polytope_3_3_5
edge center counts
4320
edge centers of the edge Component, tally
{{1, 3600}, {6, 120}}
face center counts
12000
face centers of the face Component, tally
{{10, 120}, {1, 10800}}

*)