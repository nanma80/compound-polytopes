#!/usr/bin/env wolframscript

numericalValue = True;

exportToPov = True;
imageSize = {600, 450};
perspective = True;

(* rangeFactor = 0.8; *)
rangeFactor = 0.6;

Needs["POVRayRender`"];
ConfigurePOVRayRender["POVRayPath" -> "C:\\Program Files\\POV-Ray\\v3.7\\bin\\pvengine64.exe"];

epsilon = 1*^-5;
iterationLimit = 20;
dataFolder = "data";
imageFolder = "output";

v3 = {0, 1, 1, 0};
v4 = {0, 0, 0, 1};

perspectiveProjection[point_, fullProjectionMatrix_, wDistance_ ] := projectToDirection[point, fullProjectionMatrix[[1;;3]]] / (1 + point.(fullProjectionMatrix[[4]])/wDistance);
projectToDirection[vertex4d_, projectionMatrix_] := projectionMatrix.vertex4d;
getEdgesFromFace[face_] := Table[{face[[i+1]], face[[Mod[i+1, Length[face]] + 1]]}, {i, 0, Length[face] - 1}];
sameVertex[point1_, point2_] := point1 == point2;
cubeVertices[dimension_] := 2 * Map[IntegerDigits[#, 2][[2;;dimension+1]]&, Table[k, {k, 2^dimension, 2 * 2^dimension - 1}]] - 1;
allPlusMinus[v_] := Union[Map[(v #)&, cubeVertices[Length[v]]]];
SameMatrix[m1_, m2_] := ( m1 == m2 ) || (Max[Max[N[Abs[m1 - m2]]]] < epsilon);
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

shapes = {
    "polytope_3_3_4"
(*     , "polytope_3_4_3"
    , "polytope_3_4_3_alt"
    , "polytope_3_3_5"
    , "polytope_4_3_3"
    , "polytope_5_3_3" *)
};

For[shapeIndex = 1, shapeIndex <= Length[shapes], shapeIndex++,
    shape = shapes[[shapeIndex]];
    polytopeCellsFileName = FileNameJoin[{dataFolder, "cells", shape <> ".wl"}];
    polytopeFileName = FileNameJoin[{dataFolder, shape <> ".wl"}];

    If[!FileExistsQ[polytopeFileName],

        Print["Generating polytope data for " <> polytopeCellsFileName];
        polytopeCells = Get[polytopeCellsFileName];
        Print[polytopeCells//Length]; (* 24 *)
        polytopeFaces = Union[ Flatten[polytopeCells, 1], SameTest -> sameCenter];
        Print[polytopeFaces//Length]; (* 96 *)
        polytopeEdges = Union[ Flatten[Map[getEdgesFromFace, polytopeFaces], 1], SameTest -> sameCenter];
        Print[polytopeEdges//Length]; (* 96 *)
        polytopeVertices = Union[ Flatten[polytopeFaces,1], SameTest -> sameVertex];
        Print[polytopeVertices//Length]; (* 24 *)

        polytopeData = <|
            "cells" -> polytopeCells,
            "faces" -> polytopeFaces,
            "edges" -> polytopeEdges,
            "vertices" -> polytopeVertices
        |>;

        Save[polytopeFileName, polytopeData];
        Print["Saved polytope data to " <> polytopeFileName];
    ];

    polytopeData = Get[polytopeFileName];
    polytopeEdges = polytopeData["edges"];
    polytopeFaces = polytopeData["faces"];
    polytopeVertices = polytopeData["vertices"];
    polytopeCells = polytopeData["cells"];

(*     Print[Length[polytopeCells]];
    Print[Length[polytopeFaces]];
    Print[Length[polytopeEdges]];
    Print[Length[polytopeVertices]]; *)

    vertex1 = polytopeVertices[[1]];
    Print[vertex1];
    normVertex = Norm[vertex1];

    faceCenters = Total/@polytopeFaces;
    normalizedFaceCenters = faceCenters/Norm[faceCenters[[1]]]*normVertex;
    faceCenter1 = normalizedFaceCenters[[1]];
    Print[faceCenter1];

    orth4d = Orthogonalize[{vertex1, faceCenter1, v3, v4}];
    If[Abs[Det[orth4d]] < epsilon, 
        Print["orth4d is singular"];
        Throw[singular];
        ];

    If[Det[orth4d] < 0, orth4d = Orthogonalize[{vertex1, faceCenter1, v4, v3}]];

    Print[orth4d[[1]]*normVertex - vertex1];
    Print[orth4d[[2]]*normVertex - faceCenter1//N];

    theta = ArcCos[faceCenter1.vertex1/Norm[vertex1]/Norm[faceCenter1]];
    rot4d = Transpose[orth4d].doubleRotation[-theta].orth4d;
    Print[rot4d//Simplify];

    newVertices = polytopeVertices.rot4d;

    ints = Intersection[normalizedFaceCenters, newVertices, SameTest -> SameMatrix];
    Print[Length/@{newVertices, normalizedFaceCenters, ints}];
];
