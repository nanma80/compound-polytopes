#!/usr/bin/env wolframscript

numericalValue = True;

shapes = {
    (* "polytope_3_3_4" *)
    (* "polytope_3_4_3_alt" *)
    "polytope_3_3_5"
    (* , "polytope_3_4_3" *)
    (* , "polytope_4_3_3" *)
    (* , "polytope_5_3_3" *)
};

exportToPov = True;
imageSize = {600, 450};
perspective = True;

(* rangeFactor = 0.8; *)
rangeFactor = 0.6;

Needs["POVRayRender`"];
ConfigurePOVRayRender["POVRayPath" -> "C:\\Program Files\\POV-Ray\\v3.7\\bin\\pvengine64.exe"];

epsilon = 1.0*^-5;
iterationLimit = 20;
dataFolder = "data";
imageFolder = "output";

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

For[shapeIndex = 1, shapeIndex <= Length[shapes], shapeIndex++,
    shape = shapes[[shapeIndex]];
    Print["shape: " <> shape];
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

    faceCenters = Mean/@polytopeFaces;
    normFaceCenter = Norm[faceCenters[[1]]];

    vertex1 = SortBy[polytopeVertices, -#[[4]]&][[1]];
    normVertex = Norm[vertex1];
    vertex1 = vertex1/Norm[vertex1]*normFaceCenter//FullSimplify;
    Print["vertex1: "];
    Print[vertex1];
    normalizedVertices = polytopeVertices/normVertex * normFaceCenter//FullSimplify;

    seedFaceCenters = {};
    componentVertices = {};
    componentEdges = {};
    componentFaces = {};
    componentCells = {};
    otherVectors = IdentityMatrix[4];

    If[shape == "polytope_3_3_4",
        seedFaceCenters = Select[faceCenters, #[[4]] == 0 && #[[3]] > 0&];
    ];

    If[shape == "polytope_3_4_3_alt",
        seedFaceCenters = Select[faceCenters, #[[4]] == 4/3 && #[[3]] == 2/3 &];
    ];

    If[shape == "polytope_3_3_5",
        seedFaceCenters = Select[faceCenters, Simplify[#[[4]]] == (1 + Sqrt[5])/2 && #[[1]] + #[[2]]*0.1 + #[[3]]*0.01 + #[[4]] * 0.001> 0 &]//Simplify;
        otherVectors = {{1, 0, 0, 0}, {0, 0, 1, 0}};
    ];

    (* Map[#[[4]]&, faceCenters]//Tally *)

    Print["Number of seedFaceCenters:"];
    Print[seedFaceCenters//Length];

    For[seedIndex = 1, seedIndex <= Length[seedFaceCenters], seedIndex++,
        Print[seedIndex];
        faceCenter1 = seedFaceCenters[[seedIndex]];
        Print[faceCenter1];
        If[shape == "polytope_3_3_5" && seedIndex == 9,
            otherVectors = {{1, 0, 0, 0}, {0, 1, 0, 0}};
        ];
        orth4dRaw = Orthogonalize[Join[{vertex1, faceCenter1}, otherVectors]]//Simplify;
        orth4d = Select[orth4dRaw, Norm[#] > epsilon&];

        If[Length[orth4d] < 4,
            Print["length of orth4d is: " <> IntegerString[Length[orth4d]]];
            Throw[singular1];
        ];

        If[Abs[Det[orth4d]] < epsilon, 
            Print["orth4d is singular. Unexpected"];
            Throw[singular2];
        ];

        If[Det[orth4d] < 0, orth4d = {orth4d[[1]], orth4d[[2]], orth4d[[4]], orth4d[[3]]}];
        Print["orth4d computed"];

        theta = ArcCos[faceCenter1.vertex1/Norm[vertex1]/Norm[faceCenter1]];
        rot4d = FullSimplify[Transpose[orth4d].doubleRotation[-theta].orth4d];

        newVertices = polytopeVertices.rot4d / normVertex * normFaceCenter//FullSimplify;
        newEdges = Map[#.rot4d&, polytopeEdges] / normVertex * normFaceCenter//FullSimplify;
        newFaces = Map[#.rot4d&, polytopeFaces] / normVertex * normFaceCenter//FullSimplify;
        newCells = Map[#.rot4d&, polytopeCells, {2}] / normVertex * normFaceCenter//FullSimplify;

        AppendTo[componentVertices, newVertices];
        AppendTo[componentEdges, newEdges];
        AppendTo[componentFaces, newFaces];
        AppendTo[componentCells, newCells];
    ];

    Print["number of components: " <> IntegerString[Length[componentVertices]] <> ". Each component contains: " <> IntegerString[Length[componentVertices[[1]]]] <> " vertices"];
    allComponentVertices = Union[Flatten[N[componentVertices], 1], SameTest -> SameMatrix];
    Print["all components done"];
    ints = Intersection[N[faceCenters], N[allComponentVertices], SameTest -> SameMatrix];
    Print[Length/@{allComponentVertices, faceCenters, ints}];

    polytopeData["componentVertices"] = componentVertices;
    polytopeData["componentEdges"] = componentEdges;
    polytopeData["componentFaces"] = componentFaces;
    polytopeData["componentCells"] = componentCells;

    Save[polytopeFileName, polytopeData];
    Print["Saved polytope component data to " <> polytopeFileName];
];
