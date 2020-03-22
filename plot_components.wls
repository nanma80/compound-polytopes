#!/usr/bin/env wolframscript

shapes = {
    (* "polytope_3_3_4" *)
    "polytope_3_4_3_alt"
    , "polytope_3_3_5"
};

rangeFactor = 0.75;
dataFolder = "data";
imageFolder = "output";
exportToPov = True;
plotAllCells = True;
imageSize = 150 * {4, 3};
(*imageSize = {1000, 750};*)

Needs["POVRayRender`"];
ConfigurePOVRayRender["POVRayPath" -> "C:\\Program Files\\POV-Ray\\v3.7\\bin\\pvengine64.exe"];

colors = Join[{Red, Blue, Green, Yellow, Magenta, Cyan, Brown, Orange, Pink, Purple}, RandomColor[10]];
getProjectionMatrix[direction_] := Simplify[Orthogonalize[{direction, {1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}}]][[2;;4]];
projectTo3d[vertex4d_] := {vertex4d[[1]], vertex4d[[2]], vertex4d[[3]]};
reflect[point_, mirror_] := FullSimplify[point - 2 * (point.mirror)/(mirror.mirror)*mirror];
doubleReflect[point_, mirror1_, mirror2_] := reflect[reflect[point, mirror1], mirror2];
getEdgesFromFace[face_] := Table[{face[[i+1]], face[[Mod[i+1, Length[face]] + 1]]}, {i, 0, Length[face] - 1}];
explodedFace[face_, explosionFactor_] := Map[(# + Mean[face] * explosionFactor)&, face];
sameVertex[point1_, point2_] := point1 == point2;
sameCenter[edge1_, edge2_] := Total[edge1] == Total[edge2];
sameCellCenter[cell1_, cell2_] := getCellCenter[cell1] == getCellCenter[cell2];
getCellCenter[cell_] := Total[Flatten[cell, 1]];

getGraphElements[vertices_, edges_, faces_, color_, thickness_, opacity_] :=
{
    {
        Opacity[opacity], FaceForm[color], Map[Polygon, Map[projectionTo3D.#&, faces, {2}]]
    },
    {
        Opacity[1], FaceForm[color], Map[Cylinder[#, thickness]&, Map[projectionTo3D.#&, edges, {2}]]
    },
    {
        Opacity[1], FaceForm[color], Map[Sphere[#, thickness * 2]&, Map[projectionTo3D.#&, vertices]]
    }
};

plotImage[graphElements_, range_, outputFileName_, imageSize_] := Module[
    {graph, lighting, cropCommand},

    lighting = {{"Point", White, {10, -10, 10}}};
    graph = Graphics3D[
        graphElements,
        ViewPoint -> {5, -3, 3}*3,
        PlotRange -> {{-range, range}, {-range, range}, {-range, range}},
        Lighting -> lighting
    ];

    cropCommand = "mogrify -crop " <> IntegerString[imageSize[[2]]] <> "x" <> IntegerString[imageSize[[2]]] <> "+" <> IntegerString[(imageSize[[1]]-imageSize[[2]])/2] <> "+0 " <> outputFileName;

    If[exportToPov,
        If[FileExistsQ[outputFileName], DeleteFile[outputFileName]];
        POVRayRender[graph, Method -> "Triangle", OutputPath -> outputFileName, ImageSize -> imageSize, RenderOptions -> "+A0.001 -J +UA"];
        Run[cropCommand];
        ,
        Export[outputFileName, graph]
    ];
    Print["Saved the image to " <> outputFileName];
];

projectionVector4D = {0, 0, 0, 1} + {0.1, 0.13, 0.61, 0.21} * 0.0001;
projectionTo3D = getProjectionMatrix[projectionVector4D];

For[opacity = 0.0, opacity <= 1, opacity += 0.5,
    Print["opacity: " <> ToString[N[opacity]]];
    For[shapeIndex = 1, shapeIndex <= Length[shapes], shapeIndex ++,
        shape = shapes[[shapeIndex]];
        Print["shape: " <> shape];
        polytopeFileName = FileNameJoin[{dataFolder, shape <> ".wl"}];

        polytopeData = Get[polytopeFileName];

        maxNorm = Norm[polytopeData["vertices"][[1]]];

        outputFolder = FileNameJoin[{imageFolder, shape}];
        If[!DirectoryQ[outputFolder], CreateDirectory[outputFolder]];


        (* original *)
        maxNorm = Norm[polytopeData["vertices"][[1,1]]];
        thickness = maxNorm * 0.02;
        range = maxNorm * rangeFactor;

        outputFileName = "original" <> 
            "_op_" <> IntegerString[Round[opacity * 100]] <> 
            "_" <> shape <> ".png";
        outputFileName = FileNameJoin[{outputFolder, outputFileName}];
        graphElements = getGraphElements[
                            polytopeData["vertices"],
                            polytopeData["edges"],
                            polytopeData["faces"],
                            LightGray,
                            thickness,
                            opacity
                        ];
        plotImage[graphElements, range, outputFileName, imageSize];

        (* components *)
        keys = Keys[polytopeData];
        For[keyIndex = 1, keyIndex <= Length[keys], keyIndex++,
            key = keys[[keyIndex]];
            For[signature = 1, signature >= -1, signature = signature - 2,
                signatureString = If[signature == 1, "1", "2"];
                If[StringContainsQ[key, "ComponentVertices" <> signatureString],
                (* If[StringContainsQ[key, "faceComponentVertices" <> signatureString]  && signatureString == "2", *)
                    Print[key];
                    prefix = StringTake[key, 13];
                    suffix = StringTake[key, -1];
                    componentCount = Length[polytopeData[key]];
                    Print["componentCount: " <> IntegerString[componentCount]];

                    maxNorm = Norm[polytopeData[key][[1,1]]];
                    thickness = maxNorm * 0.02;
                    range = maxNorm * rangeFactor;
                    allGraphElements = {};
                    For[cIndex = 1, cIndex <= componentCount, cIndex++,
                        graphElements = getGraphElements[
                            polytopeData[prefix <> "Vertices" <> suffix][[cIndex]],
                            polytopeData[prefix <> "Edges" <> suffix][[cIndex]],
                            polytopeData[prefix <> "Faces" <> suffix][[cIndex]],
                            colors[[cIndex]],
                            thickness,
                            opacity
                        ];
                        AppendTo[allGraphElements, graphElements];
                    ];

                    outputFileName = prefix <> "_" <> signatureString <> 
                        "_op_" <> IntegerString[Round[opacity * 100]] <> 
                        "_" <> shape <> ".png";
                    outputFileName = FileNameJoin[{outputFolder, outputFileName}];
                    plotImage[Flatten[allGraphElements, 1], range, outputFileName, imageSize];
                ];
            ];
        ];
    ];
];
