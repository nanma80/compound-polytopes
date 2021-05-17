#!/usr/bin/env wolframscript


rangeFactor = 0.75;
dataFolder = "data";
imageFolder = "output";
exportToPov = True;
imageSize = 240 * {4, 3};
epsilon = 0.000001;

shapes = {
    "compound_120_5_cells_in_120_cell",
    "compound_120var_a_5_cells_in_120_cell",
    "compound_720_5_cells_in_120_cell",
    "compound_840_5_cells_in_120_cell"
};

opacity = 1;

phi = GoldenRatio;
viewPoint = {2, -2, 2}/Sqrt[3];
targetViewPoint = {1, -phi, 0};
rotation3D = RotationMatrix[{targetViewPoint, viewPoint}];

Needs["POVRayRender`"];
ConfigurePOVRayRender["POVRayPath" -> "C:\\Program Files\\POV-Ray\\v3.7\\bin\\pvengine64.exe"];

getProjectionMatrix[direction_] := Simplify[Orthogonalize[{direction, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}}]][[2;;4]];

plotImage[graphElements_, range_, outputFileName_, imageSize_] := Module[
    {graph, lighting, cropCommand},

    lighting = {{"Point", White, viewPoint * 10}};
    graph = Graphics3D[
        graphElements,
        ViewPoint -> viewPoint * 1000,
        PlotRange -> {{-range, range}, {-range, range}, {-range, range}},
        Lighting -> lighting
    ];

    cropCommand = "mogrify -crop " <> IntegerString[imageSize[[2]]] <> "x" <> IntegerString[imageSize[[2]]] <> "+" <> IntegerString[(imageSize[[1]]-imageSize[[2]])/2] <> "+0 " <> outputFileName;

    If[exportToPov,
        If[FileExistsQ[outputFileName], DeleteFile[outputFileName]];
        POVRayRender[graph, Method -> "Triangle", OutputPath -> outputFileName, ImageSize -> imageSize, RenderOptions -> "+A0.001 -J +UA"];
        Run[cropCommand];
        Print[cropCommand];
        ,
        Export[outputFileName, graph]
    ];
    Print["Saved the image to " <> outputFileName];
];

projectionVector4D = {1, 0, 0, 0} + {0.1, 0.13, 0.61, 0.21} * 0.001;
projectionTo3D = rotation3D.getProjectionMatrix[projectionVector4D];

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

For[shapeIndex = 1, shapeIndex <= Length[shapes], shapeIndex ++,
    shape = shapes[[shapeIndex]];

    dataFileName = FileNameJoin[{dataFolder, shape <> "_componentwise_vertices.wl"}];

    Print["Reading data for " <> shape <> " from " <> dataFileName];

    inputData = Get[dataFileName];

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

    components = Map[get5Cells, inputData];
    facesPerComponent = Map[ Flatten[#, 1]&, components];
    faceCountPerComponent = Length[facesPerComponent[[1]]];
    faces = Flatten[facesPerComponent, 1];
    maxNorm = Norm[faces[[1]][[1]]];

    outputFolder = FileNameJoin[{imageFolder, "compound_5_cells"}];
    If[!DirectoryQ[outputFolder], CreateDirectory[outputFolder]];

    thickness = maxNorm * 0.02;
    range = maxNorm * rangeFactor;

    outputFileName = shape <> ".png";
    outputFileName = FileNameJoin[{outputFolder, outputFileName}];

    graphElements = {
            Table[
                {
                    Opacity[opacity], 
                    FaceForm[
                        Hue[
                            Floor[(k - 1)/faceCountPerComponent]/Length[components]
                        ]
                    ],
                    Polygon[ Map[projectionTo3D.#& , faces[[k]] ] ]
                },
                {k, 1, Length[faces]}
            ]
        };

    plotImage[graphElements, range, outputFileName, imageSize];
];
