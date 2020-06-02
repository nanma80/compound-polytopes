#!/usr/bin/env wolframscript

numericalValue = True;

epsilon = 1.0*^-5;
dataFolder = "data";

<< Quaternions`

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

compound675x16CellsFileName = FileNameJoin[{dataFolder, "compound_675_16_cells_in_120_cell_componentwise_vertices.wl"}];
compound675x16Cells = Get[compound675x16CellsFileName];

compound75x16CellsFileName = FileNameJoin[{dataFolder, "compound_75_16_cells_in_120_cell_componentwise_vertices.wl"}];
compound75x16Cells = Get[compound75x16CellsFileName];

compound600x16CellsFileName = FileNameJoin[{dataFolder, "compound_600_16_cells_in_120_cell_componentwise_vertices.wl"}];
compound600x16Cells = Get[compound600x16CellsFileName];

compound300x16CellsFileName = FileNameJoin[{dataFolder, "compound_300_16_cells_in_120_cell_componentwise_vertices.wl"}];
compound300x16Cells = Get[compound300x16CellsFileName];

compound75varx16CellsFileName = FileNameJoin[{dataFolder, "compound_75var_a_16_cells_in_120_cell_componentwise_vertices.wl"}];
compound75varx16Cells = Get[compound75varx16CellsFileName];

compound5x600CellsFileName = FileNameJoin[{dataFolder, "compound_5_600_cells_in_120_cell_componentwise_vertices.wl"}];
compound5x600CellsData = Get[compound5x600CellsFileName];

compound5x600Cells = Map[
    Union[Simplify[Flatten[#, 2]], SameTest -> sameVertex]&,
    compound5x600CellsData
];

Length/@compound5x600Cells
(* all 120 *)

compound10x600CellsFileName = FileNameJoin[{dataFolder, "compound_10_600_cells_in_120_cell_componentwise_vertices.wl"}];
compound10x600CellsData = Get[compound10x600CellsFileName];

compound10x600Cells = Map[
    Union[Simplify[Flatten[#, 2]], SameTest -> sameVertex]&,
    compound10x600CellsData
];

Length/@compound10x600Cells
(* all 120 *)

sameSet[array1_, array2_] := Length[ Intersection[array1, array2, SameTest -> sameVertex] ] == Length[array1];
swap[v_] := v[[{1, 2, 4, 3}]];

indices675x5 =
    Table[
        Select[
            Range[Length[compound675x16Cells]],
            Length[Intersection[compound675x16Cells[[#]], compound5x600Cells[[k]], SameTest -> sameVertex]]==Length[compound675x16Cells[[1]]] &
        ],
    {k, 1, 5}
    ];
(* 
In[198]:= Length/@indices675x5                                                                            
Out[198]= {75, 75, 75, 75, 75}
 *)

indices675x10 =
    Table[
        Select[
            Range[Length[compound675x16Cells]],
            Length[Intersection[compound675x16Cells[[#]], compound10x600Cells[[k]], SameTest -> sameVertex]]==Length[compound675x16Cells[[1]]] &
        ],
    {k, 1, Length[compound10x600Cells]}
    ];
(* 
In[199]:= Length/@indices675x10                                                                           
Out[199]= {75, 75, 75, 75, 75, 75, 75, 75, 75, 75}
 *)

Table[
    Length[Intersection@@indices675x10[[{i, j}]]],
    {i, 10}, {j, 10}
]//MatrixForm

(* 
Out[191]//MatrixForm= 75   0    0    0    0    3    3    3    3    3

                      0    75   0    0    0    3    3    3    3    3

                      0    0    75   0    0    3    3    3    3    3

                      0    0    0    75   0    3    3    3    3    3

                      0    0    0    0    75   3    3    3    3    3

                      3    3    3    3    3    75   0    0    0    0

                      3    3    3    3    3    0    75   0    0    0

                      3    3    3    3    3    0    0    75   0    0

                      3    3    3    3    3    0    0    0    75   0

                      3    3    3    3    3    0    0    0    0    75
 *)


Table[
    Length[Intersection@@indices675x5[[{i, j}]]],
    {i, 5}, {j, 5}
]//MatrixForm

(* 
Out[193]//MatrixForm= 75   0    0    0    0

                      0    75   0    0    0

                      0    0    75   0    0

                      0    0    0    75   0

                      0    0    0    0    75
 *)

indices600x5 =
    Table[
        Select[
            Range[Length[compound600x16Cells]],
            Length[Intersection[compound600x16Cells[[#]], compound5x600Cells[[k]], SameTest -> sameVertex]]==Length[compound600x16Cells[[1]]] &
        ],
    {k, 1, 5}
    ];
(* 
In[196]:= Length/@indices600x5                                                                            
Out[196]= {60, 60, 60, 60, 60}
 *)

Table[
    Length[Intersection@@indices600x5[[{i, j}]]],
    {i, 5}, {j, 5}
]//MatrixForm

(* 
Out[200]//MatrixForm= 60   0    0    0    0

                      0    60   0    0    0

                      0    0    60   0    0

                      0    0    0    60   0

                      0    0    0    0    60
 *)

indices600x10 =
    Table[
        Select[
            Range[Length[compound600x16Cells]],
            Length[Intersection[compound600x16Cells[[#]], compound10x600Cells[[k]], SameTest -> sameVertex]]==Length[compound600x16Cells[[1]]] &
        ],
    {k, 1, Length[compound10x600Cells]}
    ];

Length/@indices600x10
(* {60, 60, 60, 60, 60, 60, 60, 60, 60, 60} *)

Table[
    Length[Intersection@@indices600x10[[{i, j}]]],
    {i, 10}, {j, 10}
]//MatrixForm

(* 
Out[202]//MatrixForm= 60   0    0    0    0    0    0    0    0    0

                      0    60   0    0    0    0    0    0    0    0

                      0    0    60   0    0    0    0    0    0    0

                      0    0    0    60   0    0    0    0    0    0

                      0    0    0    0    60   0    0    0    0    0

                      0    0    0    0    0    60   0    0    0    0

                      0    0    0    0    0    0    60   0    0    0

                      0    0    0    0    0    0    0    60   0    0

                      0    0    0    0    0    0    0    0    60   0

                      0    0    0    0    0    0    0    0    0    60
 *)                      


indices75x5 =
    Table[
        Select[
            Range[Length[compound75x16Cells]],
            Length[Intersection[compound75x16Cells[[#]], compound5x600Cells[[k]], SameTest -> sameVertex]]==Length[compound75x16Cells[[1]]] &
        ],
    {k, 1, 5}
    ];

{15, 15, 15, 15, 15}

indices75x10 =
    Table[
        Select[
            Range[Length[compound75x16Cells]],
            Length[Intersection[compound75x16Cells[[#]], compound10x600Cells[[k]], SameTest -> sameVertex]]==Length[compound75x16Cells[[1]]] &
        ],
    {k, 1, 10}
    ];
{15, 15, 15, 15, 15, 15, 15, 15, 15, 15}


indices75varx5 =
    Table[
        Select[
            Range[Length[compound75varx16Cells]],
            Length[Intersection[compound75varx16Cells[[#]], compound5x600Cells[[k]], SameTest -> sameVertex]]==Length[compound75varx16Cells[[1]]] &
        ],
    {k, 1, 5}
    ];
{15, 15, 15, 15, 15}

Table[
    Length[Intersection@@indices75varx5[[{i, j}]]],
    {i, 5}, {j, 5}
]//MatrixForm
(* diagonal matrix with 15 *)

indices75varx10 =
    Table[
        Select[
            Range[Length[compound75varx16Cells]],
            Length[Intersection[compound75varx16Cells[[#]], compound10x600Cells[[k]], SameTest -> sameVertex]]==Length[compound75varx16Cells[[1]]] &
        ],
    {k, 1, 10}
    ];
{15, 15, 15, 15, 15, 0, 0, 0, 0, 0}

Table[
    Length[Intersection@@indices75varx10[[{i, j}]]],
    {i, 10}, {j, 10}
]//MatrixForm

(* my choice of 300 aligns with the latter 5 600-cells *)
indices300x5 =
    Table[
        Select[
            Range[Length[compound300x16Cells]],
            Length[Intersection[compound300x16Cells[[#]], compound5x600Cells[[k]], SameTest -> sameVertex]]==Length[compound300x16Cells[[1]]] &
        ],
    {k, 1, 5}
    ];
Length/@indices300x5
{0, 0, 0, 0, 0}

indices300x10 =
    Table[
        Select[
            Range[Length[compound300x16Cells]],
            Length[Intersection[compound300x16Cells[[#]], compound10x600Cells[[k]], SameTest -> sameVertex]]==Length[compound300x16Cells[[1]]] &
        ],
    {k, 1, 10}
    ];
Length/@indices300x10
{0, 0, 0, 0, 0, 60, 60, 60, 60, 60}


