(* Content-type: application/vnd.wolfram.cdf.text *)

(*** Wolfram CDF File ***)
(* http://www.wolfram.com/cdf *)

(* CreatedBy='Mathematica 10.0' *)

(*************************************************************************)
(*                                                                       *)
(*  The Mathematica License under which this file was created prohibits  *)
(*  restricting third parties in receipt of this file from republishing  *)
(*  or redistributing it by any means, including but not limited to      *)
(*  rights management or terms of use, without the express consent of    *)
(*  Wolfram Research, Inc. For additional information concerning CDF     *)
(*  licensing and redistribution see:                                    *)
(*                                                                       *)
(*        www.wolfram.com/cdf/adopting-cdf/licensing-options.html        *)
(*                                                                       *)
(*************************************************************************)

(*CacheID: 234*)
(* Internal cache information:
NotebookFileLineBreakTest
NotebookFileLineBreakTest
NotebookDataPosition[      1064,         20]
NotebookDataLength[     19237,        638]
NotebookOptionsPosition[     18602,        597]
NotebookOutlinePosition[     19023,        613]
CellTagsIndexPosition[     18980,        610]
WindowFrame->Normal*)

(* Beginning of Notebook Content *)
Notebook[{

Cell[CellGroupData[{
Cell["Reduced density matrix (RDM) demo", "Section"],

Cell["\<\
FermiFab toolbox demonstration file: calculate the 2-body reduced density \
matrix of the outer product of two Chromium 5D symmetry eigenstates.\
\>", "Text"],

Cell[TextData[{
 "Copyright (c) 2008-2015, Christian B. Mendl\nAll rights reserved.\n",
 ButtonBox["http://christian.mendl.net",
  BaseStyle->"Hyperlink",
  ButtonData->{
    URL["http://christian.mendl.net"], None},
  ButtonNote->"http://christian.mendl.net"]
}], "Text"],

Cell[TextData[{
 "This program is free software; you can redistribute it and/or\nmodify it \
under the terms of the Simplified BSD License\n",
 ButtonBox["http://www.opensource.org/licenses/bsd-license.php",
  BaseStyle->"Hyperlink",
  ButtonData->{
    URL["http://www.opensource.org/licenses/bsd-license.php"], None},
  ButtonNote->"http://www.opensource.org/licenses/bsd-license.php"]
}], "Text"],

Cell[TextData[{
 "Reference:\n\tChristian B. Mendl\n\tThe FermiFab toolbox for fermionic \
many-particle quantum systems\n\tComput. Phys. Commun.182,1327-1337 (2011) ",
 ButtonBox["DOI",
  BaseStyle->"Hyperlink",
  ButtonData->{
    URL["http://dx.doi.org/10.1016/j.cpc.2011.01.028"], None},
  ButtonNote->"http://dx.doi.org/10.1016/j.cpc.2011.01.028"],
 "\n\tpreprint ",
 ButtonBox["http://arxiv.org/abs/1103.0872",
  BaseStyle->"Hyperlink",
  ButtonData->{
    URL["http://arxiv.org/abs/1103.0872"], None},
  ButtonNote->"http://arxiv.org/abs/1103.0872"]
}], "Text"],

Cell[BoxData[{
 RowBox[{
  RowBox[{"SetDirectory", "[", 
   RowBox[{"NotebookDirectory", "[", "]"}], "]"}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"<<", "fermifab.m"}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"fermifabML", "=", 
   RowBox[{"Install", "[", 
    RowBox[{
    "\"\<mlink/fermifabML/\>\"", "<>", "$SystemID", "<>", 
     "\"\</fermifabML\>\""}], "]"}]}], ";"}]}], "Input"],

Cell[BoxData[{
 RowBox[{
  RowBox[{
   RowBox[{
    SubscriptBox["quantL", "Cr"], "=", 
    RowBox[{"{", 
     RowBox[{"2", ",", "0", ",", "1", ",", "2"}], "}"}]}], ";"}], " ", 
  RowBox[{"(*", " ", 
   RowBox[{
    RowBox[{"L", " ", "quantum", " ", "number", " ", "of", " ", "3", "d"}], 
    ",", 
    RowBox[{"4", "s"}], ",", 
    RowBox[{"4", "p"}], ",", 
    RowBox[{"4", "d", " ", "orbitals"}]}], " ", 
   "*)"}]}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{
   RowBox[{
    SubscriptBox["num", "filled"], "=", 
    RowBox[{"{", 
     RowBox[{"2", ",", "2", ",", "6", ",", "2", ",", "6"}], "}"}]}], ";"}], 
  " ", 
  RowBox[{"(*", " ", 
   RowBox[{
    RowBox[{
    "particle", " ", "numbers", " ", "in", " ", "fully", " ", "occupied", " ",
      "orbitals", " ", "1", "s"}], ",", 
    RowBox[{"2", "s"}], ",", 
    RowBox[{"2", "p"}], ",", 
    RowBox[{"3", "s"}], ",", 
    RowBox[{"3", "p"}]}], " ", "*)"}]}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{
   SubscriptBox["num", "orbs"], "=", 
   RowBox[{"Join", "[", 
    RowBox[{
     SubscriptBox["num", "filled"], ",", 
     RowBox[{"2", 
      RowBox[{"(", 
       RowBox[{
        RowBox[{"2", 
         SubscriptBox["quantL", "Cr"]}], "+", "1"}], ")"}]}]}], "]"}]}], 
  ";"}]}], "Input"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{
    RowBox[{
     RowBox[{
      RowBox[{
       RowBox[{
        RowBox[{
         RowBox[{"subshell", " ", "names"}], ";"}], " ", "'"}], "~"}], "'"}], 
     " ", "stands", " ", "for", " ", "spin", " ", "down", " ", 
     RowBox[{"(", 
      RowBox[{"otherwise", " ", "up"}], ")"}]}], ",", " ", 
    RowBox[{"i", ".", "e", "."}], ",", " ", 
    RowBox[{
     RowBox[{"dm", "~", " ", "corresponds"}], " ", "to", " ", 
     RowBox[{
      SubscriptBox["d", "m"], "\[DownArrow]"}]}]}], " ", "*)"}], 
  "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    RowBox[{
     SubscriptBox["subshells", "0"], "=", 
     RowBox[{"StringSplit", "[", "\"\<s s~\>\"", "]"}]}], ";"}], 
   "\[IndentingNewLine]", 
   RowBox[{
    RowBox[{
     SubscriptBox["subshells", "1"], "=", 
     RowBox[{"StringSplit", "[", "\"\<pz pz~ px px~ py py~\>\"", "]"}]}], 
    ";"}], "\[IndentingNewLine]", 
   RowBox[{
    RowBox[{
     SubscriptBox["subshells", "2"], "=", 
     RowBox[{
     "StringSplit", "[", "\"\<d0 d0~ dz dz~ dm dm~ dx dx~ dy dy~\>\"", 
      "]"}]}], ";"}], "\[IndentingNewLine]", 
   RowBox[{
    SubscriptBox["shells", "active"], "=", 
    RowBox[{"Join", "[", 
     RowBox[{
      RowBox[{
       RowBox[{
        RowBox[{"\"\<3\>\"", "<>", "#"}], "&"}], "/@", 
       SubscriptBox["subshells", "2"]}], ",", 
      RowBox[{
       RowBox[{
        RowBox[{"\"\<4\>\"", "<>", "#"}], "&"}], "/@", 
       SubscriptBox["subshells", "0"]}], ",", 
      RowBox[{
       RowBox[{
        RowBox[{"\"\<4\>\"", "<>", "#"}], "&"}], "/@", 
       SubscriptBox["subshells", "1"]}], ",", 
      RowBox[{
       RowBox[{
        RowBox[{"\"\<4\>\"", "<>", "#"}], "&"}], "/@", 
       SubscriptBox["subshells", "2"]}]}], "]"}]}], "\[IndentingNewLine]", 
   RowBox[{
    RowBox[{
     SubscriptBox["shells", "all"], "=", 
     RowBox[{"Flatten", "[", 
      RowBox[{
       RowBox[{"Function", "[", 
        RowBox[{"n", ",", 
         RowBox[{
          RowBox[{"Function", "[", 
           RowBox[{"l", ",", 
            RowBox[{
             RowBox[{
              RowBox[{
               RowBox[{"ToString", "[", "n", "]"}], "<>", "#"}], "&"}], "/@", 
             SubscriptBox["subshells", "l"]}]}], "]"}], "/@", 
          RowBox[{"Range", "[", 
           RowBox[{"0", ",", 
            RowBox[{"Min", "[", 
             RowBox[{"2", ",", 
              RowBox[{"n", "-", "1"}]}], "]"}]}], "]"}]}]}], "]"}], "/@", 
       RowBox[{"Range", "[", "4", "]"}]}], "]"}]}], ";"}], 
   "\[IndentingNewLine]", 
   RowBox[{
    RowBox[{
     SubscriptBox["shells", "spatial"], "=", 
     RowBox[{
      SubscriptBox["shells", "all"], "\[LeftDoubleBracket]", 
      RowBox[{"1", ";;", 
       RowBox[{"-", "1"}], ";;", "2"}], "\[RightDoubleBracket]"}]}], 
    ";"}]}]}]], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{"\<\"3d0\"\>", ",", "\<\"3d0~\"\>", ",", "\<\"3dz\"\>", 
   ",", "\<\"3dz~\"\>", ",", "\<\"3dm\"\>", ",", "\<\"3dm~\"\>", 
   ",", "\<\"3dx\"\>", ",", "\<\"3dx~\"\>", ",", "\<\"3dy\"\>", 
   ",", "\<\"3dy~\"\>", ",", "\<\"4s\"\>", ",", "\<\"4s~\"\>", 
   ",", "\<\"4pz\"\>", ",", "\<\"4pz~\"\>", ",", "\<\"4px\"\>", 
   ",", "\<\"4px~\"\>", ",", "\<\"4py\"\>", ",", "\<\"4py~\"\>", 
   ",", "\<\"4d0\"\>", ",", "\<\"4d0~\"\>", ",", "\<\"4dz\"\>", 
   ",", "\<\"4dz~\"\>", ",", "\<\"4dm\"\>", ",", "\<\"4dm~\"\>", 
   ",", "\<\"4dx\"\>", ",", "\<\"4dx~\"\>", ",", "\<\"4dy\"\>", 
   ",", "\<\"4dy~\"\>"}], "}"}]], "Output"]
}, Open  ]],

Cell[CellGroupData[{

Cell["\<\
Define two Fermi states (these are Chromium 5D symmetry eigenstates)\
\>", "Subsection"],

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{
   "encode", " ", "shell", " ", "names", " ", "as", " ", "Fermi", " ", "bit", 
    " ", "patterns"}], " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{"EncodeShells", "[", "sn_", "]"}], ":=", 
   RowBox[{"Ket", "[", 
    RowBox[{"CoordsToFermi", "[", 
     RowBox[{
      RowBox[{"Range", "[", 
       RowBox[{"Total", "[", 
        SubscriptBox["num", "filled"], "]"}], "]"}], "\[Union]", 
      RowBox[{"(", 
       RowBox[{
        RowBox[{"Total", "[", 
         SubscriptBox["num", "filled"], "]"}], "+", 
        RowBox[{"Flatten", "[", 
         RowBox[{
          RowBox[{
           RowBox[{"Position", "[", 
            RowBox[{
             SubscriptBox["shells", "active"], ",", "#"}], "]"}], "&"}], "/@", 
          RowBox[{"StringSplit", "[", "sn", "]"}]}], "]"}]}], ")"}]}], "]"}], 
    "]"}]}]}]], "Input"],

Cell[BoxData[{
 RowBox[{
  RowBox[{
   SubscriptBox["\[Psi]", "1"], "=", 
   RowBox[{
    RowBox[{
     FractionBox["1", 
      SqrtBox["10"]], 
     RowBox[{"EncodeShells", "[", "\"\<3d0 3dm 3dx 4s 4s~ 4dx\>\"", "]"}]}], 
    "-", 
    RowBox[{
     FractionBox["1", 
      SqrtBox["10"]], 
     RowBox[{"EncodeShells", "[", "\"\<3d0 3dm 3dy 4s 4s~ 4dy\>\"", "]"}]}], 
    "-", 
    RowBox[{
     FractionBox["1", 
      SqrtBox["10"]], 
     RowBox[{"EncodeShells", "[", "\"\<3d0 3dz 3dx 4s 4s~ 4dy\>\"", "]"}]}], 
    "-", 
    RowBox[{
     FractionBox["1", 
      SqrtBox["10"]], 
     RowBox[{"EncodeShells", "[", "\"\<3d0 3dz 3dy 4s 4s~ 4dx\>\"", "]"}]}], 
    "+", 
    RowBox[{
     SqrtBox[
      FractionBox["3", "10"]], 
     RowBox[{"EncodeShells", "[", "\"\<3dz 3dm 3dx 4s 4s~ 4dy\>\"", "]"}]}], 
    "-", 
    RowBox[{
     SqrtBox[
      FractionBox["3", "10"]], 
     RowBox[{
     "EncodeShells", "[", "\"\<3dz 3dm 3dy 4s 4s~ 4dx\>\"", "]"}]}]}]}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{
   SubscriptBox["\[Psi]", "1"], "=", 
   RowBox[{"(", 
    RowBox[{
     RowBox[{"(", 
      RowBox[{"List", "@@", 
       SubscriptBox["\[Psi]", "1"]}], ")"}], "/.", 
     RowBox[{
      RowBox[{"c_", " ", 
       RowBox[{"Ket", "[", "f_", "]"}]}], "\[RuleDelayed]", 
      RowBox[{"{", 
       RowBox[{"c", ",", "f"}], "}"}]}]}], ")"}]}], ";"}]}], "Input"],

Cell[BoxData[{
 RowBox[{
  RowBox[{
   SubscriptBox["\[Psi]", "2"], "=", 
   RowBox[{
    RowBox[{
     FractionBox[
      RowBox[{"-", "1"}], 
      RowBox[{"2", " ", 
       SqrtBox["7"]}]], 
     RowBox[{"EncodeShells", "[", "\"\<3d0 4s 4s~ 4pz 4px 4dy\>\"", "]"}]}], 
    "-", 
    RowBox[{
     FractionBox["1", 
      RowBox[{"2", " ", 
       SqrtBox["7"]}]], 
     RowBox[{"EncodeShells", "[", "\"\<3d0 4s 4s~ 4pz 4py 4dx\>\"", "]"}]}], 
    "+", 
    RowBox[{
     FractionBox["1", 
      SqrtBox["21"]], "2", " ", 
     RowBox[{"EncodeShells", "[", "\"\<3dm 4s 4s~ 4px 4py 4dz\>\"", "]"}]}], 
    "+", 
    RowBox[{
     FractionBox["1", 
      RowBox[{"2", " ", 
       SqrtBox["21"]}]], 
     RowBox[{"EncodeShells", "[", "\"\<3dm 4s 4s~ 4pz 4px 4dy\>\"", "]"}]}], 
    "-", 
    RowBox[{
     FractionBox["1", 
      RowBox[{"2", " ", 
       SqrtBox["21"]}]], 
     RowBox[{"EncodeShells", "[", "\"\<3dm 4s 4s~ 4pz 4py 4dx\>\"", "]"}]}], 
    "+", 
    RowBox[{
     FractionBox["1", 
      RowBox[{"2", " ", 
       SqrtBox["21"]}]], 
     RowBox[{"EncodeShells", "[", "\"\<3dx 4s 4s~ 4px 4py 4dy\>\"", "]"}]}], 
    "+", 
    RowBox[{
     FractionBox["1", 
      SqrtBox["21"]], 
     RowBox[{"EncodeShells", "[", "\"\<3dx 4s 4s~ 4pz 4px 4dz\>\"", "]"}]}], 
    "+", 
    RowBox[{
     FractionBox["1", 
      SqrtBox["7"]], 
     RowBox[{"EncodeShells", "[", "\"\<3dx 4s 4s~ 4pz 4py 4d0\>\"", "]"}]}], 
    "-", 
    RowBox[{
     FractionBox["1", 
      SqrtBox["21"]], 
     RowBox[{"EncodeShells", "[", "\"\<3dx 4s 4s~ 4pz 4py 4dm\>\"", "]"}]}], 
    "-", 
    RowBox[{
     FractionBox["1", 
      RowBox[{"2", " ", 
       SqrtBox["21"]}]], 
     RowBox[{"EncodeShells", "[", "\"\<3dy 4s 4s~ 4px 4py 4dx\>\"", "]"}]}], 
    "+", 
    RowBox[{
     FractionBox["1", 
      SqrtBox["7"]], 
     RowBox[{"EncodeShells", "[", "\"\<3dy 4s 4s~ 4pz 4px 4d0\>\"", "]"}]}], 
    "+", 
    RowBox[{
     FractionBox["1", 
      SqrtBox["21"]], 
     RowBox[{"EncodeShells", "[", "\"\<3dy 4s 4s~ 4pz 4px 4dm\>\"", "]"}]}], 
    "+", 
    RowBox[{
     FractionBox["1", 
      SqrtBox["21"]], 
     RowBox[{"EncodeShells", "[", "\"\<3dy 4s 4s~ 4pz 4py 4dz\>\"", "]"}]}], 
    "-", 
    RowBox[{
     FractionBox["2", 
      SqrtBox["21"]], " ", 
     RowBox[{"EncodeShells", "[", "\"\<3dz 4s 4s~ 4px 4py 4dm\>\"", "]"}]}], 
    "+", 
    RowBox[{
     FractionBox["1", 
      RowBox[{"2", " ", 
       SqrtBox["21"]}]], 
     RowBox[{"EncodeShells", "[", "\"\<3dz 4s 4s~ 4pz 4px 4dx\>\"", "]"}]}], 
    "+", 
    RowBox[{
     FractionBox["1", 
      RowBox[{"2", " ", 
       SqrtBox["21"]}]], 
     RowBox[{
     "EncodeShells", "[", "\"\<3dz 4s 4s~ 4pz 4py 4dy\>\"", "]"}]}]}]}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{
   SubscriptBox["\[Psi]", "2"], "=", 
   RowBox[{"(", 
    RowBox[{
     RowBox[{"(", 
      RowBox[{"List", "@@", 
       SubscriptBox["\[Psi]", "2"]}], ")"}], "/.", 
     RowBox[{
      RowBox[{"c_", " ", 
       RowBox[{"Ket", "[", "f_", "]"}]}], "\[RuleDelayed]", 
      RowBox[{"{", 
       RowBox[{"c", ",", "f"}], "}"}]}]}], ")"}]}], ";"}]}], "Input"]
}, Open  ]],

Cell[CellGroupData[{

Cell[TextData[{
 "Calculate the 2-body RDM of ",
 Cell[BoxData[
  FormBox[
   RowBox[{
    RowBox[{
     RowBox[{"\[VerticalSeparator]", 
      SubscriptBox["\[Psi]", "1"]}], "\[RightAngleBracket]"}], 
    RowBox[{"\[LeftAngleBracket]", 
     RowBox[{
      SubscriptBox["\[Psi]", "2"], "\[VerticalSeparator]"}]}]}], 
   TraditionalForm]], "None"]
}], "Subsection"],

Cell[BoxData[
 RowBox[{
  RowBox[{
   SubscriptBox["\[Rho]", "2"], "=", 
   RowBox[{"RDM", "[", 
    RowBox[{
     SubscriptBox["\[Psi]", "1"], ",", 
     SubscriptBox["\[Psi]", "2"], ",", "2"}], "]"}]}], ";"}]], "Input"],

Cell[BoxData[
 RowBox[{
  RowBox[{"SlaterToPureString", "[", 
   RowBox[{"orbnames_", ",", "s_"}], "]"}], ":=", 
  RowBox[{"StringJoin", "@@", 
   RowBox[{
    RowBox[{"Flatten", "[", 
     RowBox[{
      RowBox[{
       RowBox[{"{", 
        RowBox[{"#", ",", "\"\< \>\""}], "}"}], "&"}], "/@", 
      RowBox[{"orbnames", "\[LeftDoubleBracket]", 
       RowBox[{"FermiToCoords", "[", "s", "]"}], "\[RightDoubleBracket]"}]}], 
     "]"}], "\[LeftDoubleBracket]", 
    RowBox[{"1", ";;", 
     RowBox[{"-", "2"}]}], "\[RightDoubleBracket]"}]}]}]], "Input"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"display", " ", "RDM"}], " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   SubscriptBox["\[Rho]", "2"], "/.", 
   RowBox[{
    RowBox[{"IXI", "[", 
     RowBox[{"ket_", ",", "bra_"}], "]"}], "\[RuleDelayed]", 
    RowBox[{"\"\<\[VerticalSeparator]\>\"", "<>", 
     RowBox[{"SlaterToPureString", "[", 
      RowBox[{
       SubscriptBox["shells", "all"], ",", "ket"}], "]"}], "<>", 
     "\"\<\[RightAngleBracket]\[LeftAngleBracket]\>\"", "<>", 
     RowBox[{"SlaterToPureString", "[", 
      RowBox[{
       SubscriptBox["shells", "all"], ",", "bra"}], "]"}], "<>", 
     "\"\<\[VerticalSeparator]\>\""}]}]}]}]], "Input"],

Cell[BoxData[
 RowBox[{
  FractionBox["\<\"\[VerticalSeparator]3d0 3dx\[RightAngleBracket]\
\[LeftAngleBracket]4pz 4py\[VerticalSeparator]\"\>", 
   SqrtBox["210"]], "+", 
  FractionBox["\<\"\[VerticalSeparator]3d0 3dy\[RightAngleBracket]\
\[LeftAngleBracket]4pz 4px\[VerticalSeparator]\"\>", 
   SqrtBox["210"]], "+", 
  FractionBox["\<\"\[VerticalSeparator]3dz 3dm\[RightAngleBracket]\
\[LeftAngleBracket]4px 4py\[VerticalSeparator]\"\>", 
   SqrtBox["70"]]}]], "Output"]
}, Open  ]]
}, Open  ]],

Cell[CellGroupData[{

Cell["Trace out spin", "Subsection"],

Cell[BoxData[
 RowBox[{
  RowBox[{
   SubscriptBox["\[Rho]s", "2"], "=", 
   RowBox[{"SpinTraceCoulomb", "[", 
    RowBox[{
     RowBox[{"(", 
      RowBox[{"List", "@@", 
       SubscriptBox["\[Rho]", "2"]}], ")"}], "/.", 
     RowBox[{
      RowBox[{"c_", " ", 
       RowBox[{"IXI", "[", 
        RowBox[{"ket_", ",", "bra_"}], "]"}]}], "\[RuleDelayed]", 
      RowBox[{"{", 
       RowBox[{"c", ",", 
        RowBox[{"{", 
         RowBox[{"ket", ",", "bra"}], "}"}]}], "}"}]}]}], "]"}]}], 
  ";"}]], "Input"],

Cell[BoxData[
 RowBox[{
  RowBox[{"BosonToString", "[", 
   RowBox[{"orbnames_", ",", "s_"}], "]"}], ":=", 
  RowBox[{"StringJoin", "@@", 
   RowBox[{
    RowBox[{"Flatten", "[", 
     RowBox[{
      RowBox[{
       RowBox[{"{", 
        RowBox[{"#", ",", "\"\< \>\""}], "}"}], "&"}], "/@", 
      RowBox[{"orbnames", "\[LeftDoubleBracket]", 
       RowBox[{"BosonToCoords", "[", "s", "]"}], "\[RightDoubleBracket]"}]}], 
     "]"}], "\[LeftDoubleBracket]", 
    RowBox[{"1", ";;", 
     RowBox[{"-", "2"}]}], "\[RightDoubleBracket]"}]}]}]], "Input"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"display", " ", "Coulomb", " ", "integrals"}], " ", "*)"}], 
  "\[IndentingNewLine]", 
  RowBox[{
   SubscriptBox["\[Rho]s", "2"], "/.", 
   RowBox[{
    RowBox[{"CoulInt", "[", 
     RowBox[{"a_", ",", "b_"}], "]"}], "\[RuleDelayed]", 
    RowBox[{"\"\<(\>\"", "<>", 
     RowBox[{"BosonToString", "[", 
      RowBox[{
       SubscriptBox["shells", "spatial"], ",", "a"}], "]"}], "<>", 
     "\"\<\[VerticalSeparator]\>\"", "<>", 
     RowBox[{"BosonToString", "[", 
      RowBox[{
       SubscriptBox["shells", "spatial"], ",", "b"}], "]"}], "<>", 
     "\"\<)\>\""}]}]}]}]], "Input"],

Cell[BoxData[
 RowBox[{
  FractionBox["\<\"(3d0 4pz\[VerticalSeparator]3dx 4py)\"\>", 
   SqrtBox["210"]], "+", 
  FractionBox["\<\"(3d0 4pz\[VerticalSeparator]3dy 4px)\"\>", 
   SqrtBox["210"]], "-", 
  FractionBox["\<\"(3dm 4px\[VerticalSeparator]3dz 4py)\"\>", 
   SqrtBox["70"]], "-", 
  FractionBox["\<\"(3dx 4pz\[VerticalSeparator]3d0 4py)\"\>", 
   SqrtBox["210"]], "-", 
  FractionBox["\<\"(3dy 4pz\[VerticalSeparator]3d0 4px)\"\>", 
   SqrtBox["210"]], "+", 
  FractionBox["\<\"(3dz 4px\[VerticalSeparator]3dm 4py)\"\>", 
   SqrtBox["70"]]}]], "Output"]
}, Open  ]]
}, Open  ]],

Cell[CellGroupData[{

Cell["Clean up", "Subsection"],

Cell[BoxData[
 RowBox[{
  RowBox[{"Uninstall", "[", "fermifabML", "]"}], ";"}]], "Input"]
}, Open  ]]
}, Open  ]]
},
WindowSize->{1440, 926},
WindowMargins->{{Automatic, 197}, {Automatic, 51}},
PrivateNotebookOptions->{"VersionedStylesheet"->{"Default.nb"[8.] -> False}},
FrontEndVersion->"10.0 for Microsoft Windows (64-bit) (July 1, 2014)",
StyleDefinitions->"Default.nb"
]
(* End of Notebook Content *)

(* Internal cache information *)
(*CellTagsOutline
CellTagsIndex->{}
*)
(*CellTagsIndex
CellTagsIndex->{}
*)
(*NotebookFileOutline
Notebook[{
Cell[CellGroupData[{
Cell[1486, 35, 52, 0, 63, "Section"],
Cell[1541, 37, 168, 3, 30, "Text"],
Cell[1712, 42, 272, 7, 68, "Text"],
Cell[1987, 51, 399, 8, 68, "Text"],
Cell[2389, 61, 568, 14, 106, "Text"],
Cell[2960, 77, 411, 12, 72, "Input"],
Cell[3374, 91, 1253, 42, 72, "Input"],
Cell[CellGroupData[{
Cell[4652, 137, 2814, 85, 152, "Input"],
Cell[7469, 224, 661, 11, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[8167, 240, 98, 2, 43, "Subsection"],
Cell[8268, 244, 890, 25, 52, "Input"],
Cell[9161, 271, 1380, 48, 128, "Input"],
Cell[10544, 321, 3107, 107, 290, "Input"]
}, Open  ]],
Cell[CellGroupData[{
Cell[13688, 433, 365, 12, 46, "Subsection"],
Cell[14056, 447, 221, 7, 31, "Input"],
Cell[14280, 456, 555, 15, 31, "Input"],
Cell[CellGroupData[{
Cell[14860, 475, 680, 17, 52, "Input"],
Cell[15543, 494, 473, 10, 51, "Output"]
}, Open  ]]
}, Open  ]],
Cell[CellGroupData[{
Cell[16065, 510, 36, 0, 43, "Subsection"],
Cell[16104, 512, 513, 17, 31, "Input"],
Cell[16620, 531, 550, 15, 31, "Input"],
Cell[CellGroupData[{
Cell[17195, 550, 643, 18, 52, "Input"],
Cell[17841, 570, 562, 13, 51, "Output"]
}, Open  ]]
}, Open  ]],
Cell[CellGroupData[{
Cell[18452, 589, 30, 0, 43, "Subsection"],
Cell[18485, 591, 89, 2, 31, "Input"]
}, Open  ]]
}, Open  ]]
}
]
*)

(* End of internal cache information *)

(* NotebookSignature rxDvziwffAsopBwVJsiGWNWS *)
