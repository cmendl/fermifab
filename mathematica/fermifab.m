(* ::Package:: *)

(* :Name: FermiFab *)

(* :Author: Christian B. Mendl
	http://christian.mendl.net *)

(* :Copyright:
	Copyright (C) 2008-2015, Christian B. Mendl
	All rights reserved.
	This program is free software; you can redistribute it and/or
	modify it under the terms of the Simplified BSD License
	http://www.opensource.org/licenses/bsd-license.php *)

(* :Summary:
	FermiFab stands for 'fermion laboratory' and is a particle physics toolbox
	for small numbers of fermions. It is intended for numerical as well as analytical
	calculations and centered around Slater determinants and the index mapping
	for creation and annihilation operators, i.e., the antisymmetry of fermions
	is inherently incorporated into the toolbox. *)

(* :References:
	Christian B. Mendl
	The FermiFab toolbox for fermionic many-particle quantum systems
	Comput. Phys. Commun. 182, 1327-1337 (2011)
	preprint http://arxiv.org/abs/1103.0872 *)

BeginPackage["FermiFab`"]


(* Utility functions *)

PurgeSparseArray[A_SparseArray]:=SparseArray[Cases[ArrayRules[A],(i_->x_):>(i->x)/;x!=0],Dimensions[A]]

Reshape[A_,dim_]:=Fold[Partition,Flatten[A],Reverse[Rest[dim]]]

(* PermuteOp[KroneckerProduct[A1,A2,A3],{3,1,2},dim1,dim2] == KroneckerProduct[A3,A1,A2] *)
PermuteOp[A_,p_,dim1_,dim2_]:=Reshape[Transpose[Reshape[A,Join[dim1,dim2]],Join[Ordering[p],Length[dim1]+Ordering[p]]],{Times@@dim1,Times@@dim2}]

SimNullSpace[A__]:=NullSpace[Total[ConjugateTranspose[#].#&/@{A}]]

BlockDiagonalMatrix[A__]:=Module[{T=IdentityMatrix[Total[Length/@{A}]],n=1,i},
	If[Or@@(Head[#]===SparseArray&/@{A}),T=SparseArray[T]]; (* return a sparse matrix if any of the inputs is sparse *)
	(i=n;;n+Length[#]-1;T[[i,i]]=#;n+=Length[#])&/@{A};T]

SymToSparseArray[expr_Plus,sym_,map_]:=SymToSparseArray[List@@expr,sym,map]
SymToSparseArray[expr_,sym_,map_]:=SparseArray[Flatten[{expr/.{0->{}}}]
	/.{x_ sym[a_,b_]:>Rule[Flatten[{Position[map,a],Position[map,b]}],x]}
	/.{sym[a_,b_]:>Rule[Flatten[{Position[map,a],Position[map,b]}],1]},Length[map]{1,1}]


(* Basic definitions *)

(* Fermions *)

FermiMapPerms[orbs_Integer,N_Integer]:=Permutations[Join[ConstantArray[1,N],ConstantArray[0,orbs-N]]]
FermiMap[orbs_Integer,N_Integer]:=Sort[FromDigits[#,2]&/@FermiMapPerms[orbs,N]]
FermiMap[orbs_List,N_List]:=Sort[FromDigits[#,2]&/@Flatten[Outer@@Join[{Flatten[Reverse[{##}]]&},Apply[FermiMapPerms,Transpose[{orbs,N}],2],{1}],Length[orbs]-1]]

InterleaveFermi[\[Psi]_List,fm_List]:=Select[Transpose[{\[Psi],fm}],First[#]!=0&]
InterleaveFermi[\[Psi]_SparseArray,fm_List]:=Cases[ArrayRules[\[Psi]],Rule[{i_},a_]:>{a,fm[[i]]}/;a!=0]

MergeFermiConfig[orbs_List,N_List]:=
	Module[{fmsrc=FermiMap[orbs,N],fmtrg=FermiMap[Total[orbs],Total[N]]},
		SparseArray[Flatten[{Position[fmtrg,fmsrc[[#]]],#}]->1&/@Range[Length[fmsrc]],{Length[fmtrg],Length[fmsrc]}]]

FermiToCoords[f_]:=Flatten[Position[Reverse[IntegerDigits[f,2]],1]]
CoordsToFermi[x_]:=Total[BitSet[0,x-1]]

SlaterToString[orbnames_,s_]:="\[VerticalSeparator]"<>StringJoin@@Flatten[{#," "}&/@orbnames[[FermiToCoords[s]]]][[1;;-2]]<>"\[RightAngleBracket]"
DisplayFermi[orbnames_,\[Psi]_]:=Plus@@(First[#]SlaterToString[orbnames,#[[2]]]&/@\[Psi])

(* Bosons *)

(* '0' serves as mode delimiter, e.g., bit pattern '111001101' means one particle in 1st mode, two in 2nd mode, zero in 3rd mode, and three in 4th mode *)
BosonMap[modes_Integer,N_Integer]:=FermiMap[modes+N-1,N]

BosonToCoords[b_]:=Module[{x},x=FermiToCoords[b];x-=Range[0,Length[x]-1]]
CoordsToBoson[x_]:=CoordsToFermi[Sort[x]+Range[0,Length[x]-1]]


(* Fermi RDM related functions *)

AnnihilSign[n_,a_]:=If[BitAnd[n,a]==a,1-2Mod[Total[DigitCount[BitAnd[2^#-1,n-a],2,1]&/@(FermiToCoords[a]-1)],2],0]
CreateSign[n_,c_]:=If[BitAnd[n,c]==0,1-2Mod[Total[DigitCount[BitAnd[2^#-1,n],2,1]&/@(FermiToCoords[c]-1)],2],0]

Begin["`Private`"]

(* reference implementation - use 'MLGenerateRDM' for speed *)
GenerateRDM[orbs_,p1_,N1_,N2_]:=Module[{p2=N1-N2+p1,fmp1,fmp2,fmN1,fmN2},
	fmp1=FermiMap[orbs,p1];
	fmp2=FermiMap[orbs,p2];
	fmN1=FermiMap[orbs,N1];
	fmN2=FermiMap[orbs,N2];
	SparseArray[Apply[Flatten[Join[{#1,#2},Position[fmN1,fmN2[[#3]]-fmp1[[#1]]+fmp2[[#2]]],{#3}]]->
		AnnihilSign[fmN2[[#3]],fmp1[[#1]]]CreateSign[fmN2[[#3]]-fmp1[[#1]],fmp2[[#2]]]&,
	Select[Flatten[Outer[{##}&,Range[Length[fmp1]],Range[Length[fmp2]],Range[Length[fmN2]]],2],
		BitAnd[fmN2[[#3]],fmp1[[#1]]]==fmp1[[#1]]\[And]BitAnd[fmN2[[#3]]-fmp1[[#1]],fmp2[[#2]]]==0&@@#&],2]]]
GenerateRDM[orbs_,p_,N_]:=GenerateRDM[orbs,p,N,N]

End[] (* end `Private` section *)

p2N[h_,orbs_,p_,N_]:=Module[{K=MLGenerateRDM[orbs,p,N]},Sum[h[[i,j]]K[[j,i]],{i,Dimensions[h][[1]]},{j,Dimensions[h][[2]]}]]

TensorOp[A_,fm_List]:=Module[{Afm},
	Afm=Outer[Det[A[[FermiToCoords[#1],FermiToCoords[#2]]]]&,fm,fm];
	If[Head[A]===SparseArray,Afm=SparseArray[Afm]];Afm]
TensorOp[A_,N_Integer]:=TensorOp[A,FermiMap[Length[A],N]]

GenerateConfig[nc_Integer,N_Integer]:=BinCounts[BosonToCoords[#],{1,nc+1}]&/@BosonMap[nc,N]
GenerateConfig[nc_Integer,N_Integer,maxC_List]:=Select[GenerateConfig[nc,N],Function[c,Select[maxC-c,#<0&]=={}]]

(* workaround for MLPutInteger64 switching integer signs *)
TwoComplement[n_]:=FromDigits[IntegerDigits[n,2]/.{1->0,0->1},2]+1/;n<0
TwoComplement[n_]:=n/;n>=0

RDM[\[Psi]_,p_Integer]:=RDM[\[Psi],\[Psi],p]
RDM[\[Psi]_,\[Chi]_,p_Integer]:=Expand[Total[Flatten[Outer[Conjugate[First[#1]]First[#2]MLSlaterRDM[Last[#2],Last[#1],p]&,\[Chi],\[Psi],1]]]]/.
	{IXI[a_,b_]:>IXI[TwoComplement[a],TwoComplement[b]]}
(* 'Expand' sums coefficients of same bra-ket states *)


(* Trace out spin *)

SpinTraceSingleBase[i_,j_]:=If[Mod[i-j,2]==0,SingInt@@(CoordsToFermi[{#}]&/@Ceiling[{i,j}/2]),0]

SpinTraceSingle[\[Rho]_]:=Expand[Total[First[#]SpinTraceSingleBase@@Flatten[FermiToCoords/@Last[#]]&/@\[Rho]]]
SpinTraceCoulomb[\[Rho]_]:=Expand[Total[First[#]MLSpinTraceCoulombBase@@Flatten[FermiToCoords/@Last[#]]&/@\[Rho]]]


EndPackage[]
