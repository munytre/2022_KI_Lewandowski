#!/usr/bin/env R
# -*- coding: utf-8 -*-
#
#  example.R
#  
#  Copyright 2020 Oscar C. Bedoya-Reina <oscar@oscar-J53kOiSr>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  

#Method to read clustering and return clusters
rtrnClstrs=function(outClstrs)
	{
	dataset2=read.csv(outClstrs, header=T, sep="	")
	cell=dataset2[,1]
	cell=as.character(cell)
	dataset=data.frame(dataset2[,-1])
	rownames(dataset)=cell
	return(t(dataset))
	}

#Wrapper to execute the code
runBisqueRNAStndAln=function(cnts_inBlk_fl,cnts_SC_inFl,clstrs_SC_inFl, 
	outFlPrprtnsBisqueRNA,sep=',')
	{
	library(BisqueRNA)
	library(Biobase)
	#Load bulk
	cnts_blk_sampleID = rtrnCntsFrmFl(cnts_inBlk_fl)
	cnts_blk = cnts_blk_sampleID$gene
	#Load reference SC
	cnts_ref_sampleID = rtrnCntsFrmFl(cnts_SC_inFl,sep=sep)
	cnts_ref = cnts_ref_sampleID$gene
	ar_clstrsOri = rtrnClstrs(clstrs_SC_inFl)
	smplNames_pData_trngData = rtrnTrngDtStrctr(cnts_SC_inFl,clstrs_SC_inFl)
	smplNames = smplNames_pData_trngData$smplNames
	pData = smplNames_pData_trngData$pData
	trngData = smplNames_pData_trngData$trngData
	#
	bulk.eset <- Biobase::ExpressionSet(assayData = cnts_blk)
	sample.ids <- colnames(cnts_ref[,match(colnames(ar_clstrsOri),
	colnames(cnts_ref))])
	sc.pheno <- data.frame(check.names=F, check.rows=F,stringsAsFactors=F,
	row.names=sample.ids,SubjectName=c(smplNames),cellType=c(ar_clstrsOri))
	#
	sc.meta <- data.frame(labelDescription=c("SubjectName","cellType"),
	row.names=c("SubjectName","cellType"))
	sc.pdata <- new("AnnotatedDataFrame",data=sc.pheno,varMetadata=sc.meta)
	#
	sc.eset <- Biobase::ExpressionSet(assayData=cnts_ref[,
	match(colnames(ar_clstrsOri),colnames(cnts_ref))],phenoData=sc.pdata)
	res <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset,sc.eset, 
	markers=NULL,use.overlap=FALSE)
	#
	ref.based.estimates <- t(res$bulk.props)
	#Write matrices
	write.csv(ref.based.estimates,outFlPrprtnsBisqueRNA)
	}

#Method to format counts
rtrnCntsFrmFl=function(data_input, sep="	")
	{
	dataset2=read.csv(data_input,header=T,sep=sep,check.names=FALSE)
	genename=dataset2[,1]
	genename=as.character(genename)
	dataset=dataset2[,-1]
	rownames(dataset)=genename
	theColNames=colnames(dataset)	
	gene = as.matrix(sapply(dataset,as.numeric))
	colnames(gene)=theColNames
	rownames(gene)=rownames(dataset)
	#
	sampleID = as.character(lapply(strsplit(colnames(gene), split="_[0-9].[a-z]"),head, n=1))
	#
	return(list('gene'=gene,'sampleID'=sampleID))
	}

#Method to data
rtrnTrngDtStrctr=function(cnts_SC_inFl,clstrs_SC_inFl,sep=',')
	{
	#Load reference SC
	cnts_ref_sampleID = rtrnCntsFrmFl(cnts_SC_inFl,sep=sep)
	cnts_ref = cnts_ref_sampleID$gene
	sampleID = cnts_ref_sampleID$sampleID
	ar_clstrsOri = rtrnClstrs(clstrs_SC_inFl)
	#assayData:
	smplNames = as.character(lapply(strsplit(as.character(colnames(ar_clstrsOri)), 
	split="_[0-9].[a-z]"),head, n=1))
	pData = data.frame(row.names=colnames(cnts_ref)[match(colnames(ar_clstrsOri), 
	colnames(cnts_ref))],sampleID=c(smplNames),cellType=c(ar_clstrsOri))
	trngData = ExpressionSet(cnts_ref[,match(colnames(ar_clstrsOri), 
	colnames(cnts_ref))], phenoData=AnnotatedDataFrame(data=pData))
	#
	return(list('smplNames'=smplNames,'pData'=pData,'trngData'=trngData))
	}

#Inputs:
##################
#1) Count in bulk-sequences file
# cnts_inBlk_fl="20220531_Oscars_tool/readcount_frmtd_inGSE127465_GRCm38frmh38.csv"
cnts_inBlk_fl="20220531_Oscars_tool/ALS_Human_EWCE/Data/tt_AH.csv"
#The format for this file is as follows:
#
#!Sample_title	Bulk-sample_1	Bulk-sample_2	...	Bulk-sample_N
#Gene_name_A	Counts_bulk-sample_1	Counts_bulk-sample_2	...	Counts_bulk-sample_N
#Gene_name_B	Counts_bulk-sample_1	Counts_bulk-sample_2	...	Counts_bulk-sample_N
#...
#Gene_name_Z	Counts_bulk-sample_1	Counts_bulk-sample_2	...	Counts_bulk-sample_N
#
#Note:
#In the example file, gene names are ENSEMBL codes but it can be gene symbols
##################

##################
#2) Counts in single cell file
# cnts_SC_inFl="20220531_Oscars_tool/GSE127465_GRCm38frmh38_inreadcount_frmtd.ensmbl.cnts.csv"
cnts_SC_inFl="20220531_Oscars_tool/ALS_Human_EWCE/Data/SS_10cat.csv"
#The format for this file is as follows:
#
#file,Cell_1,Cell_2,...,Cell_N
#Gene_name_A,Counts_cell-sample_1,Counts_cell-sample_2,...,Counts_cell-sample_N
#Gene_name_B,Counts_cell-sample_1,Counts_cell-sample_2,...,Counts_cell-sample_N
#...
#Gene_name_Z,Counts_cell-sample_1,Counts_cell-sample_2,...,Counts_cell-sample_N
#
#Note:
#In the example file, gene names are ENSEMBL codes but it can be gene symbols
##################

##################
#3) Cluster labels
# clstrs_SC_inFl="20220531_Oscars_tool/GSE127465_human.immgen.clstrsLbld.csv.csv"
clstrs_SC_inFl="20220531_Oscars_tool/ALS_Human_EWCE/Data/SS_names.csv"
#The format for this file is as follows:
#
#cell	cluster
#Cell_1	cell_type_name_for_cell_1
#Cell_2	cell_type_name_for_cell_2
#...
#Cell_N	cell_type_name_for_cell_N
##################

#4) Output file for proportion
outFlPrprtnsBisqueRNA="20220531_Oscars_tool/READOUT/READOUT_HUMAN_AH.csv"

#Execute
runBisqueRNAStndAln(cnts_inBlk_fl,cnts_SC_inFl,clstrs_SC_inFl,outFlPrprtnsBisqueRNA)
