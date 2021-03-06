#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr

"""
Classes to interface in a limited way (presets) with 
R's packages 'topGO' and 'GOstats' via rpy2.
"""

#TODO: Re-write these classes without using the 'robject.globalenv'
#      because that leads to issues, when "crossing" objects from 
#      different methods. Add some input validation.

# ==============================================================================

class KEGGer(object):


    def __init__(self, keggmap, alpha, organism):

        self.mappings = keggmap
        self.alpha = alpha
        self.organism = organism

        KEGGdb = importr('KEGG.db')
        GOstats = importr('GOstats')
        GSEABase = importr('GSEABase')

    def perform_kegg_enrichment(self, targetspath):

        """
        Uses the R GOstats package to calculate KEGG pathway enrichment for a
        suplied list of genes/transcripts of interest, saving results to file.
        """

        robjects.globalenv["mappings"] = self.mappings
        robjects.globalenv["alpha"] = self.alpha
        robjects.globalenv["organism"] = self.organism


        read_target_group_of_interest(targetspath)

        handle_result_saving(targetspath, "keggenrich/")

        robjects.r('''tryCatch(
                    {
                      KEGGmap <- read.csv(mappings, sep="\t", colClasses=
                                    c("character", "character"), header=FALSE)
                      colnames(KEGGmap) <- c("path_id", "gene_id")
                      keggFrame=KEGGFrame(KEGGmap,organism=organism)
                      gsc <- GeneSetCollection(keggFrame, setType = 
                                    KEGGCollection())
                      universe <- unique(KEGGmap[,2])
                      universe <- unlist(as.character(universe))

                      kparams <- GSEAKEGGHyperGParams(
                                  name="Custom GSEA based params",
                                  geneSetCollection=gsc,
                                  geneIds = targetset,
                                  universeGeneIds = universe,
                                  pvalueCutoff = alpha,
                                  testDirection = "over")

                      kOver <- hyperGTest(kparams)
                      enrichRes <- summary(kOver)

                      write.table(enrichRes, file = outname, sep = "\t", 
                         row.names = FALSE, quote = FALSE)

                    }, error = function(e){ print( paste(basename, 
                    " failed to have any KEGG pathways mapped to its members!") 
                    ) }
                    )
        ''')


# ==============================================================================

class GOrich(object):

    # NOTE: alpha parameter is not currently been used, here.
    #       Check possibily of using it as a filter later.
    def __init__(self, gomap, ontology='BP', alpha=0.01):

        self.mappings = gomap
        self.ontology = ontology
        self.alpha = alpha

        topGO = importr('topGO')
        JSONlite = importr('jsonlite')

    def perform_go_enrichment(self, targetspath):

        """
        Uses the R topGO package to calculate GO term enrichment in suplied
        lists of genes/transcripts of interest, saving results to file.
        """

        robjects.globalenv["mappings"] = self.mappings
        robjects.globalenv["ontology"] = self.ontology
        robjects.globalenv["alpha"] = self.alpha

        read_target_group_of_interest(targetspath)
        handle_result_saving(targetspath, "goenrich/"+self.ontology)

        robjects.r('''
        capture.output(id2go <- readMappings(file = mappings), file="/dev/null")
        transcriptNames <- names(id2go)

        interestingGenes <- factor(as.integer( 
                        transcriptNames %in% as.character(unlist(targetset)) ), levels=c(0,1) )
        names(interestingGenes) <- transcriptNames

        success = FALSE
        ''')

        robjects.r('''
                      tryCatch(
                      {
                      capture.output(
                        GOdata <- new("topGOdata", ontology = ontology,
                        allGenes = interestingGenes, annot = annFUN.gene2GO, 
                        gene2GO = id2go)
                      , file="/dev/null" )
                      success <- TRUE
                      }, 
                      error = function(e){ 
                                print( paste(basename, 
                                "failed to have any GO", ontology, 
                                "terms mapped to its members!") )
                                success <- FALSE 
                                }
                      )
                   ''')

        robjects.r('''
                if (success == TRUE){
                    capture.output(
                        resultFisher <- runTest(GOdata, algorithm = "classic",
                        statistic = "fisher")
                    , file="/dev/null" )

                    capture.output(
                        resultFisher.Elim <- runTest(GOdata, algorithm = "elim",
                        statistic = "fisher")
                    , file="/dev/null" )

                    enrichRes <- GenTable(GOdata, classicFisher = resultFisher,
                               elimFisher = resultFisher.Elim, 
                               orderBy = "elimFisher", ranksOf ="classicFisher",
                               topNodes = 30)

                    write.table(enrichRes, file = outname, sep = "\t", 
                                row.names = FALSE, quote = FALSE)
                }
                   ''')


    def extract_full_annot_mapping(self, genelist, pathout):

        """
        Extracts list of GO term annotations (ancestors included)
        into a psedo-json (.js w/ 'annotmap' var file.
        """

        robjects.globalenv["mappings"] = self.mappings

        robjects.globalenv["targetset"] = robjects.StrVector(genelist)
        robjects.globalenv["pathout"] = pathout

        robjects.r('''
        capture.output(id2go <- readMappings(file = mappings), file="/dev/null")
        transcriptNames <- names(id2go)

        interestingGenes <- factor(as.integer( 
                 transcriptNames %in% as.character(unlist(targetset)) ), levels=c(0,1) )
        names(interestingGenes) <- transcriptNames

        success = FALSE
        onts <- c("BP", "MF", "CC")
        comboList <- list()
        ''')

        robjects.r('''
                  tryCatch(
                        {
                        for (ont in onts){
                            capture.output(
                                GOdata <- new("topGOdata", ontology = ont,
                                              allGenes = interestingGenes, 
                                              annot = annFUN.gene2GO, 
                                              gene2GO = id2go)
                            , file="/dev/null" )

                            go2ids <- genesInTerm(GOdata, usedGO(GOdata))
                            for (go in names(go2ids)){
                                comboList[[go]] <- go2ids[[go]]
                            }
                        }
                        success <- TRUE
                        }, 

                        error = function(e){
                                  print(e) 
                                  print( paste(basename, 
                                  "failed to have any GO", ont, 
                                  "terms mapped to its members!") )
                                  success <- FALSE 
                                  }
                  )
              ''')

        robjects.r('''
                if (success == TRUE){
                    jsonout = toJSON(comboList, pretty = TRUE)
                    fh <- file(pathout)
                    writeLines(paste("var annotmap = ", jsonout), fh)
                    close(fh)
                }
                   ''')


# accessory py2r helpers
# ------------------------------------------------------------------------------

def read_target_group_of_interest(targetspath):

    """
    Reads a tsv file, extracts the first column (while assuming they
    are transcript/gene identifiers) and sets them as a vector in R
    global environament.
    """

    genes = []
    with open(targetspath) as fh:
        for line in fh:
            token = line.split("\t")
            genes.append(token[0].strip())

    robjects.globalenv["targetset"] = robjects.StrVector(genes)


def handle_result_saving(targetspath, enrich):

    """
    Creates (if necessary) output folder to save enrichment results
    and 'returns' full path to the results file being saved.
    """

    # SAVE IT TO DISK
    basename = os.path.splitext(os.path.basename(targetspath))[0]
    foldername = os.path.dirname(targetspath)

    robjects.globalenv["basename"] = basename

    outname = os.path.join(foldername, enrich, basename+"_enrichment.tsv")
    robjects.globalenv["outname"] = outname

    # if not exists, create folder to save enrichment results
    savepath = os.path.join(foldername, enrich)
    try: 
        os.makedirs(savepath)
    except OSError:
        if not os.path.isdir(savepath):
            raise