#!/usr/bin/python
# -*- coding: utf-8 -*-

import os, re
import pandas as pd

pd.set_option('display.max_colwidth', -1)
pd.options.mode.chained_assignment = None  # default='warn'


class Slicer(object):

    def __init__(self, mapfile):

        self.descmap = self.makeTranscript121Mapping(mapfile)


    def split_up_down_results(self, filepath, outdir):

        """
        Parses baySeq result files produced by into 
        Down and Up-regulated lists of transcripts.
        NOTE: Works only for two expression groups
        """

        name = os.path.splitext(os.path.basename(filepath))[0]
        raw = pd.read_csv(filepath, sep="\t")


        sliceDOWN = raw.loc[raw['ordering'] == '1>2']
        listDOWN = sliceDOWN['annotation'].tolist()
        outDown = os.path.join(outdir, name+'_down.txt')
        print outDown
        with open(os.path.join(outdir, name+'_down.txt'), 'w') as fh:
            for d in listDOWN:
                fh.write(d+"\n")


        sliceUP = raw.loc[raw['ordering'] == '2>1']
        listUP = sliceUP['annotation'].tolist()
        outUp = os.path.join(outdir, name+'_up.txt')
        print outUp
        with open(outUp, 'w') as fh:
            for u in listUP:
                fh.write(u+"\n")



    # use for simple 2 column 1-to-1 (key:value) annotation files
    def makeTranscript121Mapping(self, mappingFile):

        names = {}

        FH = open(mappingFile, "r")
        for line in FH:
            res = line.split("\t")
            names[res[0].strip()] = res[1].strip()
        FH.close()

        return names



    def writeoutAdditionalAnnotation(self, inpath, outpath):

        """
        Creates new copy of file with added annotations provided
        via mapping dictionary.
        """

        with open(inpath, "r") as fh:
            lines = fh.readlines()

        out = open(outpath, "w")
        
        for line in lines:
            frags = line.split("\t")
            try:
                out.write(line.strip('\n')+"\t"+
                            self.descmap[frags[0].strip()]+"\n")
            except:
                out.write(line.strip('\n')+"\t\n")
        out.close()





    def read_annotation_file(self, directory, f):

        """
        Reads and parses a tsv annotation file into "two columns" split:
        identifier(key):annotations(values). -- all available annotations.
        Note: Expects whitespace, comma or semi-colon as separator even on
              unannotated entries.
        """

        annotDict = {}

        patt = re.compile('([^,;\s]+)[ ,;\t\v]+(.*)', re.IGNORECASE)

        with open(os.path.join(directory,f), 'r') as fh:
            for line in fh:

                m = patt.search(line)
                if m:
                    annotDict.setdefault(m.group(1), []).append(m.group(2).strip())

        return annotDict






    # WARNING: Currently relies on a strict folder structure and file name pattern
    #          if no enrichment results are being shown probably that structure is
    #          not being met (assuming enrichment data is actually available).
    def process_enrichment_values(self, directory, basename, alpha):

        """
        Processes KEGG/GO enrichment results files (produced by GOstats and topGO)
        into a dictionary.
        """
        

        try:
            bp = self.read_tsv(os.path.join(directory, 'goenrich/BP', 
                                  basename+'_enrichment.tsv'), "GO Biological Process")
        except:
            bp = None
            
            
        try:
            mf = self.read_tsv(os.path.join(directory, 'goenrich/MF', 
                                  basename+'_enrichment.tsv'), "GO Molecular Function")
        except:
            mf = None
            
        try:
            cc = self.read_tsv(os.path.join(directory, 'goenrich/CC',
                                  basename+'_enrichment.tsv'), "GO Cellular Component")
        except:
            cc = None
            
            
        try:
            kegg = self.read_tsv(os.path.join(directory, 'keggenrich', 
                                    basename+'_enrichment.tsv'), "KEGG Pathways")
        except:
            kegg = None


        extra = {'bp':bp, 'mf':mf, 'cc': cc, 'kegg':kegg, 'alpha': alpha}

        return extra


    def read_tsv(self, tsvpath, ont):

        """
        Reads enrichment results tsv files (as produced either by GOstats or
        the topGO R packages) into pandas dataframes (including the headers).
        """

        try:
            if ont != "KEGG Pathways":
                df = pd.read_csv(tsvpath, sep="\t")
            else:
                df = pd.read_csv(tsvpath, sep="\t", converters={'KEGGID': lambda x: str(x)})
        except:
            f = os.path.splitext(os.path.basename(tsvpath))[0]
            print "No %s enrichment found for %s!" % (ont, f)
            raise

        return df