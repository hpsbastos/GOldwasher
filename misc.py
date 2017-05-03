#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
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



    def writeoutAdditionalAnnotation(self, directory, out, readFile):

        """
        Creates new copy of file with added annotations provided
        via mapping dictionary.
        """

        OUTname = os.path.join(directory, out, 
                               os.path.splitext(readFile)[0]+".tsv")

        IN = open(os.path.join(directory, readFile), "r")
        lines = IN.readlines()
        IN.close()

        OUT = open(OUTname, "w")
        
        for line in lines:
            frags = line.split("\t")
            try:
                OUT.write(line.strip('\n')+"\t"+
                            self.descmap[frags[0].strip()]+"\n")
            except:
                OUT.write(line.strip('\n')+"\t\n")
        OUT.close()