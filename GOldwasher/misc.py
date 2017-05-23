#!/usr/bin/python
# -*- coding: utf-8 -*-

import os, re, json
import pandas as pd

pd.set_option('display.max_colwidth', -1)
pd.options.mode.chained_assignment = None  # default='warn'


class Slicer(object):

    def __init__(self, mapfile):

        self.descmap = self.read_annotation_file(mapfile, False)


    # useful for simple 2 column 1-to-1 (key:value) annotation files
    def read_annotation_file(self, f, list_=True):

        '''
        Reads and parses a tsv annotation file into 'TWO columns' split:
        identifier(key):annotations(values). -- all available annotations.
        Note: Expects _whitespace_, _comma_ or _semi-colon_ as separator.
        '''

        # TODO: increment/add input sanitization

        annotDict = {}
        patt = re.compile('([^,;\s]+)[,;\t\v]+(.*)', re.IGNORECASE)

        try:
            with open(f, 'r') as fh:
                for line in fh:
                    m = patt.search(line)
                    if m:
                        annotDict.setdefault(m.group(1), 
                                                []).append(m.group(2).strip())
                    else:
                        # considers that row has an identifier only
                        annotDict.setdefault(line.strip(), []).append('')
        except IOError:
            print 'Could not read file:', f
            raise

        if list_ == False:
            for key in annotDict:
                annotDict[key] = ' '.join(annotDict[key])

        return annotDict


    def writeout_additional_annotation(self, inpath, outpath):

        '''
        Creates new copy of file with added annotations provided
        via mapping dictionary.
        '''

        try:
            with open(inpath, 'r') as fh:
                lines = fh.readlines()

            with open(outpath, 'w') as out:
                for line in lines:
                    frags = line.split('\t')
                    try:
                        #print self.descmap[frags[0]]
                        out.write(line.strip('\n')+'\t'+
                                    self.descmap[frags[0].strip()]+'\n')
                    except KeyError:
                        out.write(line.strip('\n')+'\t\n')
                        raise
        except IOError as err:
            print 'Could not read/write file:', err.filename
            raise


    # NOTE: maybe make it more 'path flexible'
    def process_enrichment_values(self, directory, basename, alpha):

        '''
        Processes KEGG/GO enrichment results files (produced by 
        GOstats and topGO) into a single dictionary.
        '''

        try:
            bp = self.read_enrichment_tsv(os.path.join(directory, 'goenrich/BP', 
                                  basename+'_enrichment.tsv'), 
                                  'GO Biological Process')
        except:
            bp = None
            
            
        try:
            mf = self.read_enrichment_tsv(os.path.join(directory, 'goenrich/MF', 
                        basename+'_enrichment.tsv'), 'GO Molecular Function')
        except:
            mf = None
            
        try:
            cc = self.read_enrichment_tsv(os.path.join(directory, 'goenrich/CC',
                        basename+'_enrichment.tsv'), 'GO Cellular Component')
        except:
            cc = None
            
            
        try:
            kegg = self.read_enrichment_tsv(os.path.join(directory,'keggenrich', 
                                basename+'_enrichment.tsv'), 'KEGG Pathways')
        except:
            kegg = None


        extra = {'bp':bp, 'mf':mf, 'cc': cc, 'kegg':kegg, 'alpha': alpha}

        return extra


    def read_enrichment_tsv(self, tsvpath, ont):

        '''
        Reads enrichment results tsv files (as produced either by GOstats or
        the topGO R packages) into pandas dataframes (including the headers).
        '''

        try:
            if ont != 'KEGG Pathways':
                df = pd.read_csv(tsvpath, sep='\t',
                    converters={'elimFisher': lambda x: parse_ev(x)})
                df[['elimFisher']] = df[['elimFisher']].apply(pd.to_numeric)
            else:
                df = pd.read_csv(tsvpath, sep='\t', 
                                    converters={'KEGGID': lambda x: str(x)})
        except:
            f = os.path.splitext(os.path.basename(tsvpath))[0]
            print 'No %s enrichment found for %s!' % (ont, f)
            raise

        return df


    def generate_support_js_file(self, annots, outdir, f, varname):

        '''
        Turns a given dictionary and creates mock-json files
        (js files with a single var) to be used by the html 
        enrichments reports.
        '''
        outfile = os.path.join(outdir, f)

        if os.path.isfile(outfile):
            return
        else:
            try:
                with open(outfile, 'w') as fh:
                    fh.write('var '+varname+' = ')
                    jsondump = json.dumps(annots)
                    fh.write(jsondump)
            except IOError:
                print 'Could not read file:', f
                raise

        return outfile







def parse_ev(val):

    '''
    Identifies and converts '< 1e-30' strings
    into numeric convertable strings
    '''

    patt = re.compile('^<\s+(.+)')
    m = patt.search(val)
    if m:
        return m.group(1)
    else:
        return val