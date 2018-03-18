#!/usr/bin/python
# -*- coding: utf-8 -*-

import networkx as nx
import pandas as pd
import os, math, subprocess, tempfile, misc

from libraries.ontology import OBOOntology

pd.set_option('display.max_colwidth', -1)
pd.options.mode.chained_assignment = None  # default='warn'

class OBOe(object):

    """
    This class is used to create SVG GO sub-graphs
    with overlaid color-coded enrichment statistics
    and/or PNG/MAP files.

    INPUT:

    - Ontology OBO file
    - topGO result files

    OUTPUT:

    - dot, svg (and optionally other graphic filetypes)

    """


    def __init__(self, obopath):

        self.obofile = obopath
        self.ontology = self.read_obo()
        self.ontology.alt2id = self.add_alt2id()


    def read_obo(self):

        ont = OBOOntology(self.obofile)
        return ont

    def add_alt2id(self):

        "Needed to handle deprecated go ids"

        alt2id = {}
        for term in self.ontology.id2term:   
            for i in range(0, len(self.ontology.id2term[term].tags())-1 ):
                tag = self.ontology.id2term[term]._format_single_tag(i)
                token = tag.split(': ')
                if token[0] == 'alt_id':
                    alt2id.setdefault(token[1].strip(), None)
                    alt2id[token[1].strip()] = term
        return alt2id


    def create_ntx_graph(self, termlist):

        gx = self.ontology.to_networkx(termlist)
        return gx

    def ntx_nodes(self, ntxgraph):

        return ntxgraph.nodes();


    def get_directed_edges(self, termlist):

        """
        Takes full list of terms comprising a GO
        sub-graph and return its directed 'is_a'
        edges and the root term (the only sink).
        """

        root = ""
        edgelist = []

        clonelist = [x for x in termlist]
        terms = []

        for term in clonelist:
            par_edges = self.ontology.parent_edges(term)
            # if term has any parent...
            if len(par_edges) > 0:
                for pe in par_edges:
                    if pe[0] == "is_a":
                        edgelist.append( ( term, pe[1].id) )
            else:
                root = term

        return edgelist, root


    def create_basic_dag(self, edgelist):

        dag = nx.DiGraph()
        dag.add_edges_from(edgelist)

        return dag


    def generate_basic_dot(self, dag):

        """
        Converts networkx graph into a dot
        file contents output (as a list).
        """

        P = nx.nx_pydot.to_pydot(dag)
        A = nx.nx_agraph.to_agraph(dag)

        temp = tempfile.NamedTemporaryFile()
        A.write(temp.name)

        dot = []
        for line in temp:
            dot.append(line.strip())

        return dot


    def read_enrichment_tsv(self, filepath):

        df = pd.read_csv(filepath, sep="\t",
                        converters={'elimFisher': 
                            lambda x: misc.parse_ev(x)})
        df[['elimFisher']] = df[['elimFisher']].apply(pd.to_numeric)
        return df


    def extract_go_terms(self, df, alpha):

        """
        Tries to extract a list of GO accessions ids at a given
        statistical significance from a topGO results dataframe.
        """

        try:
            fatia = df["GO.ID"].loc[df['elimFisher'] < float(alpha)]
            subset = [x for x in fatia]
            return subset
        except:
            return False

    # -------------------------------------------------------------------------


    def enhance_dot(self, outname, nodes, listdotedges, df, alpha, root):

        #temp = tempfile.NamedTemporaryFile(delete=False)

        outfile = outname+".dot"
        temp = open(outfile, "w")

        sig = self.extract_go_terms(df, alpha)

        # dot file header
        temp.write('strict digraph  {\n graph[rankdir="BT"];\n')

        added = [] # stores list of nodes that are customized in the dot

        for index, row in df.iterrows():
            if row["elimFisher"] < alpha:
                # set root term as sink
                if row["GO.ID"] == root:
                    # get term name here    
                    temp.write('"'+row["GO.ID"]+'" '+'[rank="sink"];\n')
                    added.append(row["GO.ID"])
                else:
                    try:
                        term_name = self.ontology[row["GO.ID"]].name
                    except KeyError:
                        term_name = self.ontology[self.ontology.alt2id[row["GO.ID"]]].name

                    truncterm = self.process_term_name(term_name)
                    trunctermplus = row["GO.ID"]+"\n"+truncterm+"\n"+ \
                                    str(row["elimFisher"])+"\n"+"("+ \
                                    str(row["Significant"])+"/"+ \
                                    str(row["Annotated"])+")"

                    color_hsv = self.get_enrichment_color(row["elimFisher"])

                    tooltip = term_name
                    temp.write('"'+row["GO.ID"]+'" [ id="'+row["GO.ID"]+
                               '" shape="note" style="filled" color="'+
                               str(color_hsv[0])+' '+str(color_hsv[1])+' '+
                               str(color_hsv[2])+'" label="'+trunctermplus+
                                                '" tooltip="'+tooltip+'"];\n')

                    added.append(row["GO.ID"])


        # Handle statistically non-significant nodes
        for go in nodes:
            if go not in added:
                if go == root:
                    # get term name here
                    term_name = self.ontology[go].name
                    temp.write('"'+go+'" '+'[rank="sink" style="filled" label="'
                                                    +go+'\\n'+term_name+'"];\n')
                else:
                    try:
                        term_name = self.ontology[go].name
                    except KeyError:
                        term_name = self.ontology[self.ontology.alt2id[go]].name
                    temp.write('"'+go+'" '+'[style="filled" label="'+go+
                                           '" tooltip="'+term_name+'"];\n')

        for edge in listdotedges[1:]:
            temp.write(edge+"\n")
        temp.close()

        return temp.name


    def export_dag_as(self, dotfile, filetype):

        """
        Exports a dag (dot) file as either an SVG or PNG/MAP
        """

        outdef = ()

        if filetype == "svg":
            outdef = (".svg", "-Tsvg")
        elif filetype == "png":
            outdef = (".png", "-Tpng")
        elif filetype == "map":
            outdef = (".map", "-Tcmapx")
        else:
            print "Filetype ", filetype, " not recognised!"
            sys.exit()

        outfile = os.path.splitext(dotfile)[0]+outdef[0]
        with open(outfile,'w') as s_out:
            subprocess.call(['dot', outdef[1], dotfile], stdout=s_out)


    def process_term_name(self, name):

        """
        Processes a 'term name' in order to obtain a 
        truncated version to better fit a graphviz node.
        """

        if len(name) > 20:

            trunc_name = ''
            token = name.split(' ')

            # if 1st word has MORE than 20 characters
            if len(token[0]) > 20:
                # RETURN the first 18 characters plus '...'
                return token[0][0:17] + '...'

            else:
                last = ""
                for e in token:
                    # if ANY subsequent SINGLE word(s) has 
                    # MORE than 20 characters...
                    if len(e) > 20:
                        diff = 20 - (len(last) + len(e))
                        return trunc_name + last + e[0:diff] + '...'

                    else:
                        # if ANY subsequent SINGLE LINE has
                        # MORE than 20 characters...
                        if (len(last) + len(e)) < 20:

                            last = last + " " + e
                            if e == token[-1]:
                                return trunc_name + last

                        else:
                            trunc_name += last + '\\n'
                            last = e
                            if e == token[-1]:
                                return trunc_name + last

                return trunc_name

        # OTHERWISE just return it... 
        else:
            return name


    def get_enrichment_color(self, pvalue):

        """
        Creates an hue for a p-value
        within a chosen HSV color range.
        """

        # hard-coded range
        # TODO: improve this, possibly make it customizable
        #       or make a few sensible presets available
        # pall_1 = (170, 240)
        max_colour = 240
        min_colour = 170

        if pvalue == 0:
            hue = max_colour/float(360)
        else:
            hue = ( min_colour + abs(math.log(pvalue)) ) / float(360)

        return (hue, 1, 1)