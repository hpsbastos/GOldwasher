#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys, os, argparse, shutil, codecs, json
from configobj import ConfigObj
from datetime import datetime as dt
import GOgalore as gg
from makimono import enricher
from makimono import toolbox


def create_dir(dirpath):

    # if not exists, create folder to save enrichment results
    try: 
        os.makedirs(dirpath)
    except OSError:
        if not os.path.isdir(dirpath):
            raise

# =============================================================================

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument('-c', '--config',
                               help="""Indicate an INI file with the variables
                                       relevant for the current GO reports""", 
                               action='store', required=True)
    parser.add_argument('-i', '--inputfolder',  help='Indicate the workfolder',
                         action='store', required=True)

    args = parser.parse_args()

    # makes dictionary of dictionaries
    config = ConfigObj(args.config)

    now = dt.now()
    config['meta']['date'] = now.strftime("%Y-%m-%d")
    config['meta']['time'] = now.strftime("%H:%M:%S")
    config.write()

    alpha = float(config['vars']['alpha'])
    org = config['vars']['organism']

    #basedir = os.path.expanduser( os.path.join('~', config['vars']['basedir']) )
    basedir = args.inputfolder

    splitsuffix = config['vars']['splitsuffix']
    enrichsuffix = config['vars']['enrichsuffix']
    reportsuffix = config['vars']['savereports']

    g_map = os.path.expanduser( os.path.join('~', config['vars']['g_map']) )
    annotDir = os.path.expanduser( 
                        os.path.join('~', config['vars']['annotmapbasedir']) ) 
    descriptionFile = config['vars']['descriptionFile']


    obopath = os.path.expanduser( os.path.join('~', config['vars']['obopath']) )


# ===================================


    splitslocation = os.path.join(basedir, splitsuffix)
    annotedlocation = os.path.join(splitslocation, enrichsuffix)



    M = gg.Maestro(os.path.join(annotDir, descriptionFile))


    # Do the (up|down) splits
    # ------------------------------------------------------------------------
#     print "Spliting into Up- and Down- regulated DE transcripts..."

#     listing = os.listdir(basedir)
#     files = [x for x in listing if os.path.isfile(os.path.join(basedir, x))]

#     splitslocation = os.path.join(basedir, splitsuffix)
#     create_dir(splitslocation)

#     for f in files:
# #        print f
#         M.split_up_down_results(os.path.join(basedir, f), splitslocation)



    splitslocation = basedir

    # Add functional descriptors
    # -------------------------------------------------------------------------
    print "Annotating split lists..."

    annotedlocation = os.path.join(splitslocation, enrichsuffix)
    create_dir(annotedlocation)


    files = [f for f in os.listdir(splitslocation) 
             if os.path.isfile(os.path.join(splitslocation, f))]

    for f in files:
        M.writeoutAdditionalAnnotation(splitslocation, enrichsuffix, f)

    # Perform GO enrichment
    # -------------------------------------------------------------------------
    print "Computing GO term enrichment..."

    targets = [os.path.join(annotedlocation, x) 
               for x in os.listdir(annotedlocation) if 
               os.path.isfile(os.path.join(annotedlocation, x))]

    BP = enricher.GOrich(g_map, "BP", alpha)
    MF = enricher.GOrich(g_map, "MF", alpha)
    CC = enricher.GOrich(g_map, "CC", alpha)

    for t in targets:
        BP.perform_go_enrichment(t)
        MF.perform_go_enrichment(t)
        CC.perform_go_enrichment(t)

    # -------------------------------------------------------------------------

    reportsout = os.path.join(basedir, reportsuffix)
    create_dir(reportsout)
    svgout = os.path.join(reportsout, "svg")
    create_dir(svgout)



    # Generate GO graphs
    # -------------------------------------------------------------------------
    print "Generating GO graphs..."


    X = gg.OBOe(obopath)

    # target orthologonal ontologies
    aspects = ["BP", "MF", "CC"]

    # -------------------------------------------------------------------------

    for a in aspects:
        targetdir = os.path.join(annotedlocation, "goenrich", a)
        targets = os.listdir(targetdir)
        for t in targets:
            if t.endswith(".tsv"):
                target = os.path.join(annotedlocation, "goenrich", a, t)
                name = os.path.splitext(target)[0]+"_"+a

    # -------------------------------------------------------------------------

            enrich_df = X.read_enrichment_tsv(target)
            goes = X.extract_go_terms(enrich_df, alpha)

            if goes != False:
                try:
                    ntx = X.create_ntx_graph(goes)
                    nodes = X.ntx_nodes(ntx) 
                    edges, root = X.get_directed_edges(nodes)
                    dag = X.create_basic_dag(edges)

                    dot = X.generate_basic_dot(dag)
                    dotname = X.enhance_dot(name, nodes, dot, enrich_df, alpha,
                                            root)

                    X.export_dag_as(dotname, "svg")
                    shutil.copy(os.path.splitext(dotname)[0]+".svg", svgout)
                except:
                    print target
                    print goes
                    print "Something failed with the graph making!"



    print "Generating final reports..."
    # ------------------------------------------------------------------------

    listing = os.listdir(annotedlocation)
    files = [x for x in listing if 
                    os.path.isfile(os.path.join(annotedlocation, x))]

    jannots = {}
    for f in files:
        annotDict = toolbox.read_annotation_file(annotedlocation, f)
        jannots["genes"] = json.dumps(annotDict.keys())

        name = os.path.splitext(f)[0]
        scaffold = gg.Templater(name, annotDict, org)

        # gets results from topGO/GOstats (filtered by p-value)
        enriched = toolbox.process_enrichment_values(annotedlocation, 
                                        os.path.splitext(f)[0], alpha)
        # adds "genes" list to the dictionary as new key value
        enriched.update(jannots)

        html = scaffold.render_main_page(enriched)

        # Write out the file...
        filename = os.path.join(reportsout, name+".html")
        with codecs.open(filename, encoding='utf-8', mode="w") as f:
            f.write(html)