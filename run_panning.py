#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys, os, argparse, shutil, codecs, json
from configobj import ConfigObj
from datetime import datetime as dt

import hoarder, oboe, misc, enricher


def set_or_default(basedir, var, default):

    if var:
        out = var
    else:
        out = os.path.join(basedir, default)
        print out
    create_dir(out)

    return out

def create_dir(dirpath):

    # if not existant already, create folder to save enrichment results
    try: 
        os.makedirs(dirpath)
    except OSError:
        if not os.path.isdir(dirpath):
            raise


def list_dir(directory, path=False):

    files = [f for f in os.listdir(directory) 
                     if os.path.isfile(os.path.join(directory, f))]

    if path:
        files = [os.path.join(directory, x) for x in files]

    return files

# -----------------------------------------------------------------------------

def annotate(inputdir, outdir):

    path = set_or_default(inputdir, outdir, 'annotated')
    files = list_dir(inputdir)

    for f in files:
        infile = os.path.join(inputdir, f)
        outfile = os.path.join(path, 
                           os.path.splitext(f)[0]+".tsv")

        M.writeoutAdditionalAnnotation(infile, outfile)

    return path


def enrich(inputdir, gmap, alpha):

    targets = list_dir(inputdir, True)

    BP = enricher.GOrich(gmap, "BP", alpha)
    MF = enricher.GOrich(gmap, "MF", alpha)
    CC = enricher.GOrich(gmap, "CC", alpha)

    for t in targets:
        BP.perform_go_enrichment(t)
        MF.perform_go_enrichment(t)
        CC.perform_go_enrichment(t)


def dotsvg(inputdir, outdir, alpha):

    path = set_or_default(inputdir, outdir, 'svg')
    X = oboe.OBOe(obopath)

    # target orthologonal ontologies
    aspects = ["BP", "MF", "CC"]

    for a in aspects:
        targetdir = os.path.join(inputdir, "goenrich", a)
        targets = os.listdir(targetdir)

        for t in targets:
            if t.endswith(".tsv"):
                target = os.path.join(inputdir, "goenrich", a, t)
                name = os.path.splitext(target)[0]+"_"+a

            enrich_df = X.read_enrichment_tsv(target)
            goes = X.extract_go_terms(enrich_df, alpha)

            if goes != False:
                try:
                    ntx = X.create_ntx_graph(goes)
                    nodes = X.ntx_nodes(ntx) 
                    edges, root = X.get_directed_edges(nodes)
                    dag = X.create_basic_dag(edges)

                    dot = X.generate_basic_dot(dag)
                    dotname = X.enhance_dot(name, nodes, dot, enrich_df, 
                                            alpha, root)

                    X.export_dag_as(dotname, "svg")
                    # # # # # # # # # # # # # # # # # # # # # # # # # # # #
                    shutil.copy(os.path.splitext(dotname)[0]+".svg", path)
                except:
                    print target
                    print goes
                    print "Something failed with the graph making!"

    return path

# =============================================================================

if __name__ == "__main__":


    parser = argparse.ArgumentParser(conflict_handler='resolve')
    subparsers = parser.add_subparsers(dest='command')


    parser.add_argument('-c', '--config',
                               help="""Indicate an INI file with the variables
                                       relevant for the current GO reports""", 
                               action='store', required=True)
    parser.add_argument('-i', '--inputdir',  help='Indicate the workfolder',
                         action='store', required=True)



    parser_A = subparsers.add_parser('SPLIT', help='''Splits baySeq results
                    files into up and down differentially expressed subset 
                    files. Applicable only for 2-group/pairwise datasets.''')
    parser_A.add_argument('-o', '--outdir', action='store')


    parser_B = subparsers.add_parser('ANNOT', help='''Annotates gene/transcript
                    lists with descriptors found on a user-supplied source file.
                                                                ''')
    parser_B.add_argument('-o', '--outdir', action='store')


    parser_C = subparsers.add_parser('ENRICH', help='''Runs input genes lists on
                    topGO R package and retrieves the enrichment results.''')



    parser_D = subparsers.add_parser('DAG', help='''Produces GO annotation DAGs
                    svg/dot files from topGO enrichment result files.''')
    parser_D.add_argument('-o', '--outdir', action='store')


    parser_E = subparsers.add_parser('REPORT', help='''Runs the pipeline from
        input gene/transcript list up to the generation of a final GO enrichment
        html report file per list.''')
    parser_E.add_argument('-o', '--outdir', action='store')



    args = parser.parse_args()

    # makes dictionary of dictionaries
    config = ConfigObj(args.config)

    now = dt.now()
    config['meta']['date'] = now.strftime("%Y-%m-%d")
    config['meta']['time'] = now.strftime("%H:%M:%S")
    config.write()

    alpha = float(config['vars']['alpha'])
    org = config['vars']['organism']



    basedir = args.inputdir


    g_map = config['sources']['g_map']
    functional_desc_file = config['sources']['functionalDesc']
    obopath = config['sources']['obofile']


# =============================================================================


    # Segregate (up|down) ids from baySeq's 2-group experiments results files 
    # -------------------------------------------------------------------------
    if args.command == 'SPLIT':
        print "Spliting into Up- and Down- regulated DE transcripts..."

        M = misc.Slicer(functional_desc_file)
        path = set_or_default(basedir, args.outdir, 'splits')
        files = list_dir(basedir)

        for f in files:
            M.split_up_down_results(os.path.join(basedir, f), path)


    # Add functional descriptors to gene lists
    # -------------------------------------------------------------------------
    if args.command == 'ANNOT':
        print "Annotating gene lists..."

        M = misc.Slicer(functional_desc_file)
        path = annotate(basedir, args.outdir)


    # Perform GO enrichment (using topGO w/ fisherElim algorithm)
    # -------------------------------------------------------------------------
    if args.command == 'ENRICH':
        print "Computing GO term enrichment..."

        enrich(basedir, g_map, alpha)


    # Generate dot/svg GO graphs (using Graphviz package)
    # -------------------------------------------------------------------------
    if args.command == 'DAG':
        print "Generating GO graphs..."

        path = dotsvg(basedir, args.outdir, alpha)



    if args.command == 'REPORT':

        M = misc.Slicer(functional_desc_file)

        print "annotating..."
        path = annotate(basedir, args.outdir)

        print "computing GO term enrichment..."
        enrich(path, g_map, alpha)

        print "generating GO graphs..."

        reportsout = os.path.join(path, 'reports')
        create_dir(reportsout)
        svgout = os.path.join(reportsout, 'svg')
        create_dir(svgout)

        svgpath = dotsvg(path, svgout, alpha)


        print "assembling final reports..."

        files = list_dir(path)

        jannots = {}
        for f in files:
            annotDict = M.read_annotation_file(path, f)
            jannots["genes"] = json.dumps(annotDict.keys())

            name = os.path.splitext(f)[0]
            scaffold = hoarder.Templater(name, annotDict, org)

            # gets results from topGO/GOstats (filtered by p-value)
            enriched = M.process_enrichment_values(path, 
                                            os.path.splitext(f)[0], alpha)
            # adds "genes" list to the dictionary as new key value
            enriched.update(jannots)

            html = scaffold.render_main_page(enriched)

            # Write out the file...
            filename = os.path.join(reportsout, name+".html")
            with codecs.open(filename, encoding='utf-8', mode="w") as f:
                f.write(html)