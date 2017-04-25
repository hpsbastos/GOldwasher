#!/usr/bin/python
# -*- coding: utf-8 -*-

import networkx as nx
import pygraphviz as pgv
import pandas as pd
import os, math, subprocess, re, tempfile, json

from ontology import OBOOntology
from jinja2 import Template
from decimal import *

pd.set_option('display.max_colwidth', -1)
pd.options.mode.chained_assignment = None  # default='warn'



class Templater(object):

    def __init__(self, title, annot, org):

        self.title = title
        self.annots = annot

        self.svgpaths = {'bp': 'svg/'+title+'_enrichment_BP.svg', 
                         'mf': 'svg/'+title+'_enrichment_MF.svg', 
                         'cc': 'svg/'+title+'_enrichment_CC.svg' }

        with open("organisms.json", "r") as f:
            inserts = json.load(f)

        try:
            self.link1 = inserts[org]['insertlink1']
            self.link2 = inserts[org]['insertlink2']
        except:
            print "Are your target organism's linkouts in 'organisms.json' ?"
            # creates fallback non-linking fake anchors
            self.link1 = """'<a href="#null" alt="'+dialogcontent[t]+'">'
                             +dialogcontent[t]+'</a> '+descript+'</br>';"""

            self.link2 = '''<a name="{{ id }}" href="#null" alt="{{ id }}">
                            {{ id }}</a>'''





    # NOTE: currently this method overrides "process_enrichment_dict()" 
    #       regarding the display of missing (None) enrichment info.
    def render_main_page(self, plus = None):

        """
        Renders the main page template for a plot and association expression
        information: gene/transcript list, annotations, (GO/KEGG) enrichment
        results.     
        """

        #                   dir: ltr;

        # TODO: make css customizing interface/options and/or add as external
        csstables = '''
          <style type="text/css">

            .etables {
                width: 1200px;
            }

            table.sortable 
            th:not(.sorttable_sorted):not(.sorttable_sorted_reverse):not(.sorttable_nosort):after { 
                content: " '''+str('\\25B4\\25BE')+'''" 
            }

            </style>
            '''



        css_resources = '''
          <link href="support/tab-content/template1/tabcontent.css" 
                  rel="stylesheet" type="text/css" />
          <link href="support/jquery-ui.css" rel="stylesheet" type="text/css"/>

            '''+csstables


#       <link href="support/jquery-ui.css" rel="stylesheet" type="text/css"/>




        js_resources = '''
          <!-- These following 2 files need to be generated per project -->
          <script type="text/javascript" src="support/annotations.js">
          </script>
          <script type="text/javascript" src="support/descriptions.js">
          </script>

          <!-- Static support (libraries) files -->
          <script type="text/javascript" src="support/sorttable.js"></script>

          <script type="text/javascript" 
          src="support/tab-content/tabcontent.js"></script>

          <script type="text/javascript" 
          src="support/svg-pan-zoom/svg-pan-zoom.js"></script>

          <script type="text/javascript" src="support/jquery-1.10.2.js"></script>
          <script type="text/javascript" src="support/jquery-ui.js"></script>

           '''

          # TODO: add switch to choose from "local" and "web" version
          # <script src="http://code.jquery.com/jquery-1.10.2.js"></script>
          # <script src="http://code.jquery.com/ui/1.11.4/jquery-ui.js"></script>


    # ----------------------------------------------------------------------------


        if plus is not None:

            bp = plus['bp']
            mf = plus['mf']
            cc = plus['cc']
            kegg = plus['kegg']

            # a small hack to pass along the significance 
            # level (filter for GO enrichments)
            alpha = plus['alpha']

            # TODO: testing...    
            #jannots = plus['annots']
            jgenes = plus['genes']

        else:
            bp = None
            mf = None
            cc = None
            kegg = None
            alpha = 0.05
            # TODO: testing...    
            #jannots = None
            jgenes = None



        # TODO: testing...
        if jgenes is not None:




            bigjson = '''
            <script type="text/javascript">



            function toggle(target) {
                var ele = document.getElementById("toggle"+target);
                var text = document.getElementById("display"+target);
                if(ele.style.display == "block") {
                        ele.style.display = "none";
                    text.innerHTML = "show";
                }
                else {
                    ele.style.display = "block";
                    text.innerHTML = "hide";
                }
            } 



            $(function() {
                $(".dialog").dialog({
                    autoOpen: false
                    });
         
                $(".open_dialog").click(function(e) {
                    e.preventDefault();
                    
                    var linkID = $(this).attr('id')
                    var selector = "target" + linkID;

                    // Define a new observer
                    var obs = new MutationObserver(function(mutations, observer) {
                      // look through all mutations that just occured
                      for(var i=0; i<mutations.length; ++i) {
                        // look through all added nodes of this mutation
                        for(var j=0; j<mutations[i].addedNodes.length; ++j) {
                          // was a child added with the selector ID? 
                          if(mutations[i].addedNodes[j].id == selector) {

                                    displayDialog(linkID);

                          }
                        }
                      }
                    });

                    // have the observer observe the document body 
                    // for changes in children nodes
                    obs.observe($("body").get(0), {
                      childList: true
                    });




                    var $div = $("<div>", {id: selector, class: "dialog", 
                                           style: "display:none;"});

                    var dialogcontent = get_genes(linkID);

                    var output = "";
                    for (t in dialogcontent){
                        var descript = get_desc(dialogcontent[t]);
                        output += '''+self.link1+'''}
                    $div.html(output);
                    $("body").append($div);

                    }); // closes "open_dialog" click event


            }); // closes the function...



            function displayDialog(targetDiv) {

              var name = document.getElementById( targetDiv ).name
              var dialog = $( document.getElementById( 'target' + targetDiv ) ).dialog({
                title: name,
                width: 800
              });
            } 


            function isInArray(value, array) {
              return array.indexOf(value) > -1;
            }


            function get_desc(id){
                var desc = descmap[id]
                return desc
            }


            function get_genes(id){ 
              var slice = annotmap[id];
              var out = [];
              for (g in slice) { 
                if (isInArray(slice[g], genes)){
                  out.push(slice[g]);
                } 
              }
              return out;
            }
            '''+"var genes = "+jgenes+" \n\n "+"</script>"


            js_resources += bigjson



        # pre-process enrichment data for display
        t_bp, t_mf, t_cc, t_kegg = self.process_enrichment_dict(bp, mf, cc, kegg, alpha)


        if bp is not None or mf is not None or cc is not None:
            goheader = '''<p style="font-size:20px; font-weight: bold;">
                                                     GO term enrichment</p>'''
        else:
            goheader = ""



        # BIOLOGICAL PROCESS table
        # ======================================================================

        if bp is not None:
            bp_info = '''
            <span style="font-size:16px; font-weight: bold;">
            Biological Process
            </span>

            <a id="displayBP" href="javascript:toggle('BP');">hide</a></br>
            <div id="toggleBP" style="width:1200px; display: block;">

            {% block table1a %}
            {{ tablegobp }}
            {% endblock %}
            <br/>

            </div>
            </br>
            '''

        else:
            bp_info = ""

        # MOLECULAR FUCNTION table
        # ======================================================================

        if mf is not None:
            mf_info = '''
            <span style="font-size:16px; font-weight: bold;">
            Molecular Function
            </span>

            <a id="displayMF" href="javascript:toggle('MF');">show</a></br>
            <div id="toggleMF" style="width:1200px; display: none;">

            {% block table1b %}
            {{ tablegomf }}
            {% endblock %}
            <br/>

            </div>
            </br>
            '''
        else:
            mf_info = ""


        # CELLULAR COMPONENT table
        # ======================================================================

        if cc is not None:
            cc_info = '''
            <span style="font-size:16px; font-weight: bold;">
            Cellular Component
            </span>

            <a id="displayCC" href="javascript:toggle('CC');">show</a></br>
            <div id="toggleCC" style="width:1200px; display: none;">

            {% block table1c %}
            {{ tablegocc }}
            {% endblock %}
            <br/>
            <br/>

            </div>
            </br>
            '''
        else:
            cc_info = ""


        # KEGG PATHWAYS table
        # ======================================================================

        if kegg is not None:
            kegg_info = '''
            <span style="font-size:20px; font-weight: bold;">
            KEGG pathways enrichment
            </span>
            {% block table2 %}
            {{ tablekegg }}
            {% endblock %}
            <br/>
            '''
        else:
            kegg_info = ""

    # --------------------------------------------------------------------------

        template = Template('''<!DOCTYPE html>
        <html lang="en">
            <head>
                <meta charset="utf-8">
                <title>{{ title }}</title>
                {{ js_resources }}
                {{ css_resources }}
                {{ script }}
            </head>
            <body>
                {{ header }}
                
                {{ svgdiv }}
                <br/>

                <br/>
                '''+goheader+bp_info+mf_info+cc_info+kegg_info+
                ''' 
                {% block table3 %}
                {{ annots }}
                {% endblock %}

            </body>
        </html>
        ''')

        html = template.render(js_resources=js_resources,
                               css_resources=css_resources,
                               title = self.title,
                               header = self.process_title(),
                               svgdiv = self.render_svg_inset(),
                               tablegobp = t_bp,
                               tablegomf = t_mf,
                               tablegocc = t_cc,
                               tablekegg = t_kegg,                           
                               annots = self.render_gene_table(),
                               )
        return html













    def process_title(self):

        """
        NOTE: placeholder!
        TODO: Refine this method...
        """

        return '<h2>'+self.title+'</h2>'




    def render_svg_inset(self):

        """
        Renders the svg inset tabs on the html report files
        """

        svgdiv = Template('''

            <div id="dagtabs" style="width: 1200px;">

                <ul class="tabs">
                    <li><a href="#view1">Biological Process</a></li>
                    <li><a href="#view2">Molecular Function</a></li>
                    <li><a href="#view3">Cellular Component</a></li>
                </ul>
                <div class="tabcontents">

                    <div id="view1">
                        <object id="bpsvg" type="image/svg+xml" 
                         data="{{ svg['bp'] }}" width="100%" height="600px">
                            Your browser does not support SVG
                        </object>
                    </div>

                    <div id="view2">
                        <object id="mfsvg" type="image/svg+xml" 
                         data="{{ svg['mf'] }}" width="100%" height="600px">
                            Your browser does not support SVG
                        </object>
                    </div>

                    <div id="view3">
                        <object id="ccsvg" type="image/svg+xml" 
                         data="{{ svg['cc'] }}" width="100%" height="600px">
                            Your browser does not support SVG
                        </object>
                    </div>

                </div>
            </div>


            <script>

                var sbp = document.getElementById("bpsvg");
                sbp.addEventListener('load', function(){
                    var panZoomArm = svgPanZoom('#bpsvg', {
                                      zoomEnabled: true,
                                      controlIconsEnabled: true
                                     });
                });

                var smf = document.getElementById("mfsvg");
                smf.addEventListener('load', function(){
                    var panZoomArm = svgPanZoom('#mfsvg', {
                                      zoomEnabled: true,
                                      controlIconsEnabled: true
                                     });
                });

                var scc = document.getElementById("ccsvg");
                scc.addEventListener('load', function(){
                    var panZoomArm = svgPanZoom('#ccsvg', {
                                      zoomEnabled: true,
                                      controlIconsEnabled: true
                                     });
                });

            </script>


        ''')

        div = svgdiv.render(svg=self.svgpaths)
        return div




    # make Jinja2 table template for gene/transcript annotation data/list
    # NOTE: Currently all available annotations per gene/transcript are
    # dumped into a single cell. TODO: ponder the best way to improve that!
    # Maybe a sparse matrix-like table via pandas dataframe? 
    # --------------------------------------------------------------------
    def render_gene_table(self):

        """
        Renders the available annotations for the genes/transcripts being
        plotted into 1st column: gene/transcript identifier; 2nd column:
        all available annotations.
        """

        tablegene = Template('''

        <br/>

            <span style="font-size:20px; font-weight: bold;">
                Genes/transcripts ({{ len1 }}/{{ total }})</span>

        <a id="displayList" href="javascript:toggle('List');">show</a>

        <div id="toggleList" style="width:1200px; display: none;">
            <table dir="ltr" width="1200" border="1">
            <thead>
                <tr>
                    <th scope="col">Identifier</th>
                    <th scope="col">Annotations</th>    
                </tr>
            </thead>
            <tbody>
                {% for id in genedata: %}
                    {% if genedata[id][0].strip() !="": %}
                    <tr>
                        <td>
                        '''+self.link2+'''
                        </td>
                        <td>
                        {{ genedata[id][0].strip() }}
                        </td>
                    </tr>
                    {% endif %}                    
                {% endfor %}
                <br/>

            </tbody>
        </table>

        </div>



        <h3> No description ({{ len2 }})</h3>
        <p>
        <a id="displayBastards" href="javascript:toggle('Bastards');">show</a>

        <div id="toggleBastards" style="width:1200px; display: none;">
        {% for id in extra %}
            '''+self.link2+'''
            {% if not loop.last %}
                , 
            {% endif %}
        {% endfor %}
        </div>

        ''')


        # DIRTY HACK

        naughty = []
        stripped = {}

        for key in self.annots:

            if self.annots[key][0].strip() == "":
                naughty.append(key)
            else:
                stripped.setdefault(key,[]).append(self.annots[key][0].strip())

        lenstrp = len(stripped)
        lennaut = len(naughty)

        # DEBUG
        # print len(self.annots)
        # print lenstrp
        # print lennaut


        if len(self.annots) > 0:
            table = tablegene.render(genedata=stripped, 
                                     extra=naughty, 
                                     len1=str(lenstrp), 
                                     len2=lennaut, 
                                     total=str(lenstrp+lennaut))
        else:
            table = Template("").render()
            print "uh oh!"

        return table





    def go_slice(self, enrichRes, alpha):

        """
        Slices up and post-processes individual GO ontology 
        enrichment results for displaying.
        """

        ontslice = enrichRes.loc[enrichRes['elimFisher'] < alpha]
        if len(ontslice) > 0:

            ontslice["Significant"] = ontslice["GO.ID"].map(str) + str("|") + \
                                      ontslice["Term"].map(str) + str("|") + \
                                      ontslice["Significant"].map(str)
            ontslice["Significant"] = ontslice["Significant"].map( 
                lambda x: \
                '''<a id="%s" name="GO: %s" class="open_dialog" 
                href="javascript:void(0)">%s</a> '''% tuple(x.split("|")) )

            ontslice["GO.ID"] = ontslice["GO.ID"].map(
                lambda x: '<a href="http://amigo.geneontology.org/amigo/term/'
                +str(x)+'" target="_blank">'+str(x)+'</a>')

            cols = ['GO.ID', 'Term', 'Annotated', 'Significant', 'elimFisher']
            restable = ontslice[cols].to_html(index=False, 
                                              classes="etables sortable", 
                                              escape=False)
            restable = self.process_go_headers(restable)

            restable = self.process_hidden_keys(restable)

        else:
            restable = self.not_found_response()

        return restable


    # TODO: Evaluate if displaying "No significant enriched terms found!"
    #       is of worth, or just leave blank... **kwargs...it?
    def process_enrichment_dict(self, gobp, gomf, gocc, kegg, alpha):

        """
        Method receives dataframes containing GO/KEGG enrichment results and
        converts them into html tables.
        """

        if gobp is not None:
            gobptable = self.go_slice(gobp, alpha)
        else:
            gobptable = self.not_found_response()


        if gomf is not None:
            gomftable = self.go_slice(gomf, alpha)
        else:
            gomftable = self.not_found_response()


        if gocc is not None:
            gocctable = self.go_slice(gocc, alpha)
        else:
            gocctable = self.not_found_response()


        if kegg is not None:
            
            kegg["Count"] =  kegg["KEGGID"].map(str) + str("|") + \
                             kegg["Term"].map(str) + str("|") + \
                             kegg["Count"].map(str) 
            kegg["Count"] = kegg["Count"].map( 
                lambda x: \
                '''<a id="%s" name="KEGG: %s" class="open_dialog" 
                    href="javascript:void(0)">%s</a> '''% tuple(x.split("|")) )

            kegg["KEGGID"] = kegg["KEGGID"].map(
                lambda x: \
                '<a href="http://www.genome.jp/dbget-bin/www_bget?pathway:map'
                +str(x)+'" target="_blank">'+str(x)+'</a>')

            # changes the column order of the KEGG enrichment results 
            cols = ['KEGGID', 'Term', 'Size', 'Count', 'Pvalue']
            keggtable = kegg[cols].to_html(index=False, 
                                           classes="etables sortable", 
                                           escape=False)

            keggtable = keggtable.replace('<th>KEGGID</th>', 
                                    '<th class="sorttable_nosort">KEGGID</th>')

        else:
            keggtable = self.not_found_response()

        return gobptable, gomftable, gocctable, keggtable


    # render this line when no GO terms found to be enriched 
    def not_found_response(self):

        tablego = Template('<br/>No significant enriched terms found!<br/>')
        table = tablego.render()
        return table

    def process_go_headers(self, table):

        out = table.replace('<th>GO.ID</th>', 
                            '<th class="sorttable_nosort">GO.ID</th>')
        return out



    def process_hidden_keys(self, table):

        """
        Necessary trick to be able to
        sort exponential numbers.
        """

        getcontext().prec = 100
        patt = re.compile("<td>([0-9-e\.]+)</td>\s+</tr>")

        for m in re.finditer(patt, table):
            d = Decimal(m.group(1))
            s = '{0:f}'.format(d)
            replacement = '<td sorttable_customkey="'+s+'">'+ \
                                    m.group(1)+'</td>\n</tr>\n'
            table = table.replace(m.group(0), replacement)

        return table


# =============================================================================


class Maestro(object):

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
            #print res[0].strip()
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



    def read_obo(self):

        ont = OBOOntology(self.obofile)
        return ont

    def create_ntx_graph(self, termlist):

        gx = self.ontology.to_networkx(termlist)
        return gx

    def ntx_nodes(self, ntxgraph):

        return ntxgraph.nodes();



    def get_directed_edges(self, termlist):

        """
        Takes full list of terms comprising a GO
        sub-graph and return its directed edges,
        and root term (the only sink).
        """

        root = ""
        edgelist = []

        clonelist = [x for x in termlist]
        terms = []

        for term in clonelist:
            par_edges = self.ontology.parent_edges(term)
            if len(par_edges) > 0:      # if it has any parent...
                for pe in par_edges:
                    #print term, pe
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

        df = pd.read_csv(filepath, sep="\t")
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

                    term_name = self.ontology[row["GO.ID"]].name
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
                    term_name = self.ontology[go].name
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
        #       or make a few presets available
        max_colour = 240
        min_colour = 170

        if pvalue == 0:
            hue = max_colour/float(360)
        else:
            hue = ( min_colour + abs(math.log(pvalue)) ) / float(360)

        return (hue, 1, 1)