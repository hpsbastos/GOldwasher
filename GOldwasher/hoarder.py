#!/usr/bin/python
# -*- coding: utf-8 -*-

import pandas as pd
import re, json, pkg_resources

from jinja2 import Template
from decimal import *

pd.set_option('display.max_colwidth', -1)
pd.options.mode.chained_assignment = None  # default='warn'


class Templater(object):

    def __init__(self, title, annot, org, g2i, links=None):

        self.title = title
        self.annots = annot

        self.svgpaths = {'bp': 'svg/'+title+'_enrichment_BP.svg', 
                         'mf': 'svg/'+title+'_enrichment_MF.svg', 
                         'cc': 'svg/'+title+'_enrichment_CC.svg' }

        if links == None:
            inserts = json.loads(pkg_resources.resource_string('GOldwasher', 
                                                            '/organisms.json'))

        else:
            with open(links, "r") as f:
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

        self.go2ids = self.read_annot_pseudo_json(g2i)



    # NOTE: currently this method overrides "process_enrichment_dict()" 
    #       regarding the display of missing (None) enrichment info.
    def render_main_page(self, plus = None):

        """
        Renders the main page template for a plot and association expression
        information: gene/transcript list, annotations, (GO/KEGG) enrichment
        results.     
        """

        css_resources = '''
          <link href="support/custom.css" rel="stylesheet" 
                                          type="text/css" />
          <link href="support/tab-content/tabcontent.css" 
                                          rel="stylesheet" type="text/css" />
          <link href="support/jquery/jquery-ui.css" 
                                          rel="stylesheet" type="text/css"/>
        '''


        js_resources = '''
          <!-- These following 2 files need to be generated per project -->
          <script type="text/javascript" src="support/annotations.js">
          </script>
          <script type="text/javascript" src="support/descriptions.js">
          </script>


          <!-- Static support (libraries) files -->
          <script src="support/sorttable.js" type="text/javascript"></script>

          <script src="support/tab-content/tabcontent.js"
          type="text/javascript"></script>

          <script  src="support/svg-pan-zoom/svg-pan-zoom.js" 
          type="text/javascript"></script>

          <script src="support/jquery/jquery-1.10.2.js" type="text/javascript">
          </script>
          <script src="support/jquery/jquery-ui.js" type="text/javascript" >
          </script>

           '''

    # -------------------------------------------------------------------------

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
            # TODO: test this/alternatives
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
                    var obs = new MutationObserver(function(mutations, 
                                                            observer) {
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

                    var output = '<table style="border:none; width:100%;">';
                    for (t in dialogcontent){
                        var descript = get_desc(dialogcontent[t]);
                        output += '''+self.link1+'''}
                    output += '</table>'
                    $div.html(output);
                    $("body").append($div);

                    }); // closes "open_dialog" click event

            }); // closes the function...


            function displayDialog(targetDiv) {

              var name = document.getElementById( targetDiv ).name
              var dialog = $( document.getElementById( 'target' + 
                                                        targetDiv ) ).dialog({
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
        u_bp, t_bp, u_mf, t_mf, u_cc, t_cc, t_kegg = self.process_enrichment_dict(bp, mf, cc,
                                                                kegg, alpha)


        if bp is not None or mf is not None or cc is not None:
            goheader = '''<p style="font-size:20px; font-weight: bold;">
                                                     GO term enrichment</p>'''
        else:
            goheader = ""



        # BIOLOGICAL PROCESS table
        # =====================================================================

        if bp is not None:
            bp_info = '''
            <span style="font-size:16px; font-weight: bold;">
            Biological Process ({{ uniq_bp }})
            </span>

            <a id="displayBP" href="javascript:toggle('BP');">hide</a></br>
            <div id="toggleBP" style="width:1200px; display:block;">

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
        # =====================================================================

        if mf is not None:
            mf_info = '''
            <span class="secthead">Molecular Function ({{ uniq_mf }})</span>

            <a id="displayMF" href="javascript:toggle('MF');">show</a></br>
            <div id="toggleMF" class="collapsedDiv">

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
        # =====================================================================

        if cc is not None:
            cc_info = '''
            <span class="secthead">Cellular Component ({{ uniq_cc }})</span>

            <a id="displayCC" href="javascript:toggle('CC');">show</a></br>
            <div id="toggleCC" class="collapsedDiv">

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
        # =====================================================================

        if kegg is not None:
            kegg_info = '''
            <span class="secthead">KEGG pathways enrichment</span>
            {% block table2 %}
            {{ tablekegg }}
            {% endblock %}
            <br/>
            '''
        else:
            kegg_info = ""

    # -------------------------------------------------------------------------

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
                               uniq_bp = u_bp,
                               tablegomf = t_mf,
                               uniq_mf = u_mf,
                               tablegocc = t_cc,
                               uniq_cc = u_cc,
                               tablekegg = t_kegg,
                               annots = self.render_gene_table(),
                               )
        return html



    def process_title(self):

        """
        NOTE: currently a passthrough placeholder!
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

            <span class="secthead">
                Genes/transcripts ({{ len1 }}/{{ total }})</span>

        <a id="displayList" href="javascript:toggle('List');">show</a>

        <div id="toggleList" class="collapsedDiv">
            <table class="genedesc">
            <thead>
                <tr>
                    <th scope="col">Identifier</th>
                    <th scope="col">Descriptor</th>
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

        <div id="toggleBastards" class="collapsedDiv">
        {% for id in extra %}
            '''+self.link2+'''
            {% if not loop.last %}
                , 
            {% endif %}
        {% endfor %}
        </div>

        ''')


        # CHEEKY HACK
        # ---------------------------------------------------------------------
        naughty = []
        stripped = {}

        for key in self.annots:

            if self.annots[key][0].strip() == "":
                naughty.append(key)
            else:
                stripped.setdefault(key,[]).append(self.annots[key][0].strip())

        lenstrp = len(stripped)
        lennaut = len(naughty)
        # ---------------------------------------------------------------------



        if len(self.annots) > 0:
            table = tablegene.render(genedata=stripped, 
                                     extra=naughty, 
                                     len1=str(lenstrp), 
                                     len2=lennaut, 
                                     total=str(lenstrp+lennaut))
        else:
            table = Template("").render()
            print "uh oh!"
            print "Something went wrong generating the gene list table!"

        return table



    def go_slice(self, enrichRes, alpha):

        """
        Slices up and post-processes individual GO ontology 
        enrichment results for displaying.
        """

        ontslice = enrichRes.loc[enrichRes['elimFisher'] < alpha]

        targs = ontslice["GO.ID"].tolist()
        ug = len(self.get_unique_genes(targs))

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

        return ug, restable



    def get_unique_genes(self, go_ids):

        """
        For a given enriched GO term list 
        return the _unique_ & _significant_ 
        genes annotated to those terms.
        """

        l = []
        for go in go_ids:
            for x in self.go2ids[go]:
                if x in self.annots.keys():
                    l.append(x)

        return list(set(l))


    def read_annot_pseudo_json(self, fpath):

        """
        Reads the created pseudo-json file
        and returns it as a dictionary
        """

        with open(fpath, 'r') as fh:
            data = fh.readlines()

        data = [x.strip() for x in data[1:]]
        outstr = ' '.join(data)
        outstr = '{'+ outstr
        td = json.loads(outstr)

        return td


    # TODO: Evaluate if displaying "No significant enriched terms found!"
    #       is of worth, or just leave blank... **kwargs...it?
    def process_enrichment_dict(self, gobp, gomf, gocc, kegg, alpha):

        """
        Method receives dataframes containing GO/KEGG enrichment results and
        converts them into html tables.
        """

        if gobp is not None:
            ubp, gobptable = self.go_slice(gobp, alpha)
        else:
            gobptable = self.not_found_response()
            ubp = None


        if gomf is not None:
            umf, gomftable = self.go_slice(gomf, alpha)
        else:
            gomftable = self.not_found_response()
            umf = None

        if gocc is not None:
            ucc, gocctable = self.go_slice(gocc, alpha)
        else:
            gocctable = self.not_found_response()
            ucc = None

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

        return ubp, gobptable, umf, gomftable, ucc, gocctable, keggtable


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
        Necessary trick to be able to sort exponential
        numbers with the sorttable library.
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