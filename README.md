

GOldwasher
=======

GOldwasher is a light wrapper for the R package 
[topGO](https://bioconductor.org/packages/release/bioc/html/topGO.html). The 
function of GOldwasher is limited to the calculation of GO term enrichment 
(via elimFisher algorithm) of target gene lists and also the generation their 
respective GO subgraph images (via 
[Graphviz](http://www.graphviz.org/) - Graph Visualization Software) 
subsuming their annotations. These elements are then incorporated together into 
html a single report file per input list which can then be  interactively 
explored. This module purpose and focus is to facilitate batch processing of 
several gene lists.

----------

**Current release:** *0.2.1 (alpha)*

Provides basic functionality and basic documentation. Methods and functions 
are working if the input does not deviate from expected (none to little input 
sanitization) but they haven't been extensively tested.


**TODO:**
  
- Write (proper) tests
- Improve documentation
- Re-think R interface mod
- Make OBOe mod independent of ontology mod
- Extend functionality (select and compare subsets)




Installation
------------

    pip install GOldwasher-0.2.1.tar.gz

or

    pip install GOldwasher

or you can just checkout the code out of here and fiddle with it as you wish.

**Requirements**:
The non-python dependencies are:

- **graphviz** http://www.graphviz.org/
- **topGO** http://bioconductor.org/packages/release/bioc/html/topGO.html (and inherently **R**)
- **GO OBO file** (http://geneontology.org/page/download-ontology)


These packages/file should be installed/downloaded for GOldwasher to work. All 
other python dependencies should be automatically resolved by pip.




Usage
-----

GOldwasher can be used as a module, making use of its methods, or more easily
it can conveniently be used from its command tool 'goldpanner'


    goldpanner [-h] -c CONFIG -i INPUTDIR
                  {SPLIT,ANNOT,ENRICH,DAG,REPORT}
 e.g.:

    goldpanner -c settings.ini -i lists/ REPORT

-**c**  
ini file with general settings using the following structure:


    [vars]
    alpha = 0.01  
    organism = phaeodactylum

    [sources]
    functionalDesc = /path/to/tabseparedfile/withIDtabFunctionalDescription.txt
    g_map = /path/to/mappings/identifier2GOaccessions.txt
    obofile = /path/to/go-basic.obo

    #linkinsets = /path/to/custom/organisms.json

**[vars]**
---

**alpha** - significance level  

**organism** - name of the organism (as key name on 'organisms.json')

**[sources]**
---

**functionalDesc** 

    path to tab-separated file holding a column of identifiers and their matching functional descriptors. 

**g_map** 

    path to tab-separated file holding a column of identifiers and a second column with their associated GO term accession numbers separated by commas.            

    e.g.:
        Phatr3_J43587.t1  GO:0006396,GO:0005622,GO:0005515 

**obofile** 

    path to the GO ontology obo file. It can be downloaded from: http://purl.obolibrary.org/obo/go/go-basic.obo


**linkinsets**

    If using organisms other than 'Arabidopis thaliana' or 'Phaeodactylum tricornutum' uncomment this variable and set it as the path to the customized 'organisms.json'. By default no cross-links are generated for unknown/unset organisms.


**-i** directory with the target lists.


  **Commands**:
  ---

    SPLIT - splits baySeq result files (pairwise conditions) into up and down-regulated lists.

    ANNOT - annotates identifiers lists with respective available functional descriptors.

    ENRICH - performs GO term enrichment on the annotated lists.

    DAG - generates color-coded GO graph image (svg format) from (topGO) enrichment results.

    REPORT - generates an interactive html GO enrichment report for each list on the input directory.  



**optional argument**:

**-o** output directory (can be used with all commands except ENRICH)







3rd party libraries
-------

Additional required 3rd party content is also bundled together with the source 
code for this program. That content is listed below along with the licenses 
under which they have been released.

> **OBO Ontology python module**  
http://pythonhosted.org/Orange-Bioinformatics/  
Copyright (c) Bioinformatics Laboratory, FRI UL  
Released under the GNU General Public License license

> **Sortable tables**  
http://www.kryogenix.org/code/browser/sorttable/  
Copyright (c) Stuart Langridge   
Released under the X11 (MIT) license  
http://www.kryogenix.org/code/browser/licence.html  

> **jQuery**  
  Copyright (c) jQuery Foundation and other contributors  
  Released under the MIT license  
  http://jquery.org/license  

> **SVGPan v3.2.9**  
  https://github.com/ariutta/svg-pan-zoom  
  Copyright (c) Andrea Leofreddi  
  The code from the SVGPan library is licensed under the following BSD license  
  https://raw.githubusercontent.com/ariutta/svg-pan-zoom/master/LICENSE  

> **Tabbed Content v2013.7.6**  
  http://www.menucool.com/tabbed-content  
  Free to use