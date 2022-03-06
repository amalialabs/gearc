import java.io.*;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import java.util.zip.GZIPInputStream;

public class Reader {
    /*
     add GO reader
     add ensembl/gaf reader -> autoconvert/scrape info to generate ensembl verison
     diff-output reader -> flexible with col selector -> make gene, fc, label format??
     */
    HashMap<String, Gene> geneMap = new HashMap<>();
    HashMap<String, Set<String>> geneToGO = new HashMap<>();


    //fixme need to be set before calling
    int numGenesTotal;
    int deGenes;

    public Reader() {
    }

    //fixme specify and rely on ensembl mappingFile; updates may contain more options
    // - also add mutable root -> perform for each root at the same time? maybe not, but add full call with all plots and per root
    public Reader(File expressionFile, File mappingFile, File oboFile, String root) {
        System.out.println("\nObo: starting");
        long time = System.currentTimeMillis();
        GO.goNodes = readOboFile(oboFile, root);
        System.out.println("\tObo time: "+(System.currentTimeMillis()-time)+" ms");

        System.out.println("\nExpressionFile: starting");
        time = System.currentTimeMillis();
        readExpressionFile(expressionFile);
        System.out.println("\tExpressionFile time: "+(System.currentTimeMillis()-time)+" ms");

        System.out.println("\nMappingFile: starting");
        time = System.currentTimeMillis();
        readMappringEnsebl(mappingFile);
        System.out.println("\tMappingFile time: "+(System.currentTimeMillis()-time)+" ms");
        postprocess();
    }

    /**
     * Has to have header
     * id      fc   fdr
     * DNAJC25-GNG10   -1.3420  0.1
     * IGKV2-28        -2.3961 0.12
     * @param expressionFile
     * @return
     */
    private void readExpressionFile(File expressionFile, boolean toConvert, HashMap<String, String> nameMapping) {
        System.out.println("Starting");
        Set<Gene> genes = new HashSet<>();
        SplittableRandom r = new SplittableRandom();    //fixme
        try (Stream<String> stream = Files.lines(Paths.get(expressionFile.getAbsolutePath()))) {
            stream.skip(1).forEach(_line -> {
                if (_line.charAt(0) != '#' && _line.charAt(0) != 'i') { //fixme -> should not have any lines with #GOxxx
                    String[] elems = _line.split("\t");
                    String gene_id;
                    if (toConvert) {
                        gene_id = elems[0];
                    } else {
                        gene_id = nameMapping.get(elems[0]);
                    }
                    double fc = Double.parseDouble(elems[1]);
                    double fdr = r.nextDouble(1); //Double.parseDouble(elems[2]); //fixme

                    Gene g = new Gene(gene_id, fc, fdr);
                    genes.add(g);
                }
            });
        } catch (IOException e) {
            throw new RuntimeException("Error reading expression file: ", e);
        }
        Functions.filter_unclear(genes).forEach(_g -> geneMap.put(_g.gene_id, _g));

        System.out.println("---------");
        System.out.println("Filter check:");
        System.out.println("Input:\t" + genes.size());
        System.out.println("Filter:\t" + geneMap.keySet().size());
        System.out.println("---------");
    }

    private void readExpressionFile(File expressionFile) {
        readExpressionFile(expressionFile, true, null);
    }

    /**
     * [Term]
     * id: GO:0000001
     * name: mitochondrion inheritance
     * namespace: biological_process
     * def: "The distribution of mitochondria, including the mitochondrial genome, into daughter cells after mitosis or meiosis, mediated by interactions between mitochondria and the cytoskeleton." [GOC:mcc, PMID:10873824, PMID:11389764]
     * synonym: "mitochondrial inheritance" EXACT []
     * is_a: GO:0048308 ! organelle inheritance
     * is_a: GO:0048311 ! mitochondrion distribution
     *
     * @param oboFile
     * @return
     */
    private HashMap<String, Node> readOboFile(File oboFile, String root) {
        HashMap<String, Node> goNodes = new HashMap<>();
        String line, id = null, name, namespace = null;
        Set<String> set = null;
        boolean good = true;
        Node node = null;

        BufferedReader br = null;
        try {br = new BufferedReader(new FileReader(oboFile));} catch (FileNotFoundException e) {}
        try {
            while ((line = br.readLine()) != null) {
                if (line.length() > 1) {
                    if (line.charAt(0) == '[' && line.charAt(2) == 'y') break;

                    else if (line.charAt(0) == 'i' && line.charAt(1) == 'd') {
                        if (!good) { }  //TODO nice if man
                        else if (set != null) {
                            if (node != null) {
                                node.setParents(set);
                            }
                            goNodes.put(id, node);
                        }
                        id = line.substring(4);
                        set = new HashSet<>();
                        good = true;
                        node = new Node();
                        node.setNode_id(id);
                    } else if (line.charAt(0) == 'n' && line.charAt(4) == ':') {
                        name = line.substring(6);
                        if (name.contains("obsolete")) good = false;
                        else node.setNode_name(name);
                    } else if (line.charAt(0) == 'n' && line.charAt(4) == 's') {
                        namespace = line.substring(11);
                        if (!namespace.equals(root)) {
                            good = false;
                        }
                    } else if (line.charAt(1) == 's' && line.charAt(2) == '_' && line.charAt(3) == 'a') {
                        if (namespace.equals(root)) {
                            set.add(line.substring(6, line.indexOf("!", 7) - 1));
                        }
                    }
                }
            }
            br.close();
            if (set != null) {
                if (node != null) {
                    node.setParents(set);
                }
            }
        } catch (IOException e) {
            throw new RuntimeException("IOexception at obo file " + oboFile.getAbsolutePath(), e);
        }
        return goNodes;
    }

    /**
     * 3. column = geneName
     * 4. column = goName
     * UniProtKB       A0A024R161      DNAJC25-GNG10           GO:0004871      GO_REF:0000038  IEA     UniProtKB-KW:KW-0807    F       Guanine nucleotide-binding protein subunit gamma
     * UniProtKB       A0A024R161      DNAJC25-GNG10           GO:0005834      GO_REF:0000002  IEA     InterPro:IPR001770|InterPro:IPR015898   C       Guanine nucleotide-binding protei
     * UniProtKB       A0A024R161      DNAJC25-GNG10           GO:0007186      GO_REF:0000002  IEA     InterPro:IPR001770|InterPro:IPR015898   P       Guanine nucleotide-binding protei
     * UniProtKB       A0A075B6P5      IGKV2-28                GO:0002250      GO_REF:0000037  IEA     UniProtKB-KW:KW-1064    P       Immunoglobulin kappa variable 2-28      KV228_H
     * @param mappingFile
     * @param gos
     * FIXME chnage ensembl rest to automatically map names
     */
    private void readMappringGAF(File mappingFile, GO gos){
        int c = 0;
        String line;
        BufferedReader br = null;
        int one, two, three, four, five;
        String gene, go;

        int lines = 2000;

        try {br = new BufferedReader(new InputStreamReader( new GZIPInputStream(new FileInputStream(mappingFile))));} catch (Exception e) {}
        try {
            while ((line = br.readLine()) != null){
                if (line.charAt(0) == '!') continue;
                one = line.indexOf("\t");
                two = line.indexOf("\t", one+1);
                three = line.indexOf("\t", two+1);
                four = line.indexOf("\t", three+1);
                five = line.indexOf("\t", four+1);
                //System.out.println(line);


                if(c++ > lines) break;



                if (four-three > 1) {
                    //System.out.println(line);
                    continue;
                }
                gene = line.substring(two+1, three);

                //todo add this check for later -> redundant info removal
//                if (geneSignif.get(gene) == null) {
//                    continue;
//                }
                go = line.substring(four+1, five);
                Node node = gos.getGoNodes().get(go);
                if (node == null) {
                    continue;
                }
                //TODO catch unmappable ids
                String geneId = "";
                try {
                    geneId = EnsemblRestClient.getGeneID("human", gene);
                } catch (Exception e) {
                    continue;
                }
//                System.out.println(gene + "\t" + geneId);
//                System.out.println("Gene\t"+gene);
//                System.out.println("Go\t"+go);
//                System.out.println(++c);
                if(geneMap.get(geneId) != null) {
                    if (node.getGenes() == null) {
                        Set<Gene> set = new HashSet<>();
                        set.add(geneMap.get(geneId));
                        node.setGenes(set);
                    } else {
                        node.getGenes().add(geneMap.get(geneId));
                    }
                }

                //todo assignment for up/downregulation targets -> will se e later how to handle
//                if (!geneSignif.get(gene)) { //TODO my bad, should be true but its just for consistency now
//                    if (goPositive.get(go) != null) {
//                        goPositive.get(go).add(gene);
//                    } else {
//                        Set<String> set = new HashSet<>();
//                        set.add(gene);
//                        goPositive.put(go, set);
//                    }
//                } else {
//                    if (goNegative.get(go) != null) {
//                        goNegative.get(go).add(gene);
//                    } else {
//                        Set<String> set = new HashSet<>();
//                        set.add(gene);
//                        goNegative.put(go, set);
//                    }
//                }

            }
            br.close();
        }
        catch (IOException e) {
            System.out.println("IOexception gaf reader");
        }
    }


    //todo get creator functions from gergely
    private void readMappringEnsebl(File ensemblMapping){
        try (Stream<String> stream = Files.lines(Paths.get(ensemblMapping.getAbsolutePath()))) {
            final int[] counterLine = {0};
            final int[] counterEntry = {0};
            stream.skip(1).forEach(_line -> { //has header
                String[] elems = _line.split("\t");
                String gene_id = elems[0];
                String gene_name = elems[1];
                counterLine[0]++;
                if(geneMap.containsKey(gene_name)){
                    geneToGO.putIfAbsent(gene_name, new HashSet<>());
                    geneToGO.get(gene_name).addAll(Stream.of(elems[2].split("\\|"))
                            .filter(_go -> GO.goNodes.containsKey(_go)).collect(Collectors.toSet()));
                    counterEntry[0]++;
                }
                //fixme take care of either gene _id or _name
            });
            System.out.println("---------");
            System.out.println("Mapping check:");
            System.out.println("Etries:\t" + counterLine[0]);
            System.out.println("Valid Entries:\t" + counterEntry[0]);
            System.out.println("Unique Entries:\t" + geneToGO.keySet().size());
            System.out.println("---------");

        } catch (IOException e) {
            throw new RuntimeException("Error during parsing of " + ensemblMapping.getAbsolutePath(), e);
        }
    }

    /**
     * Propagates gene counts to all parent GOs.
     *
     * Collects all relevant entries, according to min/max criteria, in relevantGo.
     * Establishes parent-child connection.
     */
    private void postprocess() {
        geneToGO.forEach((key, value) -> {
            for (String go_id : value) {
                GO.getGoNodes().get(go_id).getGenes().add(geneMap.get(key));
            }
        });

        Set<Node> gos = new HashSet<>(GO.goNodes.values());
        Set<Node> notPropagated;

        System.out.println("Pre propagating posneg size\t"+gos.size());
        do {
            notPropagated = new HashSet<>();
            Set<Node> finalNotPropagated = notPropagated;
            gos.forEach(_go -> propagate(_go, finalNotPropagated));
            notPropagated = finalNotPropagated;
        } while ((gos = notPropagated).size() > 0);


//        posneg = new HashSet<>();
//
//        goParent.keySet().forEach(x -> {
//            if (goParent.get(x) != null) {
//                if(goPositive.get(x) != null || goNegative.get(x) != null) {
//                    goPositive.putIfAbsent(x, new HashSet<>());
//                    goNegative.putIfAbsent(x, new HashSet<>());
//                    if (goPositive.get(x).size() + goNegative.get(x).size() >= minSize && goPositive.get(x).size() + goNegative.get(x).size() <= maxSize) {
//                        relevantGo.add(x);
//
//                        //relevantGo.put(x, new HashSet<>(goParent.get(x)));
//                    }
//                }
////                for (String parent : goParent.get(x)) { //TODO eval if this is useful
////                    if (goChild.get(parent) == null) {
////                        Set<String> se = new HashSet<>();
////                        se.add(x);
////                        goChild.put(parent, se);
////                    } else {
////                        goChild.get(parent).add(x);
////                    }
////                }
//            }
//            goPositive.putIfAbsent(x, new HashSet<>());
//            goNegative.putIfAbsent(x, new HashSet<>());
//            posneg.addAll(goPositive.get(x));
//            posneg.addAll(goNegative.get(x));
//        });
//
//        for (String gene : posneg) {
//            if (geneSignif.get(gene)) geneSignifCount++;
//        }
//        System.out.println("Number of relevant genes\t" + posneg.size());
        //TODO need step for relevant genes
    }

    private void propagate(Node n, Set<Node> notPropagated) {
        for (Node parent : GO.getGoNodes().get(n.node_id).getParentNodes()) {
            parent.getGenes().addAll(n.getGenes());
        }
        notPropagated.addAll(GO.getGoNodes().get(n.node_id).getParentNodes());
    }

    public static void main(String[] args) {
        File expression = new File("/home/birinci/GOEnrichment/simul_exp_go_bp_ensembl.tsv");
        File mappingEnsembl = new File("/home/birinci/GOEnrichment/goa_human_ensembl.tsv");
        File oboFile = new File("/home/birinci/GOEnrichment/go.obo");
        String root = "biological_process";

        Reader r = new Reader(expression, mappingEnsembl, oboFile, root);

        GO.getGoNodes().values().forEach(_node -> {
            if(_node.getGenes() != null && _node.getGenes().size() > 3000) {
                System.out.println(_node.node_id);
                System.out.println(_node.getNode_name());
                System.out.println(_node.getGenes().size());
            }
        });
    }
}
