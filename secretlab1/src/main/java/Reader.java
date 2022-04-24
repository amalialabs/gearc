import java.io.*;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class Reader {
    /*
     add GO reader
     add ensembl/gaf reader -> autoconvert/scrape info to generate ensembl verison
     diff-output reader -> flexible with col selector -> make gene, fc, label format??
     */
    HashMap<String, Gene> allGenes = new HashMap<>();
    HashMap<String, Gene> geneMap = new HashMap<>();
    HashMap<String, Set<String>> geneToGO = new HashMap<>();
    HashMap<String, String> geneID2Name = new HashMap<>();
    HashMap<String, String> geneName2ID = new HashMap<>();

    //fixme ~~specify and rely on ensembl mappingFile~~; updates may contain more options
    // - also add mutable root -> perform for each root at the same time? maybe not, but add full call with all plots and per root
    public Reader(File expressionFile, File mappingFile, File oboFile, String root) {
        System.out.println("\nObo: starting");
        long time = System.currentTimeMillis();
        GO.goNodes = readOboFile(oboFile, root);
        System.out.println("\tObo time: " + (System.currentTimeMillis() - time) + " ms");


        System.out.println("\nMappingFile: starting");
        time = System.currentTimeMillis();
        readMappringEnsebl(mappingFile);
        System.out.println("\tMappingFile time: " + (System.currentTimeMillis() - time) + " ms");


        System.out.println("\nExpressionFile: starting");
        time = System.currentTimeMillis();
        readExpressionFile(expressionFile);
        System.out.println("\tExpressionFile time: " + (System.currentTimeMillis() - time) + " ms");


        System.out.println("\nPostprocessing: starting");
        time = System.currentTimeMillis();
        postprocess();
        System.out.println("\tPostprocessing time: " + (System.currentTimeMillis() - time) + " ms");
    }

    /**
     * Has to have header
     * id      fc   fdr
     * DNAJC25-GNG10   -1.3420  0.1
     * IGKV2-28        -2.3961 0.12
     *
     * @param expressionFile
     * @return
     */
    private void readExpressionFile(File expressionFile, boolean isGeneID) {
        Set<Gene> genes = new HashSet<>();
        try (Stream<String> stream = Files.lines(Paths.get(expressionFile.getAbsolutePath()))) {
            stream.skip(1).forEach(_line -> {
                if (_line.charAt(0) != '#' && _line.charAt(0) != 'i') {
                    String[] elems = _line.split("\t");
                    String gene_id;
                    if (isGeneID) {
                        gene_id = elems[0];
                    } else {
                        gene_id = geneName2ID.get(elems[0]);
                    }
                    double fc = Double.parseDouble(elems[1]);
                    double fdr = Double.parseDouble(elems[2]);

                    Gene g = new Gene(gene_id, geneID2Name.get(gene_id), fc, fdr);
                    genes.add(g);
                }
            });
        } catch (IOException e) {
            throw new RuntimeException("Error reading expression file: ", e);
        }
        genes.forEach(_g -> allGenes.put(_g.gene_id, _g));
        genes.forEach(_g -> geneMap.put(_g.gene_id, _g));
//        Functions.filter_unclear(genes).forEach(_g -> geneMap.put(_g.gene_id, _g));

        System.out.println("---------");
        System.out.println("Filter check:");
        System.out.println("Input:\t" + genes.size());
        System.out.println("Filter:\t" + geneMap.keySet().size());
        System.out.println("---------");
    }

    private void readExpressionFile(File expressionFile) {
        readExpressionFile(expressionFile, true);
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
        try {
            br = new BufferedReader(new FileReader(oboFile));
        } catch (FileNotFoundException e) {
            throw new RuntimeException("Missing obo input file: ", e);
        }
        try {
            while ((line = br.readLine()) != null) {
                if (line.length() > 1) {
                    if (line.charAt(0) == '[' && line.charAt(2) == 'y') break;

                    else if (line.charAt(0) == 'i' && line.charAt(1) == 'd') {
                        if (good && set != null) {
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


    private void readMappringEnsebl(File ensemblMapping) {
        try (Stream<String> stream = Files.lines(Paths.get(ensemblMapping.getAbsolutePath()))) {
            final int[] counterLine = {0};
            final int[] counterEntry = {0};
            stream.skip(1).forEach(_line -> { //has header
                String[] elems = _line.split("\t");
                String gene_id = elems[0];
                String gene_name = elems[1];
                geneID2Name.put(gene_id, gene_name);
                geneName2ID.put(gene_name, gene_id);
                counterLine[0]++;
                geneToGO.putIfAbsent(gene_id, new HashSet<>());
                geneToGO.get(gene_id).addAll(Stream.of(elems[2].split("\\|"))
                        .filter(_go -> GO.goNodes.containsKey(_go)).collect(Collectors.toSet()));
                counterEntry[0]++;

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
     * <p>
     * Collects all relevant entries, according to min/max criteria, in relevantGo.
     * Establishes parent-child connection.
     */
    private void postprocess() {
        geneToGO.keySet().removeIf(gene_id -> !geneMap.containsKey(gene_id));

        geneToGO.forEach((key, value) -> {
            for (String go_id : value) {
                GO.getGoNodes().get(go_id).getGenes().add(geneMap.get(key));
            }
        });

        Set<Node> gos = new HashSet<>(GO.goNodes.values());
        Set<Node> notPropagated;

        do {
            notPropagated = new HashSet<>();
            Set<Node> finalNotPropagated = notPropagated;
            gos.forEach(_go -> propagate(_go, finalNotPropagated));
            notPropagated = finalNotPropagated;
        } while ((gos = notPropagated).size() > 0);
    }

    private void propagate(Node n, Set<Node> notPropagated) {
        for (Node parent : GO.getGoNodes().get(n.node_id).getParentNodes()) {
            parent.getGenes().addAll(n.getGenes());
        }
        notPropagated.addAll(GO.getGoNodes().get(n.node_id).getParentNodes());
    }

    public static void main(String[] args) {
        File expression = new File("/home/birinci/Projects/secretlab/data/empires.out");
        File mappingEnsembl = new File("/home/birinci/Projects/secretlab/data/goa_human_ensembl.tsv");
        File oboFile = new File("/home/birinci/Projects/secretlab/data/go.obo");
        String root = "biological_process";

        Reader r = new Reader(expression, mappingEnsembl, oboFile, root);

        GO.getGoNodes().values().forEach(_node -> {
            if (_node.getGenes() != null && _node.getGenes().size() > 3000) {
                System.out.println(_node.node_id);
                System.out.println(_node.getNode_name());
                System.out.println(_node.getGenes().size());
            }
        });
    }
}
