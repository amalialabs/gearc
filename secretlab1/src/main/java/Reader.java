import net.minidev.json.parser.ParseException;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Stream;
import java.util.zip.GZIPInputStream;

/**
 * @author hadziahmetovic on 23.01.22
 */
public class Reader {
    /*
     add GO reader
     add ensembl/gaf reader -> autoconvert/scrape info to generate ensembl verison
     diff-output reader -> flexible with col selector -> make gene, fc, label format??
     */

    HashMap<String, Gene> geneMap;
    //fixme need to be set before calling
    int numGenesTotal;
    int deGenes;


    /**
     * Has to have header
     * id      fc   fdr      signif
     * DNAJC25-GNG10   -1.3420  0.1 false
     * IGKV2-28        -2.3961 0.12  false
     * @param expressionFile
     * @return
     */
    public Set<Gene> readExpressionFile(File expressionFile, boolean toConvert, HashMap<String, String> nameMapping) {
        Set<Gene> genes = new HashSet<>();

        SplittableRandom r = new SplittableRandom();
        try (Stream<String> stream = Files.lines(Paths.get(expressionFile.getAbsolutePath()))) {
            stream.skip(1).forEach(_line -> {
                if (_line.charAt(0) != '#' && _line.charAt(0) != 'i') {
//                    System.out.println(_line);
                    String[] elems = _line.split("\t");
                    String gene_id;
                    if (toConvert) {
                        gene_id = elems[0];
                    } else {
                        gene_id = nameMapping.get(elems[0]);
                    }
                    double fc = Double.parseDouble(elems[1]);
                    double fdr = r.nextDouble(1); //Double.parseDouble(elems[2]);

                    Gene g = new Gene(gene_id, fc, fdr);
                    genes.add(g);
                }
            });
        } catch (IOException e) {
            throw new RuntimeException("Error reading expression file: ", e);
        }
        return genes;
    }

    public Set<Gene>readExpressionFile(File expressionFile) {
        return readExpressionFile(expressionFile, true, null);
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
    public GO readOboFile(File oboFile, String root) {
        GO gos = new GO();
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
                            gos.getGoNodes().put(id, node);
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
            System.out.println("IOexception at obo file");
            throw new RuntimeException(e);
        }

        return gos;
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
    public void readMappringGAF(File mappingFile, GO gos){
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
//    private void readMappringEnsebl(){
//
//        String line, name;
//        BufferedReader br = null;
//        int one, two;
//        boolean b;
//
//        try {br = new BufferedReader(new FileReader(mapping));} catch (FileNotFoundException e) {}
//        try {
//            line = br.readLine();
//            while ((line = br.readLine()) != null) {
//                one = line.indexOf("\t");
//                two = line.indexOf("\t", one + 1);
//                if (two - one < 2) continue;
//                name = line.substring(one + 1, two);
//                if (!geneSignif.keySet().contains(name)) {
//                    continue;
//                }
//                //currentSet = new HashSet<>();
//                currentGene = name;
//
//                b = geneSignif.get(name);
//                getGOs(line.substring(two + 1), b);
//                //geneGO.put(name, currentSet);
//            }
//            br.close();
//        }
//        catch (IOException e) {
//            System.out.println("IOexception ensembl reader");
//        }
//    }
//
//    private void getGOs (String line, boolean b) {
//        try {
//            while (true) {
//                if (!allGO.contains(line.substring(0, 10))) {
//                    line = line.substring(11);
//                    continue;
//                }
//                //currentSet.add(line.substring(0, 10));
//                if (!b) {   //TODO dont ask why, it works
//                    if (goPositive.get(line.substring(0, 10)) == null) {
//                        Set<String> set123 = new HashSet<>();
//                        set123.add(currentGene);
//                        goPositive.put(line.substring(0, 10), set123);
//                    } else {
//                        goPositive.get(line.substring(0, 10)).add(currentGene);
//                    }
//                } else {
//                    if (goNegative.get(line.substring(0, 10)) == null) {
//                        Set<String> set123 = new HashSet<>();
//                        set123.add(currentGene);
//                        goNegative.put(line.substring(0, 10), set123);
//                    } else {
//                        goNegative.get(line.substring(0, 10)).add(currentGene);
//                    }
//                }
//                line = line.substring(11);
//            }
//        } catch (Exception e) {}
//    }
}
