import joptsimple.OptionParser;
import joptsimple.OptionSet;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

import static java.util.Arrays.asList;
import static org.junit.Assert.*;

public class Handler {

    public enum expected_change {LOW, AVERAGE, HIGH};

    public static void main(String[] args) {
        OptionParser optionParser = new OptionParser() {{
        accepts("genelist").withRequiredArg().required().ofType(String.class).describedAs("diffexp output");
        accepts("expectedChange").withRequiredArg().ofType(expected_change.class).defaultsTo(expected_change.AVERAGE).describedAs("expected change in diffexp data");
        accepts("FC").withRequiredArg().ofType(double.class).defaultsTo(1.0).describedAs("log2foldchange cutoff");
        accepts("FDR").withRequiredArg().ofType(double.class).defaultsTo(0.01).describedAs("adjusted p-value cutoff");
        acceptsAll( asList( "h", "?" ), "show help" ).forHelp();
        accepts("html");
        accepts("quantile").withRequiredArg().ofType(double.class).defaultsTo(0.95).describedAs("desired robustness of GO enrichment results");
        accepts("n").withRequiredArg().ofType(Integer.class).defaultsTo(5).describedAs("number of enrichment iterations");
        accepts("root").withRequiredArg().ofType(String.class).defaultsTo("biological_process").describedAs("GO DAG root");
        accepts("mapping").withRequiredArg().ofType(String.class).defaultsTo("/data/goa_human_ensembl.tsv").describedAs("ENSEMBL gene mapping");
        accepts("obo").withRequiredArg().ofType(String.class).defaultsTo("/data/go.obo").describedAs("GO DAG");
        accepts("out").withRequiredArg().ofType(String.class).defaultsTo("/out/").describedAs("output directory"); //LATER defaults to installation path
        }};

        OptionSet params = optionParser.parse(args);
        assertTrue(params.has("genelist"));

        System.out.println("Params:");
        params.specs().forEach(spec -> System.out.println(spec + "\t" + params.valueOf(spec)));

        try {
            if (params.has("h") || params.has("?")) {
                optionParser.printHelpOn(System.out);
                System.exit(0);
            }
        } catch (IOException e) {
            throw new RuntimeException("Could not display help page.", e);
        }

        double FDR_cutoff = (double) params.valueOf("FDR");
        double FC_cutoff = (double) params.valueOf("FC");
        System.out.println("FDR cutoff: " + FDR_cutoff + "\nFC cutoff: " + FC_cutoff);

        GO gos = null;
        File obo = new File((String) params.valueOf("obo"));
        File mapping = new File((String) params.valueOf("mapping"));
        File expression = new File((String) params.valueOf("genelist"));
        String root = (String) params.valueOf("root");
        Reader r = new Reader(expression, mapping, obo, root, FDR_cutoff, FC_cutoff);

        String outdir = (String) params.valueOf("out");


        Functions.score_genes(new HashSet<>(r.geneMap.values()), FDR_cutoff, FC_cutoff);

        Enum expected_change = (Enum) params.valueOf("expectedChange");
        Functions.assign_genes_to_sets(new HashSet<>(r.geneMap.values()), expected_change);

        Set<Gene> genesTotal = new HashSet<>();
        Set<Gene> genesSignif = new HashSet<>();
        Set<Gene> genesTotalFiltered = new HashSet<>();
        for (Node node : GO.getGoNodes().values()) {
            if (node.getGenes() != null) {
                for (Gene gene : node.getGenes()) {
                    if (gene.is_significant) {
                        genesSignif.add(gene);
                    }
                    if (!gene.unclear) {
                        genesTotalFiltered.add(gene);
                    }
                }
                genesTotal.addAll(node.getGenes());
            }
        }
        int numGenesTotal = genesTotal.size();
        int deGenes = genesSignif.size();
        int numGenesTotalFiltered = genesTotalFiltered.size();
        System.out.println("#used genes: " + numGenesTotalFiltered + "\tof which diff genes: " + deGenes);
        Enrichment en = new Enrichment(numGenesTotalFiltered, deGenes);

        Set<Gene> allGenes = new HashSet<>();
        allGenes.addAll(r.allGenes.values());
        Enrichment en_standard = new Enrichment(numGenesTotal, deGenes);
        HashMap<Node, Double> standard_node2fdr = en_standard.enrich(allGenes, gos);  //standard enrichment way


        Result result = new Result(FDR_cutoff);  //alternative
        int numIter = (int) params.valueOf("n");
        for (int i = 0; i < numIter; i++) {
            Set<Gene> sampled = Functions.sample_genes(new HashSet<>(r.geneMap.values()), 0.2);
            result.gather_runs(en.enrich(sampled, gos), false);
        }

//        double percent = Functions.extend_flex_set(new HashSet<>(r.geneMap.values())); //alternative, will do only if percent > 0,2
//        System.out.println("extending flex set up to " + Functions.df.format(percent*100)+  "%");       //TODO
//        if (percent > 0.2) {
//            for (int i = 0; i < numIter; i++) {
//                Set<Gene> sampled = Functions.sample_genes(new HashSet<>(r.geneMap.values()), percent);
//                result.gather_runs(en.enrich(sampled, gos), true);
//            }
//        }

        double quantile = (double) params.valueOf("quantile");
        Set<Node> robust_gos = result.getXquantileGOnodes(quantile);

        Plots plots = new Plots(outdir, r.allGenes.values(), robust_gos, FDR_cutoff, FC_cutoff, result, standard_node2fdr);
        plots.runPlots(r);

        if (params.has("out")) {
            result.writeRobustGOs(robust_gos, outdir, quantile);
            result.writeMeanGOs(robust_gos, outdir);
            result.writeStandardGOs(standard_node2fdr, outdir);
            result.writeRankDifferences(standard_node2fdr, robust_gos, outdir);
        } else {
            result.printRobustGOs(robust_gos);
            result.printStandardGOs(standard_node2fdr);
        }
    }
}
