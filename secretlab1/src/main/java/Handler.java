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
        System.out.println(FDR_cutoff + "  " + FC_cutoff);

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
        System.out.println(numGenesTotalFiltered + "\t" + deGenes); //TODO
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

        double percent = Functions.extend_flex_set(new HashSet<>(r.geneMap.values())); //alternative, will do only if percent > 0,2
        System.out.println(percent);       //TODO
        if (percent > 0.2) {
            for (int i = 0; i < numIter; i++) {
                Set<Gene> sampled = Functions.sample_genes(new HashSet<>(r.geneMap.values()), percent);
                result.gather_runs(en.enrich(sampled, gos), true);
            }
        }

        Set<Node> robust_gos = result.getXquantileGOnodes(0.95);

        Plots plots = new Plots(outdir, r.allGenes.values(), robust_gos, (double) FDR_cutoff, FC_cutoff, result, standard_node2fdr);
        System.out.println("plotting unclear_genes_BARPLOT");
        plots.unclear_genes_BARPLOT(r.allGenes.values());
        System.out.println("plotting sig_genes_VOLCANO");
        plots.sig_genes_VOLCANO();
        System.out.println("plotting gene_categories_BARPLOT");
        plots.gene_categories_BARPLOT();
        System.out.println("plotting weigthed_genes_CUMULATIVE");
        plots.weighted_genes_CUMULATIVE();
        System.out.println("plotting genes_sets_PIECHART");
        plots.genes_set_PIECHART();
        System.out.println("plotting go_fdrs_mean_vs_quantiole_SCATTER");
        plots.go_fdrs_mean_vs_quantile_SCATTER();
        System.out.println("plotting selected_gos_fdr_distrib_BOXPLOT");
        plots.selected_gos_fdr_distrib_BOXPLOT();
        System.out.println("plotting selected_gos_rob_vs_extend_BOXPLOT");
        plots.selected_gos_rob_vs_extend_BOXPLOT();
        System.out.println("plotting gos_quntile_vs_mean_fdr_JITTER");
        plots.gos_quantile_vs_mean_fdr_JITTER();
        System.out.println("plotting gos_standard_vs_robust_vs_extend_BARPLOT");
        plots.gos_standard_vs_robust_vs_extended_BARPLOT();
        System.out.println("plotting gos_standard_vs_robust_vs_extend_VENN");
        plots.gos_standard_vs_robust_vs_extended_VENN();
        System.out.println("plotting fdr cutoff finding");
        plots.fdr_cutoff_finding_CUMULATIVE();
        System.out.println("plotting fc cutoff finding");
        plots.fc_cutoff_finding_CUMULATIVE();
        System.out.println("plotting flexset extension");
        plots.flexset_extension_CURVE();
        System.out.println("plotting expected change distrib");
        plots.expected_change_BARPLOT();
        System.out.println("plotting num sig gos comparison");
        plots.num_sig_gos_BARPLOT();

        if (params.has("out")) {
            result.writeRobustGOs(robust_gos, outdir);
            result.writeStandardGOs(standard_node2fdr, outdir);
        } else {
            result.printRobustGOs(robust_gos);
            result.printStandardGOs(standard_node2fdr);
        }
    }
}
