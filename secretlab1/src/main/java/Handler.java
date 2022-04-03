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
        accepts("root").withRequiredArg().ofType(String.class).defaultsTo("biological_process").describedAs("GO DAG root");
        accepts("mapping").withRequiredArg().ofType(String.class).defaultsTo("/data/goa_human_ensembl.tsv").describedAs("ENSEMBL gene mapping");
        accepts("obo").withRequiredArg().ofType(String.class).defaultsTo("/data/go.obo").describedAs("GO DAG");
        accepts("out").withRequiredArg().ofType(String.class).defaultsTo("/out/").describedAs("output directory"); //LATER defaults to installation path
        }};

        OptionSet params = optionParser.parse(args);
        assertTrue(params.has("genelist"));

        System.out.println("Params:"); //FIXME doesnt print options but only some java-id
        params.specs().forEach(spec -> {
            System.out.println(spec + "\t" + params.valueOf(spec));
        });

        try {
            if (params.has("h") || params.has("?")) {
                optionParser.printHelpOn(System.out);
                System.exit(0);
            }
        } catch (IOException e) {
            throw new RuntimeException("Could not display help page.", e);
        }

        GO gos = null;
        File obo = new File((String) params.valueOf("obo"));
        File mapping = new File((String) params.valueOf("mapping"));
        File expression = new File((String) params.valueOf("genelist"));
        String root = (String) params.valueOf("root");
        Reader r = new Reader(expression, mapping, obo, root);

        String outdir = (String) params.valueOf("out");

        Functions.score_genes(new HashSet<>(r.geneMap.values()));

        Enum expected_change = (Enum) params.valueOf("expectedChange");
        Functions.assign_genes_to_sets(new HashSet<>(r.geneMap.values()), expected_change);

        Set<Gene> genesTotal = new HashSet<>();
        Set<Gene> genesSignif = new HashSet<>();
        for (Node node : GO.getGoNodes().values()) {
            if (node.getGenes() != null) {
                for (Gene gene : node.getGenes()) {
                    if (gene.is_significant) {
                        genesSignif.add(gene);
                    }
                }
                genesTotal.addAll(node.getGenes());
            }
        }
        int numGenesTotal = genesTotal.size();
        int deGenes = genesSignif.size();
        System.out.println(numGenesTotal + "\t" + deGenes);
        Enrichment en = new Enrichment(numGenesTotal, deGenes);


        Result result = new Result();  //alternative
        for (int i = 0; i < 1000; i++) {
            Set<Gene> sampled = Functions.sample_genes(new HashSet<>(r.geneMap.values()), 0.2);
            result.gather_runs(en.enrich(sampled, gos));
        }

        double percent = Functions.extend_flex_set(new HashSet<>(r.geneMap.values())); //alternative, will do only if percent > 0,2
        System.out.println(percent);
        if (percent > 0.2) {
            for (int i = 0; i < 1000; i++) {
                Set<Gene> sampled = Functions.sample_genes(new HashSet<>(r.geneMap.values()), percent);
                result.gather_runs(en.enrich(sampled, gos));
            }
        }

        Set<Node> robust_gos = result.getXquantileGOnodes(0.95);

        Plots plots = new Plots(outdir, r.geneMap.values(), robust_gos, (double) params.valueOf("FDR"), (double) params.valueOf("FC"), result);
        plots.unclear_genes_BARPLOT(r.allGenes.values());
        plots.sig_genes_VOLCANO();
        plots.gene_categories_BARPLOT();
        plots.weighted_genes_CUMULATIVE();
        plots.genes_set_PIECHART();

        if (params.has("out")) {
            String filepath = (String) params.valueOf("out");
            result.writeRobustGOs(robust_gos, filepath);
        } else {
            result.printRobustGOs(robust_gos);
        }
    }
}
