import joptsimple.OptionParser;
import joptsimple.OptionSet;
import java.io.File;
import java.io.IOException;
import java.util.*;

import static java.util.Arrays.asList;
import static org.junit.Assert.*;

public class Handler {

    public enum expected_change {LOW, AVERAGE, HIGH};

    public static void main(String[] args) {
        OptionParser optionParser = new OptionParser() {{
        accepts("genelist").withRequiredArg().required().ofType(File.class).describedAs("diffexp output");
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

        System.out.println(params.toString());

        try {
            optionParser.printHelpOn(System.out); //FIXME print only if --help not always
        } catch (IOException e) {
            throw new RuntimeException("Could not display help page.", e);
        }

        GO gos = null;
        File obo = new File((String) params.valueOf("obo"));
        File mapping = new File((String) params.valueOf("mapping"));
        //File expression = (File) params.valueOf("genelist");
        File expression = new File("/data/simul_exp_go_bp_ensembl.tsv"); //for testing only
        String root = (String) params.valueOf("root");
        Reader r = new Reader(expression, mapping, obo, root);

        String outdir = (String) params.valueOf("out");
        Plots plots = new Plots(outdir, r.geneMap.values(), (double) params.valueOf("FDR"), (double) params.valueOf("FC"));

        //plots.unclear_genes_BARPLOT(r.geneMap.values()); //FIXME must be called before genes are filtered!!!
        plots.sig_genes_VOLCANO();
        plots.gene_categories_BARPLOT();

        Functions.score_genes(new HashSet<>(r.geneMap.values()));

        Enum expected_change = (Enum) params.valueOf("expectedChange");
        Functions.assign_genes_to_sets(new HashSet<>(r.geneMap.values()), expected_change);


        int numGenesTotal = GO.getGoNodes().values().stream()
                .filter(_node -> _node.getGenes() != null).mapToInt(_node -> _node.getGenes().size()).sum();
        //todo need to remove "null"genes from Nodes
        Set<Gene> ge = new HashSet<>();
        GO.getGoNodes().values().stream().forEach(_node -> {
            if (_node.getGenes() != null) ge.addAll(_node.getGenes());
        });
        int deGenes = (int) ge.stream().filter(_go -> _go.is_significant).count();
        System.out.println(numGenesTotal + "\t" + deGenes);
        Enrichment en = new Enrichment(numGenesTotal, deGenes);

        Result result = new Result();
        //alternative
        for (int i = 0; i < 1; i++) { //LATER 1000
            Functions.sample_genes(new HashSet<>(r.geneMap.values()), 0.2);
            result.gather_runs(en.enrich(new HashSet<>(r.geneMap.values()), gos));  //fixme change method call
        }

        //alternative
        //will do only if percent > 0,2
        double percent = Functions.extend_flex_set(new HashSet<>(r.geneMap.values()));
        System.out.println(percent);
        for (int i = 0; i < 1; i++) { //LATER 1000
            Functions.sample_genes(new HashSet<>(r.geneMap.values()), percent);
            result.gather_runs(en.enrich(new HashSet<>(r.geneMap.values()), gos));
        }

        if (params.has("out")) { //FIXME probably always true because of default --> fix it
            String filepath = (String) params.valueOf("out");
            result.writeRobustGOs(result.getXquantileGOnodes(0.95), filepath);
        } else {
            result.printRobustGOs(result.getXquantileGOnodes(0.95));
        }

    }

}
