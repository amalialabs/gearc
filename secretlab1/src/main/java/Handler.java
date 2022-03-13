import joptsimple.OptionParser;
import joptsimple.OptionSet;
import java.io.File;
import java.io.IOException;
import java.util.*;

import static java.util.Arrays.asList;
import static org.junit.Assert.*;

public class Handler {

    public enum expected_change {LOW, AVERAGE, HIGH};

    //todo parse arguments
    //todo insert all optional args
    public static void main(String[] args) {
        OptionParser optionParser = new OptionParser() {{
        accepts("genelist").withRequiredArg().required().ofType(File.class).describedAs("diffexp output");
        accepts("expectedChange").withRequiredArg().ofType(expected_change.class).defaultsTo(expected_change.AVERAGE).describedAs("expected change in diffexp data");
        accepts("FC").withRequiredArg().ofType(double.class).defaultsTo(1.0).describedAs("log2foldchange cutoff");
        accepts("pval").withRequiredArg().ofType(double.class).defaultsTo(0.01).describedAs("adjusted p-value cutoff");
        acceptsAll( asList( "h", "?" ), "show help" ).forHelp();
        accepts("out").withRequiredArg().ofType(String.class).describedAs("output file path"); //LATER defaults to installation path
        }};

        OptionSet params = optionParser.parse(args);
        assertTrue(params.has("genelist"));

        File genelist = (File) params.valueOf("genelist");
        try {
            optionParser.printHelpOn(System.out);
        } catch (IOException e) {
            new RuntimeException("Could not display help page.", e);
        }

        //TODO this all moved to Reader as part of constructor
//        //read Input, init gene object thereby -> hashset<gene>
        Reader r = new Reader();
        GO gos = null;

        //define_FDR_and_FC_cutoff_interval()
        // now also performed directly while reading the expression input in the reader
        //        Functions.filter_unclear(set).forEach(_g -> genes.put(_g.gene_id, _g));
        /**
         * @see Reader#readExpressionFile(File, boolean, HashMap)
         */


        // readmappingGAF after reduction in names
        // --> done in reader constructor now

        //score_genes()
        Functions.score_genes(new HashSet<>(r.geneMap.values()));

        //assign_genes_to_sets()
        //later -> make expected change as optional user param
        Functions.assign_genes_to_sets(new HashSet<>(r.geneMap.values()), expected_change.AVERAGE);


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
        for (int i = 0; i < 1; i++) { //later
            Functions.sample_genes(new HashSet<>(r.geneMap.values()), 0.2);
            result.gather_runs(en.enrich(new HashSet<>(r.geneMap.values()), gos));  //fixme change method call
        }


        //alternative
        //will do only if percent > 0,2
        double percent = Functions.extend_flex_set(new HashSet<>(r.geneMap.values()));
        System.out.println(percent);
        for (int i = 0; i < 1; i++) { //later
            Functions.sample_genes(new HashSet<>(r.geneMap.values()), percent);
            result.gather_runs(en.enrich(new HashSet<>(r.geneMap.values()), gos));
        }
        //later add user param
        if (params.has("out")) {
            String filepath = (String) params.valueOf("out"); //FIXME file name or just path???
            result.writeRobustGOs(result.getXquantileGOnodes(0.95), filepath);
        } else {
            result.printRobustGOs(result.getXquantileGOnodes(0.95));
        }

        //plots

    }

}
