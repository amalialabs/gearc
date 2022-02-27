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
        }};

        OptionSet params = optionParser.parse(args);
        assertTrue(params.has("genelist"));

        File genelist = (File) params.valueOf("genelist");
        //TODO
        //read Input, init gene object thereby -> hashset<gene>
        Reader r = new Reader();
        r.readExpressionFile(new File("/home/birinci/GOEnrichment/simul_exp_go_bp_ensembl.tsv"));
        GO gos = r.readOboFile(new File("/home/birinci/GOEnrichment/go.obo"), "biological_process");

        //define_FDR_and_FC_cutoff_interval()
        HashMap<String, Gene> genes = new HashMap<>();
        Set<Gene> set = r.readExpressionFile(new File("/home/birinci/GOEnrichment/simul_exp_go_bp_ensembl.tsv"));
        System.out.println(set.size());
        Functions.filter_unclear(set).forEach(_g -> genes.put(_g.gene_id, _g));
        r.geneMap = genes;
        System.out.println(genes.values().size());

        // readmappingGAF after reduction in names
        r.readMappringGAF(new File("/home/birinci/GOEnrichment/goa_human.gaf.gz"), gos);
        gos.getGoNodes().values().forEach(_node -> {
            if(_node.getGenes() != null && _node.getGenes().size() > 2) {
                System.out.println(_node.node_id);
                System.out.println(_node.getNode_name());
                System.out.println(_node.getGenes().size());
            }
        });

        //score_genes()
        Functions.score_genes(new HashSet<>(r.geneMap.values()));

        //assign_genes_to_sets()
        //later -> make expected change as optional user param
        Functions.assign_genes_to_sets(new HashSet<>(r.geneMap.values()), expected_change.AVERAGE);


        int numGenesTotal = gos.getGoNodes().values().stream()
                .filter(_node -> _node.getGenes() != null).mapToInt(_node -> _node.getGenes().size()).sum();
        //todo need to remove "null"genes from Nodes
        Set<Gene> ge = new HashSet<>();
        gos.getGoNodes().values().stream().forEach(_node -> {
            if (_node.getGenes() != null) ge.addAll(_node.getGenes());
        });
        int deGenes = (int) ge.stream().filter(_go -> _go.is_significant).count();
        System.out.println(numGenesTotal + "\t" + deGenes);
        Enrichment en = new Enrichment(numGenesTotal, deGenes);

        Result result = new Result();
        //alternative
        for (int i = 0; i < 1; i++) { //later
            Functions.sample_genes(new HashSet<>(r.geneMap.values()), 0.2);
            result.gather_runs(en.enrich(new HashSet<>(r.geneMap.values()), gos));
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
        result.writeRobustGOs(result.getXquantileGOnodes(0.95));
        //plots

        try {
            optionParser.printHelpOn(System.out);
        } catch (IOException e) {
            new RuntimeException("Could not display help page.", e);
        }
    }

}
