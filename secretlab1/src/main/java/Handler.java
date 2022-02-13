import joptsimple.OptionParser;
import joptsimple.OptionSet;
import java.io.File;
import java.io.IOException;
import static java.util.Arrays.asList;
import static org.junit.Assert.*;

public class Handler {

    public enum expected_change {LOW, AVERAGE, HIGH};

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
        //define_FDR_and_FC_cutoff_interval()
        //filter_unclear()
        //score_genes()
        //assign_genes_to_sets()
        // for (i in 1:1000)
        //      sample_genes(0.2)
        //      enrich()
        // percent = extend_flex_set()
        // do following only if percent > 0.2
        // for (i in 1:1000)
        //      sample_genes(percent)
        //      enrich()
        //gather enrich results in table
        //plots

        try {
            optionParser.printHelpOn(System.out);
        } catch (IOException e) {
            new RuntimeException("Could not display help page.", e);
        }
    }

}
