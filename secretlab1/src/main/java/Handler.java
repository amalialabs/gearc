import joptsimple.OptionParser;
import joptsimple.OptionSet;

import java.io.File;

import static org.junit.Assert.*;

public class Handler {

    public static void main(String[] args) {
        OptionParser optionParser = new OptionParser();
        optionParser.accepts("genelist").withRequiredArg().ofType(File.class);

        OptionSet options = optionParser.parse("-genelist");
        System.out.println(options.has("genelist"));

        optionParser.accepts("expectedChange").withOptionalArg(); //of type enum {high, average, low}
        optionParser.accepts("FC").withOptionalArg().ofType(Double.class);
        optionParser.accepts("pval").withOptionalArg().ofType(Double.class);

        System.out.println(options.valueOf("genelist"));
        File genelist = (File) options.valueOf("genelist");
        //print help for calling cmd
    }

}
