import joptsimple.OptionParser;
import joptsimple.OptionSet;

import java.io.File;
import java.util.ArrayList;


public class Handler {

    public static void main(String[] args) {
        //parse parameters#
        //print help for calling cmd
        OptionParser optionParser = new OptionParser("");
        OptionSet options = optionParser.parse("-genelist");

        optionParser.accepts("genelist").withRequiredArg();

        optionParser.accepts("expectedChange").withOptionalArg();
        optionParser.accepts("FC").withOptionalArg();
        optionParser.accepts("pval").withOptionalArg();


    }

}
