import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

public class Result {

    public HashMap<String, ArrayList<Double>> GOnode2FDRruns = new HashMap<>();

    /*
    collects FDR values for every 1000 runs to each GO node
    @param map of GO node to FDR value
     */
    public void gather_runs(HashMap<String, Double> GOnode2FDR) {
        GOnode2FDR.keySet().forEach(_node -> {
            if (GOnode2FDRruns.containsKey(_node)) {
                GOnode2FDRruns.get(_node).add(GOnode2FDR.get(_node));
            } else {
                ArrayList<Double> tmp = new ArrayList<>();
                tmp.add(GOnode2FDR.get(_node));
                GOnode2FDRruns.put(_node, tmp);
            }
        });
    }

    /*
    sorts FDR values ascending for every GO node
    @param map of GO nodes to lists of FDR values
    @return same map but with sorted FDR values
     */
    public HashMap<String, ArrayList<Double>> sortFDR(HashMap<String, ArrayList<Double>> GOnode2fdrs) {
        GOnode2fdrs.keySet().forEach(_node -> Collections.sort(GOnode2fdrs.get(_node)));
        return(GOnode2fdrs);
    }

    /*
    selects GO nodes which are significant in 95% (quantile) of 1000 runs
    @param specified quantile in double format eg 0.95
    @return list of robust GO nodes to specified quantile
     */
    public void getXquantileGOnodes(double quantile) {
        //TODO
    }

}
