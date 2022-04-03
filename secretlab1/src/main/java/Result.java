import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

public class Result {

    public HashMap<Node, ArrayList<Double>> GOnode2FDRruns = new HashMap<>();
    public HashMap<Node, ArrayList<Double>> GOnode2FDRrunsExtend = new HashMap<>();

    /*
    collects FDR values for every 1000 runs to each GO node, sorts FDRS each time
    @param map of GO node to FDR value
     */
    public void gather_runs(HashMap<Node, Double> GOnode2FDR, boolean extend) {
        if (extend) {
            GOnode2FDR.keySet().forEach(_node -> {
                if (GOnode2FDRrunsExtend.containsKey(_node)) {
                    GOnode2FDRrunsExtend.get(_node).add(GOnode2FDR.get(_node));
                } else {
                    ArrayList<Double> tmp = new ArrayList<>();
                    tmp.add(GOnode2FDR.get(_node));
                    GOnode2FDRrunsExtend.put(_node, tmp);
                }
            });
        } else {
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
        sortFDR();
    }

    /*
    sorts FDR values ascending for every GO node
     */
    public void sortFDR() {
        GOnode2FDRruns.keySet().forEach(_node -> Collections.sort(GOnode2FDRruns.get(_node)));
        GOnode2FDRrunsExtend.keySet().forEach(_node -> Collections.sort(GOnode2FDRrunsExtend.get(_node)));
    }

    /*
    selects GO nodes which are significant in 95% (quantile) of 1000 runs
    @param specified quantile in double format eg 0.95
    @return list of robust GO nodes to specified quantile
     */
    public Set<Node> getXquantileGOnodes(double quantile) {
        double FDR_cutoff = 0.05; //LATER may be defined by user too
        Node node = GOnode2FDRruns.entrySet().iterator().next().getKey();
        int num_runs = GOnode2FDRruns.get(node).size();
        int position = (int) quantile*num_runs;
        Set<Node> robustGOs = GOnode2FDRruns.keySet().stream().
                filter(_go -> GOnode2FDRruns.get(_go).get(position)<=FDR_cutoff).
                collect(Collectors.toSet());
        return(robustGOs);
    }

    /*
    computes mean FDR value of all runs for given GO node
    @param GO node
    @return double mean FDR over all runs
     */
    public double getMeanFDRofGO(Node GO) {
        ArrayList<Double> fdrs = GOnode2FDRruns.get(GO);
        double sum = fdrs.stream().mapToDouble(a -> a).sum();
        return(sum/fdrs.size());
    }

    /*
    prints robust GOs of a specified quantile
    @param list of GOs of quantile
     */
    public void printRobustGOs(Set<Node> robustGOs) {
        //LATER output in html table -> need optional argument in handler
        System.out.println("GOnode\tmeanFDR");
        robustGOs.forEach(_go -> System.out.println(_go.node_id + "\t" + getMeanFDRofGO(_go)));
        System.out.println(robustGOs.size() + " GOs in total.");
    }

    public void writeRobustGOs(Set<Node> robustGOs, String outdir) {
        try {
            BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outdir, "/robust_GOs.tsv")));
            bw.write("GOnode\tmeanFDR\n");
            for (Node n : robustGOs) {
                bw.write(n.node_id + "\t" + getMeanFDRofGO(n) + "\n");
            }
            bw.close();
        } catch (IOException e) {
            throw new RuntimeException("could not init file ", e);
        }
    }

}
