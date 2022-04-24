import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.*;
import java.util.stream.Collectors;

public class Result {

    public HashMap<Node, ArrayList<Double>> GOnode2FDRruns = new HashMap<>();
    public HashMap<Node, ArrayList<Double>> GOnode2FDRrunsExtend = new HashMap<>();

    public double FDR_cutoff;

    DecimalFormat df = new DecimalFormat("0.####E0");

    public Result (double fdr) {
        this.FDR_cutoff = fdr;
    }

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
        //LATER output in html table
        System.out.println("GOnode\tmeanFDR");
        robustGOs.forEach(_go -> System.out.println(_go.node_id + "\t" + df.format(getMeanFDRofGO(_go))));
        System.out.println(robustGOs.size() + " GOs in total.");
    }

    public void writeRobustGOs(Set<Node> robustGOs, String outdir) {
        String file = outdir + "/robust_GOs.tsv";
        System.out.println(file);
        File f = new File(file);
        BufferedWriter bw;
        try {
            bw = new BufferedWriter(new FileWriter(f));
            bw.write("GOnode\tmeanFDR\n");
            for (Node n : robustGOs) {
                bw.write(n.node_id + "\t" + df.format(getMeanFDRofGO(n)) + "\n");
            }
            bw.close();
            System.out.println("normally printed robust_GOs.tsv file " + f.getAbsolutePath());
        } catch (IOException e) {
            throw new RuntimeException("could not init file ", e);
        }
    }

    public void printStandardGOs(HashMap<Node, Double> gos) {
        gos = gos.entrySet().stream()
                .sorted(Map.Entry.comparingByValue())
                .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue,
                        (e1, e2) -> e1, LinkedHashMap::new));
        HashMap<Node, Double> standardGOs = gos;
        System.out.println("GOnode\tFDR");
        int sig = (int) standardGOs.keySet().stream().filter(_go -> standardGOs.get(_go) <= this.FDR_cutoff).count();
        standardGOs.keySet().forEach(_go -> {
            if (standardGOs.get(_go) <= this.FDR_cutoff) {
                System.out.println(_go.node_id + "\t" + df.format(standardGOs.get(_go)));
            }
        }
        );
        System.out.println(sig + "/" + gos.size() + " enriched GOs");
    }

    public void writeStandardGOs(HashMap<Node, Double> standardGOs, String outdir) {
        standardGOs = standardGOs.entrySet().stream()
                .sorted(Map.Entry.comparingByValue())
                .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue,
                        (e1, e2) -> e1, LinkedHashMap::new));
        try {
            BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outdir, "standard_GOs.tsv")));
            bw.write("GOnode\tFDR\n");
            for (Node n : standardGOs.keySet()) {
                bw.write(n.node_id + "\t" + df.format(standardGOs.get(n)) + "\n");
            }
            bw.close();
        } catch (IOException e) {
            throw new RuntimeException("could not init file ", e);
        }
    }

}
