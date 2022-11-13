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

//    DecimalFormat df = new DecimalFormat("0.####E0");
    public String format(double num) {
        return String.format("%.3e", num);
    }

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
        Node node = GOnode2FDRruns.entrySet().iterator().next().getKey();
        int num_runs = GOnode2FDRruns.get(node).size();
        int min_num_sig = (int) Math.round(num_runs * quantile);
        Set<Node> robustGOs = GOnode2FDRruns.keySet().stream().
                filter(_go -> GOnode2FDRruns.get(_go).stream().filter(_r -> _r<= this.FDR_cutoff).count()>=min_num_sig).
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

    public double getWorstFDRofGO(Node GO, double quantile) {
        ArrayList<Double> fdrs = GOnode2FDRruns.get(GO);
        Collections.sort(fdrs);
        double worst_quantile_fdr = fdrs.get((int) Math.round(quantile*fdrs.size()));
        return(worst_quantile_fdr);
    }

    /*
    prints robust GOs of a specified quantile
    @param list of GOs of quantile
     */
    public void printRobustGOs(Set<Node> robustGOs) {
        //LATER output in html table
        System.out.println("GOnode\tmeanFDR");
        robustGOs.forEach(_go -> System.out.println(_go.node_id + "\t" + format(getMeanFDRofGO(_go))));
        System.out.println(robustGOs.size() + " GOs in total.");
    }

    public void writeRobustGOs(Set<Node> robustGOs, String outdir, double quantile) {
        File f = new File(outdir + File.separator + "robust_GOs_worstQuantileFDR.tsv");
        BufferedWriter bw;
        try {
            bw = new BufferedWriter(new FileWriter(f));
            bw.write("GOnode\tquantileFDR\n");
            for (Node n : robustGOs) {
                bw.write(n.node_id + "\t" + format(getWorstFDRofGO(n, quantile)) + "\n");
            }
            bw.close();
        } catch (IOException e) {
            throw new RuntimeException("could not init file ", e);
        }
    }

    public void writeMeanGOs(Set<Node> robustGOs, String outdir) {
        File f = new File(outdir + File.separator + "robust_GOs_meanFDR.tsv");
        BufferedWriter bw;
        try {
            bw = new BufferedWriter(new FileWriter(f));
            bw.write("GOnode\tmeanFDR\n");
            for (Node n : robustGOs) {
                bw.write(n.node_id + "\t" + format(getMeanFDRofGO(n)) + "\n");
            }
            bw.close();
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
                System.out.println(_go.node_id + "\t" + format(standardGOs.get(_go)));
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
        File f = new File(outdir + File.separator + "standard_GOs.tsv");
        try {
            BufferedWriter bw = new BufferedWriter(new FileWriter(f));
            bw.write("GOnode\tFDR\n");
            for (Node n : standardGOs.keySet()) {
                bw.write(n.node_id + "\t" + format(standardGOs.get(n)) + "\n");
            }
            bw.close();
        } catch (IOException e) {
            throw new RuntimeException("could not init file ", e);
        }
    }

    public void writeRankDifferences(HashMap<Node, Double> standardGOs, Set<Node> robustGO, String outDir) {
        Map<Node, Double> robustGOs = robustGO.stream().collect(Collectors.toMap(Node::getSelf, Node::getBhFDR));

        List<Node> standard = standardGOs.entrySet().stream().sorted(Comparator.comparingDouble(Map.Entry::getValue))
                .map(Map.Entry::getKey).collect(Collectors.toList());
        List<Node> robust = robustGOs.entrySet().stream().sorted(Comparator.comparingDouble(Map.Entry::getValue))
                .map(Map.Entry::getKey).collect(Collectors.toList());
        Set<Node> joint = new HashSet<>();
        joint.addAll(standard);
        joint.retainAll(robust);

        File f = new File(outDir, "GO_rank_diff.tsv");
        try {
            BufferedWriter bw = new BufferedWriter(new FileWriter(f));
            bw.write("GOnode\tstandard_FDR\trobust_FDR\tstandard_Rank\trobust_Rank\trank_diff\n");
            for (Node n : joint) {
                bw.write(n.node_id + "\t" + format(standardGOs.get(n)) + "\t" + format(robustGOs.get(n))
                        + "\t" + standard.indexOf(n) + "\t" + robust.indexOf(n) + "\t" + (standard.indexOf(n) - robust.indexOf(n)) + "\n");
            }
            bw.close();
        } catch (IOException e) {
            throw new RuntimeException("could not init file ", e);
        }
    }
}
