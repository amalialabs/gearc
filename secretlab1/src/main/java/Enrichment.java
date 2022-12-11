import org.apache.commons.math3.distribution.HypergeometricDistribution;
import java.util.*;
import java.util.stream.Collectors;

public class Enrichment {

    int numGenesTotal;
    int deGenes;

    public Enrichment(int numGenesTotal, int deGenes) {
        this.numGenesTotal = numGenesTotal;
        this.deGenes = deGenes;
    }

    public HashMap<Node, Double> enrich(Set<Gene> sampled_genes, GO gos) {
        HashSet<Gene> foreground_genes = new HashSet<>();
        HashSet<Gene> background_genes = new HashSet<>();
        sampled_genes.stream().forEach(_gene -> {
            if (_gene.is_significant) {
                foreground_genes.add(_gene);
            } else {
                background_genes.add(_gene);
            }
        });

        HashMap<Node, Double> node2FDR = new HashMap<>();

        // LATER vielleicht beschränkt man die nodes bzgl Anzahl Gene also nodes mit #genes > x und < y
        for (Node node : gos.getGoNodes().values()) {
            double pvalue = 1.0;
            if (node.getGenes() != null && node.getRelevantGenesSize() != 0) {
                pvalue = hypergeometric(node, sampled_genes);
            }
            node2FDR.put(node, pvalue);
        }
        bhAdjust(node2FDR);
        return(node2FDR);
    }

    public double hypergeometric(Node node, Set<Gene> sampled_genes) {
        int numMeasuredGenesInSet = (int) node.getGenes().stream().filter(_g -> sampled_genes.contains(_g) && !_g.unclear).count();
        int numSignificantGenesInSet = (int) node.getGenes().stream().filter(_g -> _g != null && _g.is_significant && sampled_genes.contains(_g)).count();
        double hg_pval = 1.0;
        if (numMeasuredGenesInSet != 0) {
            hg_pval = new HypergeometricDistribution(numGenesTotal, deGenes,
                    numMeasuredGenesInSet).upperCumulativeProbability(numSignificantGenesInSet);
        }
        return(hg_pval);
    }

    public void bhAdjust(HashMap<Node, Double> pval) {
        List<Map.Entry<Node, Double>> list = pval.entrySet().stream()
                .sorted(Comparator.comparingDouble(Map.Entry::getValue)).collect(Collectors.toList());
        for (int i = list.size() - 1; i >= 0; i--) {
            Node node = list.get(i).getKey();
            if (i == list.size() - 1) {
                node.setBhFDR(pval.get(node));
            } else {
                double unadjust = pval.get(node);
                double left = pval.get(list.get(i + 1).getKey());
                double right = (list.size() / (i + 1.0)) * unadjust;
                node.setBhFDR(Math.min(left, right));
            }
        }
    }
}
