import org.apache.commons.math3.distribution.HypergeometricDistribution;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;

public class Enrichment {

    //fixme need to be set before calling
    int numGenesTotal;
    int deGenes;

    public HashMap<Node, Double> enrich(HashSet<Gene> sampled_genes, GO gos) {
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
        //
        // lATER vielleicht beschränkt man die nodes bzgl Anzahl Gene also nodes mit #genes > x und < y
        for (Node node : gos.getGoNodes().values()) {
            double pvalue = hypergeometric(node);
            node2FDR.put(node, pvalue);
        }
        benjamini(node2FDR); // sollte die FDRs überschreiben
        return(node2FDR);
    }

    public double hypergeometric(Node node) {
        int numMeasuredGenesInSet = node.getGenes().size();
        int numSignificantGenesInSet = (int) node.getGenes().stream().filter(_g -> _g.is_significant).count();
        double hg_pval = new HypergeometricDistribution(numGenesTotal, deGenes,
                numMeasuredGenesInSet).upperCumulativeProbability(numSignificantGenesInSet);
        return(hg_pval);
    }

    public void benjamini(HashMap<Node, Double> pval) {
        Double[] pvalues = new Double[pval.size()];
        String[] keys = new String[pval.size()];
        int counter = 0;
        for (Map.Entry<Node, Double> entry : pval.entrySet()) {
            pvalues[counter] = entry.getValue();
            keys[counter] = entry.getKey().node_id;
            counter++;
        }
        //FIXME entweder eigene BH Klasse oder nur eine kleine Funktion hier
//        BenjaminiHochbergFDR bh = new BenjaminiHochbergFDR(pvalues);
//        int[] index = bh.createIndex(pvalues);
//        bh.calculate();
//        Double[] result = bh.getAdjustedPvalues();
//        for (int i = 0; i < result.length; i++) {
//            benjaminipvalue.put(keys[index[i]], result[i]);
//        }
    }

}
