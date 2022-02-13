import java.util.HashSet;

public class Enrichment {

    public void enrich(HashSet<Gene> sampled_genes) {
        HashSet<Gene> foreground_genes = new HashSet<>();
        HashSet<Gene> background_genes = new HashSet<>();
        sampled_genes.stream().forEach(_gene -> {
            if (_gene.is_significant) {
                foreground_genes.add(_gene);
            } else {
                background_genes.add(_gene);
            }
        });

        //TODO start enrichment
        //gather GO node - FDR as result
    }

}
