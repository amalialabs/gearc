import java.util.HashSet;
import java.util.stream.Collectors;

public class Functions {

    public double FDR_cutoff = 0.05;
    public double FC_cutoff = 1.0;

    /*
    @param hashmap og gene <-> FC, FDR
    @return same hashmap without entries of unclear genes
     */
    public void filter_unclear(HashSet<Gene> gene2FCandFDR) {
        gene2FCandFDR.stream().
                filter(gene -> gene.fdr <= FDR_cutoff && (gene.fc >= FC_cutoff || gene.fc <= -FC_cutoff)).
                collect(Collectors.toSet());
        // TODO hier fehlt noch das von Gergely zum abschätzen der FDR von FC element von (a,b)
    }

}
