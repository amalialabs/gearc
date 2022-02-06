public class Gene {

    String gene_id;
    double weighted_score;
    double fdr;
    double fc;
    boolean is_significant;

    public Gene(String gene_id, double fc, double fdr, boolean is_significant) {
        this.gene_id = gene_id;
        this.fc = fc;
        this.fdr = fdr;
        this.is_significant = is_significant;
    }
}
