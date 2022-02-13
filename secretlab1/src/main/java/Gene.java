public class Gene {

    String gene_id;
    double weighted_score;
    double fdr;
    double fc;
    boolean is_significant;
    public enum corresponding_set {SIG_CORE, FLEX, SIGNON_CORE};
    corresponding_set set;

    public Gene(String gene_id, double fc, double fdr, boolean is_significant) {
        this.gene_id = gene_id;
        this.fc = fc;
        this.fdr = fdr;
        this.is_significant = is_significant;
    }

    public double get_FDR_value() {
        return this.fdr;
    }

    public double get_FC_value() {
        return this.fc;
    }

    public double get_weighted_score() {
        return this.weighted_score;
    }

}
