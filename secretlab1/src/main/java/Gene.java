public class Gene {

    String gene_id;
    String gene_name;
    double weighted_score;
    double fdr;
    double fc;
    boolean is_significant;
    boolean unclear;
    boolean not_signif;
    public enum corresponding_set {SIG_CORE, FLEX, SIGNON_CORE};
    corresponding_set set;

    public Gene(String gene_id, double fc, double fdr) {
        this.gene_id = gene_id;
        this.fc = fc;
        this.fdr = fdr;
        this.is_significant = Math.abs(fc) >= Functions.FC_cutoff && fdr <= Functions.FDR_cutoff;
        this.not_signif = Math.abs(fc) < Functions.FC_cutoff && fdr > Functions.FDR_cutoff;
        this.unclear = !is_significant && !not_signif;
    }

    public Gene(String gene_id, String gene_name, double fc, double fdr) {
        this.gene_id = gene_id;
        this.gene_name = gene_name;
        this.fdr = fdr;
        this.fc = fc;
        this.is_significant = Math.abs(fc) >= Functions.FC_cutoff && fdr <= Functions.FDR_cutoff;
        this.not_signif = Math.abs(fc) < Functions.FC_cutoff && fdr > Functions.FDR_cutoff;
        this.unclear = !is_significant && !not_signif;
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

    @Override
    public String toString() {
        return "Gene{" +
                "gene_id='" + gene_id + '\'' +
                ", weighted_score=" + weighted_score +
                ", fdr=" + fdr +
                ", fc=" + fc +
                ", is_significant=" + is_significant +
                ", set=" + set +
                '}';
    }
}
