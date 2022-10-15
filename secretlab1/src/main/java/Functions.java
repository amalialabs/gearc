import org.apache.commons.math3.analysis.function.Gaussian;

import java.text.DecimalFormat;
import java.util.*;
import java.util.stream.Collectors;

public class Functions {

    static DecimalFormat df = new DecimalFormat("#.###");

    /*
    calculates FDR and FC interval based on input distribution
    takes minimum of #genes to extend upper bound (either from sig genes up or from non-sig genes down),
    then centers lower bound around CUTOFF
    @param set of gene objects (gene, FC, FDR)
    @return double[] interval with lower and upper bound
     */
    public static double[] define_FDR_and_FC_cutoff_interval(Set<Gene> gene2FDRandFC, String type, double FDR_cutoff, double FC_cutoff) {
        Set<Gene> tmp = gene2FDRandFC.stream().filter(_g -> !_g.unclear).collect(Collectors.toSet());
        int num_sig_genes;
        int num_nonsig_genes;
        if (type.equals("FDR")) {
            num_sig_genes = (int) tmp.stream().filter(_gene -> !_gene.unclear && _gene.fdr<=FDR_cutoff).count();
            num_nonsig_genes = (int) tmp.stream().filter(_gene -> !_gene.unclear && _gene.fdr>FDR_cutoff).count();
        } else {
            num_sig_genes = (int) tmp.stream().filter(_gene -> !_gene.unclear && Math.abs(_gene.fc)>=FC_cutoff).count();
            num_nonsig_genes = (int) tmp.stream().filter(_gene -> !_gene.unclear && Math.abs(_gene.fc)<FC_cutoff).count();
        }
        int num_one_percent_ontop = Math.max((int) (0.01 * num_sig_genes), 5); //LATER at least 5? more genes
        int num_five_percent = (int) (0.05 * num_nonsig_genes);
        ArrayList<Gene> sorted_genes;
        if (type.equals("FDR")) {
            sorted_genes = (ArrayList<Gene>) tmp.stream().
                    sorted(Comparator.comparing(Gene::get_FDR_value)).collect(Collectors.toList());
        } else { //FC
            sorted_genes = (ArrayList<Gene>) tmp.stream().sorted(Comparator.comparing(Gene::get_abs_FC_value)).
                            collect(Collectors.toList());
            Collections.reverse(sorted_genes);
        }
        int idx_last_sig_gene = num_sig_genes-1;
        int idx_extended = idx_last_sig_gene + Math.min(num_one_percent_ontop , num_five_percent);

        double diff_to_center = type.equals("FDR") ? Math.abs(FDR_cutoff-sorted_genes.get(idx_extended).fdr) : Math.abs(FC_cutoff-Math.abs(sorted_genes.get(idx_extended).fc)); //centered interval around cutoff
        double upper_bound = type.equals("FDR") ? sorted_genes.get(idx_extended).fdr : FC_cutoff+diff_to_center;
        double lower_bound = type.equals("FDR") ? FDR_cutoff-diff_to_center : Math.abs(sorted_genes.get(idx_extended).fc);
        if (type.equals("FDR")) {
            assert lower_bound > 0.0 && upper_bound <= 1.0 && upper_bound > FDR_cutoff && lower_bound < FDR_cutoff;
        } else {
            assert lower_bound > 0.0 && upper_bound > FC_cutoff && lower_bound < FC_cutoff;
        }
        String temp = df.format(lower_bound);
        temp = temp.replace(",", ".");
        lower_bound = Double.parseDouble(temp);
        temp = df.format(upper_bound);
        temp = temp.replace(",", ".");
        upper_bound = Double.parseDouble(temp);
        double[] interval = {lower_bound, upper_bound};
        return interval;
    }

    /*
    filters unclear genes, leaving only significantly regulated or non-regulated genes
    @param set of gene object (gene, FC, FDR)
    @return same set without entries of unclear genes
     */
    public static Set<Gene> filter_unclear(Set<Gene> gene2FCandFDR, double FDR_cutoff, double FC_cutoff) {
        Set<Gene> sig_genes = gene2FCandFDR.stream().
                filter(gene -> gene.fdr <= FDR_cutoff && Math.abs(gene.fc) >= FC_cutoff).
                collect(Collectors.toSet());
        sig_genes.addAll(gene2FCandFDR.stream().filter(_gene -> _gene.fdr > FDR_cutoff &&
                Math.abs(_gene.fc) < FC_cutoff).collect(Collectors.toSet())); // Gene wo nur 1 sig ist von FDR/FC fliegen raus
        //LATER Gergely: FDR(Gene in [-1,1])<=0.05
        return sig_genes;
    }

    /*
    computes weighted score per gene based on FDR and FC
    @param set of gene object (gene, FC, FDR)
    @return same set with now computed values
    LATER remove return??
     */
    public static Set<Gene> score_genes(Set<Gene> gene2FCandFDR, double FDR_cutoff, double FC_cutoff) {
        double[] FDR_interval = define_FDR_and_FC_cutoff_interval(gene2FCandFDR, "FDR", FDR_cutoff, FC_cutoff);
        double fdr_interval_width = FDR_interval[1]-FDR_interval[0];
        double[] FC_interval = define_FDR_and_FC_cutoff_interval(gene2FCandFDR, "FC", FDR_cutoff, FC_cutoff);
        double fc_interval_width = FC_interval[1]-FC_interval[0];
        System.out.println("fdr   " + FDR_interval[0] + "   " + FDR_interval[1]);
        System.out.println("fc    " + FC_interval[0] + "    " + FC_interval[1]);
        gene2FCandFDR.forEach(_obj -> {
            if (_obj.fdr <= FDR_interval[1] && Math.abs(_obj.fc) >= FC_interval[0]) {
                double fdr_score = _obj.fdr<=FDR_interval[0] ? 1.0 : (_obj.fdr - FDR_interval[0]) / fdr_interval_width;
                double fc_score = Math.abs(_obj.fc)>=FC_interval[1] ? 1.0 : (Math.abs(_obj.fc) - FC_interval[0]) / fc_interval_width;
                _obj.weighted_score = 0.5*(fdr_score+fc_score);
            } else {
                _obj.weighted_score = 0.0;
            }
        });
        return gene2FCandFDR;
    }

    /*
    assign genes into core-sig, flex. core-sig_nonreg according to expected change
    @param set of gene objects (gene, weighted_score), enum expected change (high, average, low)
    @return set of gene objects with defined set belonging
     */
    public static Set<Gene> assign_genes_to_sets(Set<Gene> gene2weightedscore, Enum expected_change) {
        ArrayList<Gene> sorted_genes = (ArrayList<Gene>) gene2weightedscore.stream().
                sorted(Comparator.comparing(Gene::get_weighted_score).reversed()).collect(Collectors.toList()); //TODO chech if ascending
        int sig_core;
        int flex;
        if (expected_change == Handler.expected_change.AVERAGE) { //30%-40%-30%
            sig_core = (int) (0.3 * gene2weightedscore.size());
            flex = (int) (0.4 * gene2weightedscore.size());
        } else if (expected_change == Handler.expected_change.HIGH) { //50%-40%-10%
            sig_core = (int) (0.5 * gene2weightedscore.size());
            flex = (int) (0.4 * gene2weightedscore.size());
        } else { //10%-40%-50%
            sig_core = (int) (0.1 * gene2weightedscore.size());
            flex = (int) (0.5 * gene2weightedscore.size());
        }
        sorted_genes.subList(0,sig_core).forEach(_gene -> _gene.set = Gene.corresponding_set.SIG_CORE);
        sorted_genes.subList(sig_core, sig_core+flex).forEach(_gene -> _gene.set = Gene.corresponding_set.FLEX);
        sorted_genes.subList(sig_core+flex, sorted_genes.size()).forEach(_gene -> _gene.set = Gene.corresponding_set.SIGNON_CORE);
        return new HashSet<>(sorted_genes);
    }

    /*
    extend proportion of flex genes as long as bonus >= penalty based on gauss distribution (x-axis) and weighted score (y-axis)
    @param set of gene objects (gene, weighted_score)
    @return percentage (optimally > 0.2) which could be tested
     */
    public static double extend_flex_set(HashSet<Gene> genes2sets) {
        ArrayList<Gene> sorted_genes = (ArrayList<Gene>) genes2sets.stream().
                sorted(Comparator.comparing(Gene::get_weighted_score)).collect(Collectors.toList());
        int idx_last_flex = sorted_genes.stream().
                filter(_g -> _g.set == Gene.corresponding_set.FLEX).
                collect(Collectors.toSet()).size()-1;
        int nsig = (int) sorted_genes.stream().filter(_g -> _g.set==Gene.corresponding_set.SIG_CORE).count();

        int gauss_mean = nsig;
        int SD = (int) Math.round(0.997*sorted_genes.size()/8); //because ND divided into 8 parts and we want nonsig to be 0
        Gaussian gauss = new Gaussian(gauss_mean, SD);

        double z = sorted_genes.get(idx_last_flex).weighted_score + 0.000001; //LATER can change to 0.01 and multiply gauss with 10000
        double bonus = 1.0;
        double penalty = 0.0;
        int idx_current_gene = idx_last_flex;
        while (bonus > penalty && idx_current_gene < sorted_genes.size()-1) {
            idx_current_gene = idx_current_gene + 1;
            bonus = gauss.value(idx_current_gene);
            penalty = Math.abs(sorted_genes.get(idx_current_gene).weighted_score-z) * (idx_current_gene-idx_last_flex);
        }
        int idx_profitable_extension = idx_current_gene - 1;
        double new_percentage = (idx_profitable_extension*1.0) / sorted_genes.size();
        return(idx_profitable_extension > idx_last_flex ? new_percentage : 0.2);  //if not at all profitable return default 20%
    }

    /*
    sample genes from each of 3 sets by 70%-20%-90%
    @param set of gene objects (gene,corresponding_set)
    @return random sampled subset of gene objects
     */
    public static Set<Gene> sample_genes(Set<Gene> gene2set, double flex_percentage) {
        HashSet<Gene> sampled_genes = new HashSet<>();
        Set<Gene> sig_core_genes = gene2set.stream().
                filter(_gene -> _gene.set == Gene.corresponding_set.SIG_CORE).collect(Collectors.toSet());
        ArrayList<Gene> sigs = new ArrayList<>(){{addAll(sig_core_genes);}};
        Collections.shuffle(sigs);
        int seventy_percent = (int) (0.7 * sigs.size());
        sampled_genes.addAll(sigs.subList(0, seventy_percent));
        Set<Gene> flex_genes = gene2set.stream().
                filter(_gene -> _gene.set == Gene.corresponding_set.FLEX).collect(Collectors.toSet());
        ArrayList<Gene> flexs = new ArrayList<>(){{addAll(flex_genes);}};
        Collections.shuffle(flexs);
        int flex_percent = (int) (flex_percentage * flexs.size());
        sampled_genes.addAll(flexs.subList(0, flex_percent));
        Set<Gene> signon_core_genes = gene2set.stream().
                filter(_gene -> _gene.set == Gene.corresponding_set.SIGNON_CORE).collect(Collectors.toSet());
        ArrayList<Gene> signons = new ArrayList<>(){{addAll(signon_core_genes);}};
        Collections.shuffle(signons);
        int ninety_percent = (int) (0.9 * signons.size());
        sampled_genes.addAll(signons.subList(0, ninety_percent));
        return sampled_genes;
    }

}
