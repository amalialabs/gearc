import java.io.*;
import java.text.DecimalFormat;
import java.util.*;
import java.util.stream.Collectors;

public class Plots {

    double FDR_cutoff;
    double FC_cutoff;
    String out_dir;
    String ranktable_path;
    String genetable_path;
    String gos;
    String gotable_path;
    String gotable_extend_path;
    String gotable_standard;
    Result res;
    HashMap<Node, Double> gos_standard;

    String root;
    List<String> plotList;

    DecimalFormat df = new DecimalFormat("0.####E0");

    public Plots(String output_dir, Collection<Gene> genes, Collection<Node> gos, double fdr, double fc, Result res, HashMap<Node, Double> gos_standard) {
        this.FDR_cutoff = fdr;
        this.FC_cutoff = fc;
        this.out_dir = output_dir;
        this.genetable_path = out_dir + File.separator + "genes.table";
        this.gos = output_dir + File.separator + "gos.table";
        this.gotable_path = output_dir + File.separator + "gos2fdrs.table";
        this.gotable_extend_path = output_dir + File.separator + "gos2fdrs_extended.table";
        this.gotable_standard = output_dir + File.separator + "gos2fdrs_standard.table";
        this.ranktable_path = output_dir + File.separator + "rank_differences.table";  //FIXME Armin hier brauch ich die table
        this.res = res;
        this.gos_standard = gos_standard;

        this.root = new File("/secretlab1/src/main/resources/rscripts").exists() ? "/secretlab1/src/main/resources/rscripts" :
                getClass().getClassLoader().getResource("/rscripts/").toString();
        System.out.println(root);

        try {
            root = output_dir;
            plotList = Arrays.asList(//TODO replace with resource file scanner
                    "unclear_genes_BARPLOT.R",
                    "sig_genes_VOLCANO.R",
                    "gene_categories_BARPLOT.R",
                    "genes_scored_CUMULATIVE.R",
                    "genes_sets_PIECHART.R",
                    "gos_fdr_mean_quantile_SCATTER.R",
                    "selected_GOs_FDR_BOXPLOT.R",
                    "selected_GOs_rob_vs_extend_BOXPLOT.R",
                    "gos_quantile_vs_mean_JITTER.R",
                    "gos_standard_vs_robust_vs_extended_BARPLOT.R",
                    "gos_standard_vs_robust_vs_extended_VENN.R",
                    "fdr_cutoff_finding_CUMULATIVE.R",
                    "fc_cutoff_finding_CUMULATIVE.R",
                    "flexset_extension_CURVE.R",
                    "num_sig_gos_comparison_BARPLOT.R",
                    "expected_change_distrib_BARPLOT.R");

//            preprocess();
        } catch (Exception e) {
            throw new RuntimeException("", e);
        }


        createGeneTable(genes);
        createGOTable(gos);
    }

    void runPlots(Reader r) {
        de_scores_rank();
        unclear_genes_BARPLOT(r.allGenes.values());
        sig_genes_VOLCANO();
        gene_categories_BARPLOT();
        weighted_genes_CUMULATIVE();
        genes_set_PIECHART();
        go_fdrs_mean_vs_quantile_SCATTER();
        selected_gos_fdr_distrib_BOXPLOT();
        selected_gos_rob_vs_extend_BOXPLOT();
        gos_quantile_vs_mean_fdr_JITTER();
        gos_standard_vs_robust_vs_extended_BARPLOT();
        gos_standard_vs_robust_vs_extended_VENN();
        fdr_cutoff_finding_CUMULATIVE();
        fc_cutoff_finding_CUMULATIVE();
        flexset_extension_CURVE();
        expected_change_BARPLOT();
        num_sig_gos_BARPLOT();
    }

    //TODO eval if scripts needed
    private void preprocess(){
        for (String filename : plotList) {
            try (PrintWriter pw = new PrintWriter(new File(out_dir, filename))) {
                String result = new BufferedReader(new InputStreamReader(getClass().getResourceAsStream("rscripts/" + filename)))
                        .lines().collect(Collectors.joining("\n"));
                pw.write(result + "\n");
            } catch (IOException e) {
                throw new RuntimeException("", e);
            }
        }
    }

    //LATER change all script params to reading automatically from /out/ in build so that all gene/go tables are not a param anymore
    public void createGeneTable(Collection<Gene> genes) {
        try {
            File dir = new File(this.out_dir);
            if (!dir.exists()) {
                dir.mkdirs();
            }
            BufferedWriter bw = new BufferedWriter(new FileWriter(new File(this.out_dir, File.separator + "genes.table")));
            bw.write("geneID\tFDR\tlog2FC\tgeneset\tweighted_score\tis_sig\tis_unclear\tnot_sig\n");
            for (Gene gene : genes) {
                bw.write(gene.gene_id + "\t" + df.format(gene.fdr) + "\t" + gene.fc + "\t" + gene.set +
                  "\t" + gene.weighted_score + "\t" + gene.is_significant + "\t" + gene.unclear + "\t" + gene.not_signif + "\n");
            }
            bw.close();
        } catch (IOException e) {
            throw new RuntimeException("could not write file. ", e);
        }
    }

    public void createGOTable(Collection<Node> gos) {
        try {
            File dir = new File(this.out_dir);
            if (!dir.exists()) {
                dir.mkdirs();
            }
            BufferedWriter bw = new BufferedWriter(new FileWriter(new File(this.out_dir, File.separator + "gos2fdrs.table")));
            for (Node go : gos) {
                bw.write(go.node_id + "\t" + go.node_name);
                for (double fdr : res.GOnode2FDRruns.get(go)) {
                    bw.write("\t" + df.format(fdr));
                }
                bw.write("\n");
            }
            bw.close();
            //TODO add later
//            bw = new BufferedWriter(new FileWriter(new File(this.out_dir, File.separator + "gos2fdrs_extended.table")));
//            for (Node go : gos) {
//                bw.write(go.node_id + "\t" + go.node_name);
//                for (double fdr : res.GOnode2FDRrunsExtend.get(go)) {
//                    bw.write("\t" + df.format(fdr));
//                }
//                bw.write("\n");
//            }
//            bw.close();
            bw = new BufferedWriter(new FileWriter(new File(this.out_dir, File.separator + "gos2fdrs_standard.table")));
            for (Node go : gos_standard.keySet()) {
                bw.write(go.node_id + "\t" + go.node_name + "\t" + df.format(gos_standard.get(go)) + "\n");
            }
            bw.close();
        } catch (IOException e) {
            throw new RuntimeException("could not write file. ", e);
        }
    }

    public void unclear_genes_BARPLOT(Collection<Gene> genes) {
        String rcommand = root + File.separator + "unclear_genes_BARPLOT.R";
        try {
            Process p = new ProcessBuilder("Rscript", rcommand, genetable_path, out_dir).inheritIO().start();
            p.waitFor();
        } catch (IOException e) {
            throw new RuntimeException("could not read/find Rscript ", e);
        } catch (InterruptedException i) {
            throw new RuntimeException("could not run subprocess ", i);
        }
    }

    public void sig_genes_VOLCANO() {
        String rcommand = root + File.separator + "sig_genes_VOLCANO.R";
        try {
            Process p = new ProcessBuilder("Rscript", rcommand, genetable_path, out_dir).inheritIO().start();
            //Process p = new ProcessBuilder("Rscript", rcommand, String.valueOf(this.FDR_cutoff), String.valueOf(this.FC_cutoff), genetable_path, out_dir).inheritIO().start();
            p.waitFor();
        } catch (IOException e) {
            throw new RuntimeException("could not read/find Rscript ", e);
        } catch (InterruptedException i) {
            throw new RuntimeException("could not run subprocess ", i);
        }
    }

    public void gene_categories_BARPLOT() {
        String rcommand = root + File.separator + "gene_categories_BARPLOT.R";
        try {
            Process p = new ProcessBuilder("Rscript", rcommand, genetable_path, out_dir).inheritIO().start();
            p.waitFor();
        } catch (IOException e) {
            throw new RuntimeException("could not read/find Rscript ", e);
        } catch (InterruptedException i) {
            throw new RuntimeException("could not run subprocess ", i);
        }
    }

    public void weighted_genes_CUMULATIVE() {
        String rcommand = root + File.separator + "genes_scored_CUMULATIVE.R";
        try {
            Process p = new ProcessBuilder("Rscript", rcommand, genetable_path, out_dir).inheritIO().start();
            p.waitFor();
        } catch (IOException e) {
            throw new RuntimeException("could not read/find Rscript ", e);
        } catch (InterruptedException i) {
            throw new RuntimeException("could not run subprocess ", i);
        }
    }

    public void genes_set_PIECHART() {
        String rcommand = root + File.separator + "genes_sets_PIECHART.R";
        try {
            Process p = new ProcessBuilder("Rscript", rcommand, genetable_path, out_dir).inheritIO().start();
            p.waitFor();
        } catch (IOException e) {
            throw new RuntimeException("could not read/find Rscript ", e);
        } catch (InterruptedException i) {
            throw new RuntimeException("could not run subprocess ", i);
        }
    }

    public void go_fdrs_mean_vs_quantile_SCATTER() {
        String rcommand = root + File.separator + "gos_fdr_mean_quantile_SCATTER.R";
        try {
            Process p = new ProcessBuilder("Rscript", rcommand, gotable_path, String.valueOf(0.95), out_dir).inheritIO().start();
            p.waitFor();
        } catch (IOException e) {
            throw new RuntimeException("could not read/find Rscript ", e);
        } catch (InterruptedException i) {
            throw new RuntimeException("could not run subprocess ", i);
        }
    }

    public void selected_gos_fdr_distrib_BOXPLOT() {
        String rcommand = root + File.separator + "selected_GOs_FDR_BOXPLOT.R";
        try {
            Process p = new ProcessBuilder("Rscript", rcommand, gotable_path, out_dir).inheritIO().start();
            p.waitFor();
        } catch (IOException e) {
            throw new RuntimeException("could not read/find Rscript ", e);
        } catch (InterruptedException i) {
            throw new RuntimeException("could not run subprocess ", i);
        }
    }

    public void selected_gos_rob_vs_extend_BOXPLOT() {
        String rcommand = root + File.separator + "selected_GOs_rob_vs_extend_BOXPLOT.R";
        try {
            Process p = new ProcessBuilder("Rscript", rcommand, gotable_path, gotable_extend_path, out_dir).inheritIO().start();
            p.waitFor();
        } catch (IOException e) {
            throw new RuntimeException("could not read/find Rscript ", e);
        } catch (InterruptedException i) {
            throw new RuntimeException("could not run subprocess ", i);
        }
    }

    public void gos_quantile_vs_mean_fdr_JITTER() {
        String rcommand = root + File.separator + "gos_quantile_vs_mean_JITTER.R";
        try {
            Process p = new ProcessBuilder("Rscript", rcommand, gotable_path, gotable_extend_path, out_dir).inheritIO().start();
            p.waitFor();
        } catch (IOException e) {
            throw new RuntimeException("could not read/find Rscript ", e);
        } catch (InterruptedException i) {
            throw new RuntimeException("could not run subprocess ", i);
        }
    }

    public void gos_standard_vs_robust_vs_extended_BARPLOT() {
        String rcommand = root + File.separator + "gos_standard_vs_robust_vs_extended_BARPLOT.R";
        try {
            Process p = new ProcessBuilder("Rscript", rcommand, gotable_path, gotable_extend_path, gotable_standard, out_dir).inheritIO().start();
            p.waitFor();
        } catch (IOException e) {
            throw new RuntimeException("could not read/find Rscript ", e);
        } catch (InterruptedException i) {
            throw new RuntimeException("could not run subprocess ", i);
        }
    }

    public void gos_standard_vs_robust_vs_extended_VENN() {
        String rcommand = root + File.separator + "gos_standard_vs_robust_vs_extended_VENN.R";
        try {
            Process p = new ProcessBuilder("Rscript", rcommand, gotable_path, gotable_extend_path, gotable_standard, out_dir).inheritIO().start();
            p.waitFor();
        } catch (IOException e) {
            throw new RuntimeException("could not read/find Rscript ", e);
        } catch (InterruptedException i) {
            throw new RuntimeException("could not run subprocess ", i);
        }
    }

    public void fdr_cutoff_finding_CUMULATIVE() {
        String rcommand = root + File.separator + "fdr_cutoff_finding_CUMULATIVE.R";
        try {
            Process p = new ProcessBuilder("Rscript", rcommand, genetable_path, out_dir).inheritIO().start();
            p.waitFor();
        } catch (IOException e) {
            throw new RuntimeException("could not read/find Rscript ", e);
        } catch (InterruptedException i) {
            throw new RuntimeException("could not run subprocess ", i);
        }
    }

    public void fc_cutoff_finding_CUMULATIVE() {
        //LATER check y-axe >= 0
        String rcommand = root + File.separator + "fc_cutoff_finding_CUMULATIVE.R";
        try {
            Process p = new ProcessBuilder("Rscript", rcommand, genetable_path, out_dir).inheritIO().start();
            p.waitFor();
        } catch (IOException e) {
            throw new RuntimeException("could not read/find Rscript ", e);
        } catch (InterruptedException i) {
            throw new RuntimeException("could not run subprocess ", i);
        }
    }

    public void flexset_extension_CURVE() {
        String rcommand = root + File.separator + "flexset_extension_CURVE.R";
        try {
            Process p = new ProcessBuilder("Rscript", rcommand, genetable_path, out_dir).inheritIO().start();
            p.waitFor();
        } catch (IOException e) {
            throw new RuntimeException("could not read/find Rscript ", e);
        } catch (InterruptedException i) {
            throw new RuntimeException("could not run subprocess ", i);
        }
    }

    public void num_sig_gos_BARPLOT() {
        String rcommand = root + File.separator + "num_sig_gos_comparison_BARPLOT.R";
        try {
            Process p = new ProcessBuilder("Rscript", rcommand, gotable_path, gotable_extend_path, gotable_standard, out_dir).inheritIO().start();
            p.waitFor();
        } catch (IOException e) {
            throw new RuntimeException("could not read/find Rscript ", e);
        } catch (InterruptedException i) {
            throw new RuntimeException("could not run subprocess ", i);
        }
    }

    public void expected_change_BARPLOT() {
        String rcommand = root + File.separator + "expected_change_distrib_BARPLOT.R";
        try {
            Process p = new ProcessBuilder("Rscript", rcommand, genetable_path, out_dir).inheritIO().start();
            p.waitFor();
        } catch (IOException e) {
            throw new RuntimeException("could not read/find Rscript ", e);
        } catch (InterruptedException i) {
            throw new RuntimeException("could not run subprocess ", i);
        }
    }

    public void de_scores_rank() {
        String rcommand = root + File.separator + "de_scores_rank_PLOT.R";
        try {
            Process p = new ProcessBuilder("Rscript", rcommand, ranktable_path, out_dir).inheritIO().start();
            p.waitFor();
        } catch (IOException e) {
            throw new RuntimeException("could not read/find Rscript ", e);
        } catch (InterruptedException i) {
            throw new RuntimeException("could not run subprocess ", i);
        }
    }
}
