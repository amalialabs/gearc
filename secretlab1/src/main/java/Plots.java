import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;

public class Plots {

    double FDR_cutoff;
    double FC_cutoff;
    String out_dir;
    String genetable_path;
    Result res;

    public Plots(String output_dir, Collection<Gene> genes, Collection<Node> gos, double fdr, double fc, Result res) {
        this.FDR_cutoff = fdr;
        this.FC_cutoff = fc;
        this.out_dir = output_dir;
        genetable_path = out_dir + File.separator + "genes.table";
        this.res = res;
        createGeneTable(genes);
        createGOTable(gos);
    }

    public void createGeneTable(Collection<Gene> genes) {
        try {
            File dir = new File(this.out_dir); //fixme Elena -> hier vom input param nehmen
            if (!dir.exists()) {
                dir.mkdirs();
            }
            BufferedWriter bw = new BufferedWriter(new FileWriter(new File(this.out_dir, "genes.table")));
            bw.write("geneID\tFDR\tlog2FC\tgeneset\tweighted_score\tis_sig\tis_unclear\tnot_sig\n");
            for (Gene gene : genes) {
                bw.write(gene.gene_id + "\t" + gene.fdr + "\t" + gene.fc + "\t" + gene.set +
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
            BufferedWriter bw = new BufferedWriter(new FileWriter(new File(this.out_dir, "gos.table")));
            bw.write("nodeID\tnodeName\tenrichScore\tstandardFDR\n");
            for (Node go : gos) {
                bw.write(go.node_id + "\t" + go.node_name + "\t" + go.enrichment_score + "\t" + go.bhFDR + "\n");
            }
            bw.close();
            bw = new BufferedWriter(new FileWriter(new File(this.out_dir, "gos2fdrs.table")));
            for (Node go : gos) {
                bw.write(go.node_id + "\t" + go.node_name);
                for (double fdr : res.GOnode2FDRruns.get(go)) {
                    bw.write("\t" + fdr);
                }
                bw.write("\n");
            }
            bw.close();
        } catch (IOException e) {
            throw new RuntimeException("could not write file. ", e);
        }
    }

    public void unclear_genes_BARPLOT(Collection<Gene> genes) {
        String rcommand = "/secretlab1/src/main/rscripts/unclear_genes_BARPLOT.R";
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
        String rcommand = "/secretlab1/src/main/rscripts/sig_genes_VOLCANO.R";
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
        String rcommand = "/secretlab1/src/main/rscripts/gene_categories_BARPLOT.R";
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
        String rcommand = "/secretlab1/src/main/rscripts/genes_scored_CUMULATIVE.R";
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
        String rcommand = "/secretlab1/src/main/rscripts/genes_sets_PIECHART.R";
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
        String rcommand = "/secretlab1/src/main/rscripts/genes_sets_PIECHART.R";
        try {
            Process p = new ProcessBuilder("Rscript", rcommand, genetable_path, String.valueOf(0.95), out_dir).inheritIO().start();
            p.waitFor();
        } catch (IOException e) {
            throw new RuntimeException("could not read/find Rscript ", e);
        } catch (InterruptedException i) {
            throw new RuntimeException("could not run subprocess ", i);
        }
    }

}
