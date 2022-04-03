import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;

public class Plots {

    double FDR_cutoff;
    double FC_cutoff;
    String out_dir;
//    String gene_ids;
//    String gene_groups;
//    String gene_fcs;
//    String gene_fdrs;
    String genetable_path;

    public Plots(String output_dir, Collection<Gene> genes, Collection<Node> gos, double fdr, double fc) {
        this.FDR_cutoff = fdr;
        this.FC_cutoff = fc;
        this.out_dir = output_dir;
        genetable_path = out_dir + File.separator + "genes.table";
        //init_r_vectors(genes);
        createGeneTable(genes);
        createGOTable(gos);
    }

//    public String list2vector(ArrayList<Object> list) {
//        StringBuilder sb = new StringBuilder();
//        sb.append("c(");
//        for (Object elem : list) {
//            sb.append(elem.toString()).append(",");
//        }
//        sb = sb.deleteCharAt(sb.length()-1);
//        sb.append(")");
//        return sb.toString();
//    }

    public void createGeneTable(Collection<Gene> genes) {
        try {
            File dir = new File(this.out_dir); //fixme Elena -> hier vom input param nehmen
            if (!dir.exists()) {
                dir.mkdirs();
            }
            BufferedWriter bw = new BufferedWriter(new FileWriter(new File(this.out_dir, "genes.table")));
            bw.write("geneID\tFDR\tlog2FC\tgeneset\tweighted_score\tis_sig\n");
            for (Gene gene : genes) {
                bw.write(gene.gene_id + "\t" + gene.fdr + "\t" + gene.fc + "\t" + gene.set +
                        "\t" + gene.weighted_score + "\t" + gene.is_significant + "\n");
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
            bw.write("nodeID\tnodeName\tenrichScore\tFDR\n");
            for (Node go : gos) {
                bw.write(go.node_id + "\t" + go.node_name + "\t" + go.enrichment_score + "\t" + go.bhFDR + "\n");
            }
            bw.close();
        } catch (IOException e) {
            throw new RuntimeException("could not write file. ", e);
        }
    }

//    public void init_r_vectors(Collection<Gene> genes) {
//        ArrayList<Object> gene_ids = new ArrayList<>();
//        ArrayList<Object> gene_groups = new ArrayList<>();
//        ArrayList<Object> gene_fc = new ArrayList<>();
//        ArrayList<Object> gene_fdr = new ArrayList<>();
//        ArrayList<Object> gene_sets = new ArrayList<>();
//        genes.stream().forEach(_gene -> {
//            if (_gene.fdr <= FDR_cutoff && Math.abs(_gene.fc) >= FC_cutoff) {
//                gene_ids.add(_gene.gene_id);
//                gene_groups.add("clear");
//            } else if (_gene.fdr > FDR_cutoff && Math.abs(_gene.fc) < FC_cutoff) {
//                gene_ids.add(_gene.gene_id);
//                gene_groups.add("clear");
//            } else {
//                gene_ids.add(_gene.gene_id);
//                gene_groups.add("unclear");
//            }
//            gene_fc.add(_gene.fc);
//            gene_fdr.add(_gene.fdr);
//            gene_sets.add(_gene.set);
//        });
//        this.gene_fcs = list2vector(gene_fc);
//        this.gene_fdrs = list2vector(gene_fdr);
//        this.gene_ids = list2vector(gene_ids);
//        this.gene_groups = list2vector(gene_groups);
//    }

    public void unclear_genes_BARPLOT(Collection<Gene> genes) {
        String rcommand = "/secretlab1/src/main/rscripts/unclear_genes_BARPLOT.R";
//        try {
//            Process p = new ProcessBuilder("Rscript", rcommand, this.gene_ids, this.gene_groups, out_dir).inheritIO().start();
//            p.waitFor();
//            //Runtime.getRuntime().exec(rscript_exe + " " + rscript_folder+"unclear_genes.R " + g + " " + grp + " " + out_dir);
//        } catch (IOException e) {
//            throw new RuntimeException("could not read/find Rscript ", e);
//        } catch (InterruptedException i) {
//            throw new RuntimeException("could not run subprocess ", i);
//        }
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

    public void selected_GOs_fdrs_BOXPLOT() {

    }

//    public static void main(String[] args) {
//        List<Gene> list = new ArrayList<>();
//        Plots plots = new Plots("/home/birinci/GOEnrichment/", list, 0, 0);
//        plots.sig_genes_VOLCANO(null);
//    }

}
