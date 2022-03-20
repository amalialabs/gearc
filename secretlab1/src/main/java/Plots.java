import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;

public class Plots {

    double FDR_cutoff;
    double FC_cutoff;
    String out_dir;
    String rscript_folder = "/secretlab1/src/main/rscripts/";  //locally in docker image

    public Plots(String output_dir) {
        this.FDR_cutoff = 0.05;
        this.FC_cutoff = 1.0;
        this.out_dir = output_dir;
    }

    public String list2vector(ArrayList<Object> list) {
        StringBuilder sb = new StringBuilder();
        sb.append("c(");
        for (Object elem : list) {
            sb.append(elem.toString()).append(",");
        }
        sb = sb.deleteCharAt(sb.length()-1);
        sb.append(")");
        return sb.toString();
    }

    public void unclear_genes_BARPLOT(Collection<Gene> genes) {
        ArrayList<Object> gene_ids = new ArrayList<>();
        ArrayList<Object> gene_unclear = new ArrayList<>();
        genes.stream().forEach(_gene -> {
            if (_gene.fdr <= FDR_cutoff && Math.abs(_gene.fc) >= FC_cutoff) {
                gene_ids.add(_gene.gene_id);
                gene_unclear.add("clear");
            } else if (_gene.fdr > FDR_cutoff && Math.abs(_gene.fc) < FC_cutoff) {
                gene_ids.add(_gene.gene_id);
                gene_unclear.add("clear");
            } else {
                gene_ids.add(_gene.gene_id);
                gene_unclear.add("unclear");
            }
        });
        String g = list2vector(gene_ids);
        String grp = list2vector(gene_unclear);
        try {
            Process p = new ProcessBuilder("Rscript ", rscript_folder+"unclear_genes_BARPLOT.R", g, grp, out_dir).inheritIO().start();
            p.waitFor();
            //Runtime.getRuntime().exec(rscript_exe + " " + rscript_folder+"unclear_genes.R " + g + " " + grp + " " + out_dir);
        } catch (IOException e) {
            throw new RuntimeException("could not read/find Rscript ", e);
        } catch (InterruptedException i) {
            throw new RuntimeException("could not run subprocess ", i);
        }
    }

    public void sig_genes_VOLCANO(Collection<Gene> genes) {
        ArrayList<Object> gene_ids = new ArrayList<>();
        ArrayList<Object> gene_fc = new ArrayList<>();
        ArrayList<Object> gene_fdr = new ArrayList<>();
        genes.stream().forEach(_gene -> {
            gene_ids.add(_gene.gene_id);
            gene_fc.add(_gene.fc);
            gene_fdr.add(_gene.fdr);
        });
        String g = list2vector(gene_ids);
        String fc = list2vector(gene_fc);
        String fdr = list2vector(gene_fdr);
        try {
            Process p = new ProcessBuilder("Rscript ", rscript_folder+"significant_genes_VOLCANO.R", g, fc, fdr, out_dir).inheritIO().start();
            p.waitFor();
        } catch (IOException e) {
            throw new RuntimeException("could not read/find Rscript ", e);
        } catch (InterruptedException i) {
            throw new RuntimeException("could not run subprocess ", i);
        }
    }

    public void gene_categories_BARPLOT() {

    }

    public void weighted_genes_CUMULATIVE() {

    }

}
