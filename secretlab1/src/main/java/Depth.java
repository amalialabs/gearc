import org.apache.commons.lang3.StringUtils;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Stream;

/**
 * @author hadziahmetovic on 2022-12-11
 */
public class Depth {
    public static void main(String[] args) {
        //TODO remove paths
        File expression = new File("/mnt/raidinput2/tmp/hadziahmetovic/GSE/results_reduced/rep2_mock8_wt8/diff_exp_outs/DESeq_hisat.reformatted");
        File mappingEnsembl = new File("/mnt/raidinput2/tmp/hadziahmetovic/GSE130342/gse/secretlab/data/goa_human_ensembl.tsv");
        File oboFile = new File("/mnt/raidinput2/tmp/hadziahmetovic/GSE130342/gse/secretlab/data/go.obo");
        String root = "biological_process";

        Reader r = new Reader(expression, mappingEnsembl, oboFile, root, 0.01, 1.0);

        List<Node> nodes = new ArrayList<>();
        Set<String> search = new HashSet<>();

        try (Stream<String> stream = Files.lines(Paths.get("/mnt/raidinput2/tmp/hadziahmetovic/GSE/gse_result_reduced_high/rep2_mock12_wt12_DESeq_hisat/standard.only"))) {
            stream.forEach(_line -> {
                search.add(_line);
                nodes.add(GO.getGoNodes().get(_line));
            });
        } catch (IOException e) {
            throw new RuntimeException("", e);
        }

        try (PrintWriter pw = new PrintWriter(new File("/mnt/raidinput2/tmp/hadziahmetovic/GSE/gse_result_reduced_high/rep2_mock12_wt12_DESeq_hisat/standard.only.eval"))) {
            for (Node n : nodes) {
                pw.println(n.node_id + "\t" + r.depthCheck(n, root) + "\t" + n.node_name + "\t" + StringUtils.join(r.iterativeParentSearch(n, search), ", "));
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }

        System.out.println("\n");

        nodes.clear();
        search.clear();

        try (Stream<String> stream = Files.lines(Paths.get("/mnt/raidinput2/tmp/hadziahmetovic/GSE/gse_result_reduced_high/rep2_mock12_wt12_DESeq_hisat/robust"))) {
            stream.forEach(_line -> {
                search.add(_line);
                nodes.add(GO.getGoNodes().get(_line));
            });
        } catch (IOException e) {
            throw new RuntimeException("", e);
        }

        try (PrintWriter pw = new PrintWriter(new File("/mnt/raidinput2/tmp/hadziahmetovic/GSE/gse_result_reduced_high/rep2_mock12_wt12_DESeq_hisat/robust.eval"))) {
            for (Node n : nodes) {
                pw.println(n.node_id + "\t" + r.depthCheck(n, root) + "\t" + n.node_name + "\t" + StringUtils.join(r.iterativeParentSearch(n, search), ","));
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }


        nodes.clear();
        search.clear();

        try (Stream<String> stream = Files.lines(Paths.get("/mnt/raidinput2/tmp/hadziahmetovic/GSE/gse_result_reduced_high/rep2_mock12_wt12_DESeq_hisat/standard"))) {
            stream.forEach(_line -> {
                search.add(_line);
                nodes.add(GO.getGoNodes().get(_line));
            });
        } catch (IOException e) {
            throw new RuntimeException("", e);
        }

        try (PrintWriter pw = new PrintWriter(new File("/mnt/raidinput2/tmp/hadziahmetovic/GSE/gse_result_reduced_high/rep2_mock12_wt12_DESeq_hisat/standard.eval"))) {
            for (Node n : nodes) {
                pw.println(n.node_id + "\t" + r.depthCheck(n, root) + "\t" + n.node_name + "\t" + StringUtils.join(r.iterativeParentSearch(n, search), ","));
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }
}
