import java.io.*;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Stream;


public class EnsemblGOMapping {

    public static void main(String[] args) {
        File geneIDsFile = new File("/home/birinci/Downloads/Homo_sapiens.GRCh38.106.gtf.geneIDs");
        EnsemblGOMapping egm = new EnsemblGOMapping(new File("/home/birinci/Downloads/"), 106, getGeneIDs(geneIDsFile));
//        getGO("ENSG00000170152").stream().forEach(_go -> System.out.println(_go));
    }

    public EnsemblGOMapping(File outputDir, int ensemblVersion, List<String> geneIDs) {
        File outputFile = new File(outputDir,"goa_human_ensembl_release_" + ensemblVersion + ".tsv");

        try (PrintWriter pw = new PrintWriter(outputFile)) {
            pw.println("ensembl_id\thgnc\tgos");
            for (String gene : geneIDs) {
                pw.print(gene + "\t");
                String geneName = EnsemblRestClient.getGeneName(gene);
                pw.print(geneName + "\t");
                Collection<String> gos = getGO(gene);
                if (gos.size() > 0) {
                    StringBuilder sb = new StringBuilder();
                    for (String elem : getGO(gene)) {
                        sb.append("|").append(elem);
                    }
                    if (sb.length() > 0) {
                        pw.println(sb.toString().substring(1));
                    } else {
                        System.out.println("Error in:\t" + gene + "\t" + geneName + "\t" + sb.toString());
                    }
                } else
                    pw.println();
            }
        } catch (IOException e) {
            throw new RuntimeException("Can not write to outputfile\t" + outputFile.getAbsolutePath(), e);
        }
    }

    private static List<String> getGeneIDs(File geneListFile) {
        List<String> list = new ArrayList<>();
        try (Stream<String> stream = Files.lines(Paths.get(geneListFile.getAbsolutePath()))) {
            stream.forEach(list::add);
        } catch (IOException e) {
            throw new RuntimeException("Error reading expression file: ", e);
        }
        return list;
    }

    /**
     * TODO evaluate if we want to filter by specific evidence level
     * Guide to evidence codes from GO:
     *     EXP - Inferred from Experiment
     *     IC - Inferred by Curator
     *     IDA- Inferred from Direct Assay
     *     IEA - Inferred from Electronic Annotation
     *     IEP - Inferred from Expression Pattern
     *     IGC - Inferred from Genomic Context
     *     IGI - Inferred from Genetic Interaction
     *     IMP - Inferred from Mutant Phenotype
     *     IPI - Inferred from Physical Interaction
     *     ISA - Inferred from Sequence Alignment
     *     ISM - Inferred from Sequence Model
     *     ISO - Inferred from Sequence Orthology
     *     ISS - Inferred from Sequence or Structural Similarity
     *     NAS - Non-traceable Author Statement
     *     ND - No biological Data available
     *     RCA - inferred from Reviewed Computational Analysis
     *     TAS - Traceable Author Statement
     *     NR - Not Recorded
     * @param geneID
     * @throws IOException
     */
    public static Collection<String> getGO(String geneID) {
        String endpoint = "/xrefs/id/" + geneID + "?external_db=GO;all_levels=1";
        String output = "";
        try {
            output = EnsemblRestClient.getContent(endpoint);
        } catch (InterruptedException | IOException e) {
            e.printStackTrace();
        }

        String[] elems = output.split("\\\"");
        Set<String> gos = new HashSet<>();
        for (String elem : elems) {
            if (elem.startsWith("GO:")) gos.add(elem);
        }
        return gos;
    }
}
