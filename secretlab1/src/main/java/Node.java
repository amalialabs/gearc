import java.util.HashSet;
import java.util.Set;

public class Node {

    String node_id;
    String node_name;
    Set<Gene> genes;
//    Set<String> geneIds;
    double enrichment_score;
    Set<String> parents;
    Set<Node> parentNodes;

    public String getNode_id() {
        return node_id;
    }

    public void setNode_id(String node_id) {
        this.node_id = node_id;
    }

    public String getNode_name() {
        return node_name;
    }

    public void setNode_name(String node_name) {
        this.node_name = node_name;
    }

    public Set<Gene> getGenes() {
        if (genes == null) {
            genes = new HashSet<>();
        }
        return genes;
    }

    public void setGenes(Set<Gene> genes) {
        this.genes = genes;
    }

    public double getEnrichment_score() {
        return enrichment_score;
    }

    public void setEnrichment_score(double enrichment_score) {
        this.enrichment_score = enrichment_score;
    }

    public Set<String> getParents() {
        return parents;
    }

    public void setParents(Set<String> parents) {
        this.parents = parents;
    }

    public Set<Node> getParentNodes() {
        if (parentNodes == null) {
            parentNodes = new HashSet<>();
            for (String parent : getParents()) {
                parentNodes.add(GO.goNodes.get(parent));
            }
        }
        return parentNodes;
    }

    public void setParentNodes(Set<Node> parentNodes) {
        this.parentNodes = parentNodes;
    }

    @Override
    public String toString() {
        return "Node{" +
                "node_id='" + node_id + '\'' +
                ", node_name='" + node_name + '\'' +
                ", genes=" + genes +
                ", enrichment_score=" + enrichment_score +
                ", parents=" + parents +
                '}';
    }
}
