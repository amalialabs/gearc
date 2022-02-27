import java.util.HashMap;

/**
 * @author hadziahmetovic on 27.02.22
 */
public class GO {


    //fixme add test to see if goParents == goNodes -> check if go file complete


    HashMap<String, Node> goNodes;

    public GO() {
        this.goNodes = new HashMap<>();
    }

    public HashMap<String, Node> getGoNodes() {
        return goNodes;
    }

    public void setGoNodes(HashMap<String, Node> goNodes) {
        this.goNodes = goNodes;
    }
}
