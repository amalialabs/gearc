import java.util.HashMap;

public class GO {

    static HashMap<String, Node> goNodes;

    public GO() {
        goNodes = new HashMap<>();
    }

    public static HashMap<String, Node> getGoNodes() {
        return goNodes;
    }

    public void setGoNodes(HashMap<String, Node> goNodes) {
        GO.goNodes = goNodes;
    }
}
