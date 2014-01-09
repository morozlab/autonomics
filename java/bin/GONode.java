import java.util.ArrayList;
import java.util.HashMap;

public class GONode extends TreeNode{

    public String GOID;
    public int numSeqs;
    public int id;
    public int abundance;
    public HashMap<String, Integer> sbIDs;

    public GONode(String name, String GOID){
        super(name);
        this.GOID = GOID;
        this.numSeqs = 0;
        this.children = new ArrayList<TreeNode>();
        this.sbIDs = new HashMap<String, Integer>();
        this.abundance = 0;
    }

    public GONode(String name, String GOID, int dbid){
        super(name);
        this.GOID = GOID;
        this.id = dbid;
        this.numSeqs = 0;
        this.children = new ArrayList<TreeNode>();
        this.sbIDs = new HashMap<String, Integer>();
        this.abundance = 0;
    }
}
