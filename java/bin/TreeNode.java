
import java.util.ArrayList;

public class TreeNode {

    public String name;
    public ArrayList<TreeNode> children;
    public TreeNode parent = null;;

    public TreeNode(){
        this.name = "";
        this.children = new ArrayList<>();
    }

    public TreeNode(String name){
        this.name = name;
        this.children = new ArrayList<>();
    }

    public void addChild(TreeNode child){
        this.children.add(child);
    }

    public ArrayList<TreeNode> getChildren(){
        return this.children;
    }
}
