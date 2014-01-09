import bioinformatics.files.FileUtils;
import bioinformatics.files.tab.DelimitedFile;
import bioinformatics.files.tab.DelimitedRecord;

import java.io.File;
import java.io.IOException;
import java.sql.*;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.HashMap;
import java.util.Set;


public class GOTree {

    private final DelimitedFile goFile;
    protected GONode tree;
    private String projectName;
    private String goPath;
    private String quantPath;
    private String uname;
    private String db_name;
    private String db_host;
    private String passwd;
    private Connection con = null;
    private Statement stmt = null;
    private HashMap<String, ArrayList<GONode>> relationships = new HashMap<>();
    private HashMap<String, HashMap<String, Integer>> goCatDB = new HashMap<>();
    private HashMap<String, ArrayList<String>> goAnnotations;
    private HashMap<Integer, Integer> seen = new HashMap<>();
    private String url = "jdbc:mysql://";
    /*    private String url = "jdbc:mysql://localhost:3306/zero_click"; */
    int iterations = 0;

    public GOTree(String projectName, String goPath, String quantPath, String uname, String passwd, String db_name,  String db_host) throws Exception{
    	this.projectName = projectName;
    	this.goPath = goPath;
    	this.quantPath = quantPath;
        this.uname = uname;
        this.db_name = db_name;
        this.db_host = db_host;
        this.passwd = passwd;
        url = url + this.db_host + ":3306/" + this.db_name;
        tree = new GONode("root", "");
        Class.forName("com.mysql.jdbc.Driver");
        try{
            this.con = DriverManager.getConnection(url, this.uname, this.passwd);
            this.stmt = this.con.createStatement();
        }
        catch(Exception e){
            System.exit(-1);
        }       
       // String[] fields = {"project_name:str", "go_id:str", "seq_id:str", "sp_acc:str", "evalue:float"};
        String[] fields = {"go_id:2:str", "seq_id:0:str", "sp_acc:1:str", "evalue:7:float"};
        int[] keyCols = {0};
        //generate a hash of annotations
        this.goFile = new DelimitedFile(goPath, fields, keyCols, "\t");
        this.goFile.hashify();
        
    }

    public GONode construct(){
    	
        //first, populate a hash with the children for each GO term
        try{
            ResultSet rs = this.stmt.executeQuery("SELECT term1_id, term2_id, acc, name FROM term2term, term"
                    + " WHERE term2term.term2_id = term.id");
            while(rs.next()){
                GONode child = new GONode(rs.getString("name"), rs.getString("acc"), rs.getInt("term2_id"));
                if(this.relationships.containsKey(rs.getString("term1_id"))){
                    this.relationships.get(rs.getString("term1_id")).add(child);
                }
                else{
                    ArrayList<GONode> tmp = new ArrayList<GONode>();
                    tmp.add(child);
                    this.relationships.put(rs.getString("term1_id"), tmp);
                }
            }
        }
        catch(Exception e){
            e.printStackTrace();
        }

        //attach the children nodes to parent nodes
        this.tree.addChild(new GONode("Biological Process", "GO:0008150", 6490));
        this.tree.addChild(new GONode("Molecular Function", "GO:0003674", 2668));
        this.tree.addChild(new GONode("Cellular Component", "GO0005575", 4396));
        Iterator<TreeNode> it = this.tree.children.iterator();
        while(it.hasNext()){
            attachChildren((GONode) it.next());
        }
        /*tree.addChild(generateSubTree(6490, "Biological Process", "GO:0008150"));
        tree.addChild(generateSubTree(2668, "Molecular Function", "GO:0003674"));
        tree.addChild(generateSubTree(4396, "Cellular Component", "GO:0005575"));
        */
        //System.out.println("Iterations: " + iterations);
        relationships = new HashMap<String, ArrayList<GONode>>();
        return tree;
    }

    private void attachChildren(GONode n) {
        //System.out.println(n.name);
        //check to see if this node has any children
        if(relationships.containsKey(String.valueOf(n.id))){
            //add all the children to the node
            Iterator<GONode> it = relationships.get(String.valueOf(n.id)).iterator();
            //System.out.println("Parent: " + n.name + " " + n.GOID);
            while(it.hasNext()){
                //iterations++;
                GONode tmp = it.next();
                //duplicate the node from the hash
                GONode child = new GONode(tmp.name, tmp.GOID, tmp.id);
                //System.out.println("Child: " + child.name + " " + child.GOID);
                n.addChild(child);
                //System.out.println("Down");
                attachChildren(child);
                //System.out.println("up");
            }
        }
    }

    public void count() {
        //populate a hash of hashes for each go ID, storing sbIDs seen for that category
        try{
        	DelimitedFile quantFile = null;
        	try{
        		quantFile = new DelimitedFile(this.quantPath, "\t");
        		if(quantFile.exists()){
        			quantFile.nameCol(0, "seq_id");
        			quantFile.nameCol(1, "abundance");
        			quantFile.hashify();
        		}
        	}
        	catch(IOException e){
        		quantFile = null;
        	}
        	Iterator<String> keys = this.goFile.data.keySet().iterator();
        	while(keys.hasNext()){
        		String key = keys.next();
        		Iterator<DelimitedRecord> records = this.goFile.data.get(key).iterator();
        		while(records.hasNext()){
        			DelimitedRecord r = records.next();
        			String seqID = r.get("seq_id").trim().split(" ")[0];
        			String goID = r.get("go_id");
        			goID = goID.split(":")[1];
        			//look up the abundance information for this sequence
        			int abundance = 1;
        			if(quantFile != null && quantFile.data.containsKey(seqID)){
        				abundance = quantFile.data.get(seqID).get(0).getInt("abundance");
        				
        			}
        			if(goCatDB.containsKey(goID)){
        				goCatDB.get(goID).put(seqID, abundance);
        			}
        			else{
        				HashMap<String, Integer> tmp = new HashMap<String, Integer>();
        				tmp.put(seqID, abundance);
        				goCatDB.put(goID, tmp);
        			}
        		}
        	}

            /*Set<String> keys = (Set<String>) goCatDB.keySet();
            Iterator it = keys.iterator();
            while(it.hasNext()){
                String key = (String) it.next();
                //System.out.println(key + " has " + goCatDB.get(key).size() + " annotations");
            }*/
        }
        catch(Exception e){
            e.printStackTrace();
        }
        //traverse the tree, keeping a list of unique sequences seen at each level
        doCount((GONode) tree.children.get(0));
        doCount((GONode) tree.children.get(1));
        doCount((GONode) tree.children.get(2));
        /*Iterator it = tree.children.iterator();
        while(it.hasNext()){
            doCount((GONode) it.next());
        }*/
        //doCount((GONode) tree.children.get(0), 1);
    }

    private void doCount(GONode t){
        String[] IDArray = t.GOID.split(":");
        String GO = "";
        if(IDArray.length != 1){
            GO = IDArray[1];
        }
        else{
            GO = t.GOID.replace("GO", "");
        	t.GOID = t.GOID.replace("GO", "GO:");
        }
        //first, make sure this node's hash has all the IDs annotated to this category
        if(goCatDB.containsKey(GO)){
            mergeHashMaps(t.sbIDs, goCatDB.get(GO));
        }

        //then, go down to children's level and do the same
        Iterator it = t.children.iterator();
        while(it.hasNext()){
            GONode child = (GONode) it.next();
            doCount(child);
            //merge child's hash map with this node's hash map
            mergeHashMaps(t.sbIDs, child.sbIDs);
            //merge the child's list of unique go/sbid pairs with this node's list
            //mergeHashMaps(t.sbGOMap, child.sbGOMap);
        }
        //store the number of sequences for this ID
        t.numSeqs = t.sbIDs.size();
        //calculate the abundance for this category
        it = t.sbIDs.keySet().iterator();
        while(it.hasNext()){
            t.abundance += t.sbIDs.get((String) it.next()).intValue();
        }
    }

    public void countUniqueAnnotations() {
        doCountUnique((GONode) tree.children.get(0));
        doCountUnique((GONode) tree.children.get(1));
        doCountUnique((GONode) tree.children.get(2));
    }

    private void doCountUnique(GONode t){
        String[] IDArray = t.GOID.split(":");
        String GO = "";
        if(IDArray.length != 1){
            GO = IDArray[1];
        }
        else{
            GO = t.GOID;
        }
        //first, make sure this node's hash has all the IDs annotated to this category
        if(goCatDB.containsKey(GO)){
            //keep track of unique go/sbid pairs for this go term
            populateGOSBMap(t.sbIDs, goCatDB.get(GO), GO);
        }
        //then, go down to children's level and do the same
        Iterator it = t.children.iterator();
        while(it.hasNext()){
            GONode child = (GONode) it.next();
            
            doCountUnique(child);
            //merge child's hash map with this node's hash map of unique annotations
            mergeHashMaps(t.sbIDs, child.sbIDs);
            //clear the contents of the child hash
            child.sbIDs.clear();
        }
        //remove the children
        t.children.clear();
        //print the results of the counting to be stored in the DB
        if(t.numSeqs > 0){
            if(!seen.containsKey(t.id)){
                System.out.println(projectName + "\t" + t.id + "\t" + t.GOID + "\t" + t.name + "\t" + t.numSeqs + "\t" + t.sbIDs.size() + "\t" + t.abundance);
                seen.put(t.id, 1);
            }
        }
    }

    public void clearNodes(){
        doClear((GONode) tree.children.get(0));
        doClear((GONode) tree.children.get(1));
        doClear((GONode) tree.children.get(2));
    }

    private void doClear(GONode t) {
        //clear this node's hash table
        t.sbIDs = new HashMap<String, Integer>();
        //clear all the children
        Iterator it = t.children.iterator();
        while(it.hasNext()){
            doClear((GONode) it.next());
        }
    }

    private void mergeHashMaps(HashMap<String, Integer> hash1, HashMap<String, Integer> hash2) {
        //method takes two hashes, and adds every element from the second to the first
        Set keys = hash2.keySet();
        Iterator it = keys.iterator();
        while(it.hasNext()){
            String key = (String) it.next();
            hash1.put(key, hash2.get(key));
        }

    }
    

    public void store(){
        
    }

    private void populateGOSBMap(HashMap<String, Integer> sbGOMap, HashMap<String, Integer> sbIDList, String term) {
        //concatenate term to each sbID in the sbIDList and then make sure they are in the sbGOMap object
        Set keys = sbIDList.keySet();
        Iterator it = keys.iterator();
        while(it.hasNext()){
            String key = (String) it.next();
            sbGOMap.put(term + "_" + key, 1);
        }
    }


}
