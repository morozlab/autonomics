
public class GOAnnotator {

    public static void main(String[] args) throws Exception{

        String projectName = args[0];
        String goFile = args[1];
        String uname = args[2];
        String passwd = args[3];
        String db_name = args[4];
        String db_host = args[5];
        String quantFile = "";
        if(args.length > 6){
        	quantFile = args[6];
        }

        GOTree t = new GOTree(projectName, goFile, quantFile, uname, passwd, db_name, db_host);
        t.construct();
        t.count();
        t.clearNodes();
        t.countUniqueAnnotations();
    
    }
}
