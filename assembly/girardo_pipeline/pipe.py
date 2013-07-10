import psycopg2
import os
import sys

errr = "\nStill processing previous job\n Aborting.\n..."
fromPath = "/var/www/output/Home/"
superComp = "girardo@128.227.70.246"
superPath = "/srv/data2/pipeline/"
assemblePath = superPath + "assembly/"
annotatePath = superPath + "annotation/"
predictionPath = superPath + "prediction/"

sqlQuery = "SELECT experiment_id,\"fastqLink\",project,assembly FROM rundb_experiment,rundb_results WHERE rundb_experiment.id=experiment_id AND status='Completed';"
dbConn = "dbname=iondb user=postgres"
def sql_update(asmbly,expid) = "UPDATE rundb_results SET assembly = " + asmbly + " WHERE experiment_id = " + expid + ";"
def running_update(onoff) = "UPDATE pipeline SET running = " + onoff + ";"


conn = psycopg2.connect(dbConn) #
cur = conn.cursor() #Opens database connection
cur.execute("SELECT * from pipeline;")
if cur.fetchone()==(True,) and not "-f" in sys.argv: print "Pipepline is still running" #If the pipeline is not running...
else:
    cur.execute(sqlQuery) #
    filter_fun = lambda x: not x[3]=='ignore'$ #
    files = filter(filter_fun,cur.fetchall()) #Retrieves new projects
    cur.execute(runningUpdate(True)) #Indicates pipeline is running
    conn.commit
    for expe in files: #For each project...
        os.popen("ssh " + superComp + " mkdir " + superPath + expe[2] #Creates the project directory
        os.popen("scp /var/www/" + expe[1] + " " + superComp + ":" + filePath + expe[2] + "_in.solexa.fastq") #send the fastq file
        os.popen("ssh " + superComp + " " + assemblePath + "assemble.py " + expe[2] #run the assembly script
        os.popen("ssh " + superComp + " " + annotatePath + "annoate.py " + expe[2] #run the annotation script
        os.popen("ssh " + superComp + " " + predictionPath + "predicts.py " + expe[2] #run the prediction script    
        os.popen("ssh " + superComp + " " + superPath + "uploading.py -s" + expe[2] #run the database upload script    
        cur.execute(sql_update("ignore",expe[0]) #Indicate project is done
        cur.execute("ssh " + superComp + " rm -r " + superPath + expe[2] #Remove the project directory
    cur.execute(runningUpdate(False)) #Pipeline is done
    conn.commit
