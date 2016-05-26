import helper
import os
import helpersetting
import time
import random
import MySQLdb
def upephelper_processing(key, codons, override):
    outpath = helpersetting.UPEPHELPER_STAGING
    data_loc = helpersetting.UPEPHELPER_DATABASE
    dbuser = helpersetting.DATABASES['default']['USER']
    dbpass = str(helpersetting.DATABASES['default']['PASSWORD'])
    dbhost = helpersetting.DATABASES['default']['HOST']
    daba = helpersetting.DATABASES['default']['DB']
    #local_version = []
    timeid = time.time()
    unid = '%012x%016x' % (int(timeid * 1000), random.randint(0, 0xFFFFFFFFFFFFFFFF))
    
    #if Refseqdb_blast_db_build_log.objects.all().exists() is False:
        #local_version.append(0)

        #print('No local database version information available.')
    #else:
        #r = Refseqdb_blast_db_build_log.objects.order_by('-input_date')[0]
        #local_version.append(r.database_version)

        #print('Latest local database version: %i' % r.database_version)
    
    remote = helper.get_NCBI_RefSeq_release()
    
    
    query= ','.join(codons)
    override_condition = "True"
    dbv = str(remote)
    helper.upep_mysql_database(unid, dbuser, dbpass, dbhost, daba, key, query, override_condition, remote)
    local_version = 0
    lv = """select * from updater_log where refseq_database = %s and success_log = 1 order by unix_timestamp(time_id_start) desc;"""
    dbcon = MySQLdb.connect(user=dbuser, passwd=dbpass, host=dbhost, db=daba)
    cursor = dbcon.cursor()
    cursor.execute(lv, (key,))
    if not cursor.rowcount:
        local_version = 0
    else:
        local = cursor.fetchone()
        local_version = local[4]
    cursor.close()
    dbcon.close()
    fn = 0
    if remote > local_version or override_condition == "True":
        home = os.getcwd()
        dbs = ['RefSeq-complete',
               'RefSeq-fungi',
               'RefSeq-invertebrate',
               'RefSeq-plant',
               'RefSeq-vertebrate_mammalian',
               'RefSeq-vertebrate_other']
        if key:
            print("Working with database " + key)

            if key in dbs:
                
                wd = helper.setup(outpath, home, key)
                fn = helper.download_db(key)
                
                compacted = helper.compact_RefSeq(wd, dbv)
                
                print("Compiling ACC and GI database for " + key)
                helper.compile_RefSeq(compacted, dbv, fn, dbuser, dbpass, dbhost, daba, timeid)                
                               
                print("Recorded log for building ACC and GI database of " + key)
                os.chdir(home)
            else:

                print("Not a defined db")
                sys.exit(1)
        #else:
            #for db in dbs:

                #print("Working with database" + db)
                #wd = helper.setup(outpath, home, db)
                #helper.download_db(db)
                #compacted = helper.compact_RefSeq(wd, dbv)
                #helper.compile_RefSeq(compacted, dbv, dbuser, dbpass, dbhost)
                
                #os.chdir(home)

        os.chdir(outpath)
        try:
            os.mkdir("../tmp/RefSeqdb"+dbv)
        except OSError:
            print("RefSeqdb %s directory exists \n Overriding directory" % dbv)
            pass
        os.system("mv "+outpath+"RefSeq* ../tmp/RefSeqdb"+dbv)
        
        for starting_codon in codons:
            _, proc_list = helper.uPEP_finder(codon=starting_codon, db_version=dbv, outpath=outpath, fn=fn)
            helper.build_blast_db(proc_list)
        # Here is where we should go into maintainence mode
        helper.finalise_update()
        dbcon = MySQLdb.connect(user=dbuser, passwd=dbpass, host=dbhost, db=daba)
        cursor = dbcon.cursor()
    
        update_log = """UPDATE updater_log SET time_id_finish = current_timestamp, success_log = 1 WHERE unique_id = %s"""
        cursor.execute(update_log, (unid,))
        dbcon.commit()
        cursor.close()
        dbcon.close()
        #for starting_codon in codons:
            #rdb_blast_log = Refseqdb_blast_db_build_log(input_date=timezone.now(), database_name = key, database_version = remote, codon = starting_codon)
            #rdb_blast_log.save()

    else:
        print("No updrade required")

