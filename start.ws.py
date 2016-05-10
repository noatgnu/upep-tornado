import tornado.ioloop
import tornado.web
import os
from tornado.httpclient import AsyncHTTPClient
from tornado import gen, websocket, httpserver
import settingmain
from concurrent.futures import ThreadPoolExecutor
from redis import Redis
from rq import Queue
import json
import time
from upep_helper import helpersetting, helpertasks
from rq.job import Job
from tornado.escape import json_encode, json_decode
import csv
from upep import upepsetting
import random
from upep import conserveduPEP, images, tblastx, accessions, whm
import MySQLdb
import subprocess
import sys


class wsUI(tornado.web.RequestHandler):
    def get(self):
        self.render("websocket.html")

class WebSocketHandler(websocket.WebSocketHandler):
    def open(self):
        pass

    @gen.coroutine
    def on_message(self, parameters):
        
        input_p = json.loads(parameters)
        print(input_p)
        q = Queue(connection=Redis())
        pool = ThreadPoolExecutor(max_workers=settingmain.MAX_WORKERS)
        
        if input_p['task'] == 1:
            yield pool.submit(self.uPEPdataprocessor, input_p['unid'], input_p['timeid'], input_p['codon'], input_p['dbv'], input_p['seqquery'], input_p['db'], input_p['min'], input_p['max'], input_p['grace'], input_p['hms'], input_p['kaks'], input_p['window'], input_p['gradi'], input_p['rhm'])
        if input_p['task'] == 0:
            self.write_message('Complete.')
            
            #weboutput(self, input_p['unid'], input_p['filter'], input_p['strict'])
        #if input_p['task'] == 0:
            #weboutput(self, input_p['unid'], input_p['filter'], input_p['strict'])
            
        
        #weboutput(unid, timeid, codon, dbv, seqquery,
                              #db_name, minm,
                              #maxm, gracem, hms, kaks,
                              #window, gradi, refheatmaps)
         
    def on_close(self):
        pass
    
    def uPEPdataprocessor(self, unique_file_id, tid, cod, db_version, seqqu, database_n, minm, maxm, gracem, heatms, KaKs_t, window_t, grad, refhm):
        
        upep_loc = upepsetting.APP_ROOT
        apps_loc = upepsetting.APPS
        lagan_loc = upepsetting.LAGAN
        
        upeplib_loc = upepsetting.UPEPLIB
        data_loc = upepsetting.DATA
        database = upepsetting.UPEPHELPER_DATABASE

        dbuser = upepsetting.DATABASES["default"]["USER"]
        dbpass = str(upepsetting.DATABASES["default"]["PASSWORD"])
        dbhost = upepsetting.DATABASES["default"]["HOST"]
        daba = upepsetting.DATABASES["default"]["DB"]

        tempfilename = os.path.join(data_loc, unique_file_id)

        codon = cod.replace("U", "t").replace("A", "a").replace("G", "g").replace("C", "c")
        
        mammaldb = os.path.join(database, db_version + upepsetting.MAMMAL + codon + ".db")
        non_mammalvertdb = os.path.join(database, db_version + upepsetting.NON_MAMMALIAN_VERTEBRATES + codon + ".db")
        invertdb = os.path.join(database, db_version + upepsetting.INVERTEBRATES + codon + ".db")
        plantdb = os.path.join(database, db_version + upepsetting.PLANTS + codon + ".db")
        fungidb = os.path.join(database, db_version + upepsetting.FUNGI + codon + ".db")
        completedb = os.path.join(database + db_version + upepsetting.COMPLETE + codon + ".db")
        dbcon = MySQLdb.connect(user=dbuser, passwd=dbpass, host=dbhost, db=daba)
        cursor = dbcon.cursor()
        codon_aminoa = {"AUG": "M", "AUA": "I", "AUC": "I", "ACG": "T", "AUU": "I", "AAG": "K", "AGG": "R", "CUG": "L",
                        "UUG": "L", "GUG": "V"}
        
        sql_insert_query_info = """INSERT INTO query_info (time_id, unique_id , seqquery, database_name, db_version, protein_name, tempfilename, initial_flag, trivial_flag, blastmacthes_flag, CDSKaKs_flag, uPEPKaKs_flag, heatmap_flag, window, Refheatmap, starting_aa) VALUES (from_unixtime("%f"), "%s", 0, 0, 0, 0, "%s", 0, 0, 0, 0, 0, 0, 0, 0, "%s");"""
        
        #print(sql_insert_query_info % (tid, unique_file_id, tempfilename, codon_aminoa[cod]))
        cursor.execute(sql_insert_query_info % (tid, unique_file_id, tempfilename, codon_aminoa[cod]))
        dbcon.commit()
        cursor.close()
        dbcon.close()
        CDSfeature = []
        calcKaKs = True
        dictionaries = dict([("Human", [[mammaldb], "Homo sapiens"]),
                             ("Mouse", [[mammaldb], "Mus musculus"]),
                             ("Mammals", [[mammaldb], None]),
                             ("Non-mammalian vertebrates", [[non_mammalvertdb], None]),
                             ("All vertebrates", [[mammaldb, non_mammalvertdb], None]),
                             ("Invertebrates", [[invertdb], None]),
                             ("Plants", [[plantdb], None]),
                             ("Fungi", [[fungidb], None]),
                             ("Complete", [[completedb], None])])
        dbstuff = db_version + upepsetting.DBRE[database_n]

        if seqqu:
            seqquery = seqqu.upper().replace("\n", "").replace("\r", "").replace(" ", "")
            if (seqquery[:3] == "NM_") or (seqquery[:3] == "XM_") or (seqquery[:3] == "GI:"):
                queryname = seqquery
                calcKaKs = True
                mRNAparams = [0, 0, 20]
                if minm:
                    mRNAparams[0] = int(minm)
                if maxm:
                    mRNAparams[1] = int(maxm)
                if gracem:
                    mRNAparams[2] = int(gracem)
                self.write_message("Started query through accession ID.")
                details = accessions.getmRNA(seqquery, dbstuff, codon, mRNAparams[0] * 3, mRNAparams[1] * 3,
                                             mRNAparams[2])
                if details:
                    alignquery = details[0]
                    blastquery = alignquery[
                                 :details[1][0] + mRNAparams[2]]  # send only the 5"UTR + grace of alignquery to tblastx.py
                    CDSfeature = details[1] + ["CDS"]

                    dbcon = MySQLdb.connect(user=dbuser, passwd=dbpass, host=dbhost, db=daba)
                    cursor = dbcon.cursor()
                    
                    sql_update_query_info = """UPDATE query_info SET seqquery = "%s", database_name = "%s", db_version = "%s", protein_name = "%s", initial_flag = 2 WHERE unique_id = "%s";"""
                    
                    cursor.execute(sql_update_query_info % (seqquery, database_n, db_version, details[3], unique_file_id))
                    dbcon.commit()
                    cursor.close()
                    dbcon.close()
                else:
                    dbcon = MySQLdb.connect(user=dbuser, passwd=dbpass, host=dbhost, db=daba)
                    cursor = dbcon.cursor()
                    sql_update_query_info = """UPDATE query_info SET seqquery = "%s", database_name = "%s", db_version = "%s", initial_flag = 1 WHERE unique_id = "%s";"""
                    cursor.execute(sql_update_query_info % (seqquery, database_n, db_version, unique_file_id))
                    dbcon.commit()
                    cursor.close()
                    dbcon.close()
                    self.write_message("Unable to find accession ID in database.")
            else:
                #print("<b>QUERY: User Entered</b>")
                self.write_message("Started query through user entered sequence.")
                queryname = "Query"
                dbcon = MySQLdb.connect(user=dbuser, passwd=dbpass, host=dbhost, db=daba)
                cursor = dbcon.cursor()
                sql_update_query_info = """UPDATE query_info SET seqquery = "Query", database_name = "%s", db_version = "%s", initial_flag = 3 WHERE unique_id = "%s";"""
                cursor.execute(sql_update_query_info % (database_n, db_version, unique_file_id))
                dbcon.commit()
                cursor.close()
                dbcon.close()

                alignquery = seqquery
                blastquery = alignquery
            ## TODO Un-hardcore these directories to make this code less terrible....
            if conserveduPEP.seqcheck(alignquery, ("A", "T", "C", "G")) is not None:
                changedseq = conserveduPEP.seqcheck(alignquery, ("A", "T", "C", "G"))
                #print(changedseq)
                #self.write_message("Checked for non-standard nucleotide.")
                if heatms == '':
                    dbcon = MySQLdb.connect(user=dbuser, passwd=dbpass, host=dbhost, db=daba)
                    cursor = dbcon.cursor()
                    sql_update_query_info = """UPDATE query_info SET heatmap_flag = 1 WHERE unique_id = "%s";"""
                    cursor.execute(sql_update_query_info % (unique_file_id))
                    dbcon.commit()
                    cursor.close()
                    dbcon.close()

                    return 0
                heatmapsize = int(heatms)
                print(heatmapsize)
                if (heatmapsize < 400) or (heatmapsize > 10000):
                    dbcon = MySQLdb.connect(user=dbuser, passwd=dbpass, host=dbhost, db=daba)
                    cursor = dbcon.cursor()
                    sql_update_query_info = """UPDATE query_info SET heatmap_flag = 2 WHERE unique_id = "%s";"""
                    cursor.execute(sql_update_query_info % (unique_file_id))
                    dbcon.commit()
                    cursor.close()
                    dbcon.close()

                    return 0
                blastmatches = tblastx.blastit(blastquery, dictionaries[database_n], data_loc)
                n = 0
                l = len(blastmatches)
                self.write_message("Blast search completed.")
                #print(blastmatches)
                if len(blastmatches) == 0:
                    dbcon = MySQLdb.connect(user=dbuser, passwd=dbpass, host=dbhost, db=daba)
                    cursor = dbcon.cursor()
                    sql_update_query_info = """UPDATE query_info SET blastmacthes_flag = 1 WHERE unique_id = "%s";"""
                    cursor.execute(sql_update_query_info % (unique_file_id))
                    dbcon.commit()
                    cursor.close()
                    dbcon.close()

                    return 0
                elif len(blastmatches) == 1:

                    dbcon = MySQLdb.connect(user=dbuser, passwd=dbpass, host=dbhost, db=daba)
                    cursor = dbcon.cursor()
                    sql_update_query_info = """UPDATE query_info SET blastmacthes_flag = 2 WHERE unique_id = "%s";"""

                    cursor.execute(sql_update_query_info % (unique_file_id))
                    dbcon.commit()
                    sql_insert_query_matchnumber = """INSERT INTO query_matchnumber (time_id, unique_id , match_number) VALUES (current_timestamp, "%s", "%i");"""
                    cursor.execute(sql_insert_query_matchnumber % (unique_file_id, len(blastmatches)))
                    dbcon.commit()
                    cursor.close()
                    dbcon.close()
                else:

                    dbcon = MySQLdb.connect(user=dbuser, passwd=dbpass, host=dbhost, db=daba)
                    cursor = dbcon.cursor()
                    sql_update_query_info = """UPDATE query_info SET blastmacthes_flag = 3 WHERE unique_id = "%s";"""
                    cursor.execute(sql_update_query_info % (unique_file_id))
                    dbcon.commit()
                    sql_insert_query_matchnumber = """INSERT INTO query_matchnumber (time_id, unique_id , match_number) VALUES (current_timestamp, "%s", "%i");"""
                    cursor.execute(sql_insert_query_matchnumber % (unique_file_id, len(blastmatches)))
                    dbcon.commit()
                    cursor.close()
                    dbcon.close()
                
                self.write_message("Analyzing each Blast result.")
                print(blastmatches[0])
                for blastmatch in blastmatches:
                    n += 1
                    hitdef = blastmatch[0][3:blastmatch[0].find("_", 6)]  ##Extract RefSeq Accession
                    self.write_message("%s" % (hitdef))
                    uPEPloc = eval(blastmatch[0][blastmatch[0].find("["):])  ##Extract uPEP location
                    if hitdef == seqquery:
                        self.write_message("- Trivial hit.")
                        dbcon = MySQLdb.connect(user=dbuser, passwd=dbpass, host=dbhost, db=daba)
                        cursor = dbcon.cursor()
                        sql_update_query_info = """UPDATE query_info SET trivial_flag = 1 WHERE unique_id = "%s";"""
                        cursor.execute(sql_update_query_info % (unique_file_id))
                        dbcon.commit()
                        sql_insert_query_trivial = """INSERT INTO query_trivial (time_id, unique_id, hitdef) VALUES (current_timestamp, "%s", "%s");"""
                        cursor.execute(sql_insert_query_trivial % (unique_file_id, hitdef))
                        dbcon.commit()
                        
                        cursor.close()
                        dbcon.close()
                        continue
                    details = accessions.getmRNA(hitdef, dbstuff, codon, 60, 300, 20)
                    if details:
                        alignref = details[0]

                    else:
                        self.write_message("Error: Internal uPEP database/Refseq database mismatch. " + hitdef)
                        sys.stderr.write("Error: Internal uPEP database/Refseq database mismatch. " + hitdef)
                        sys.stderr.flush()
                        continue

                    dbcon = MySQLdb.connect(user=dbuser, passwd=dbpass, host=dbhost, db=daba)
                    cursor = dbcon.cursor()
                    sql_insert_query_non_trvial = """INSERT INTO query_non_trivial (time_id, unique_id, protein_name, input_sequence, target_organism, input_hitdef, input_starting_position, input_ending_position, alignment, target_sequence, target_hitdef, target_starting_position, target_ending_position, upep_kaks, cds_kaks) VALUES (current_timestamp, "%s", "%s", "%s", "%s", "%s", "%i", "%i", "%s", "%s", "%s", "%i", "%i", 0, 0);"""
                    cursor.execute(sql_insert_query_non_trvial % (unique_file_id, details[3], blastmatch[5], details[4], queryname, blastmatch[2][0], blastmatch[2][1], blastmatch[7], blastmatch[6], hitdef, uPEPloc[0] + blastmatch[3][0] - 1, uPEPloc[0] + blastmatch[3][1] - 1))

                    dbcon.commit()
                    cursor.close()
                    dbcon.close()
                    ### GET DETAILS FOR KaKs -- yn00 ####
                    
                    if KaKs_t == 1 and calcKaKs:
                        self.write_message("- Calculating KaKs.")
                        dbcon = MySQLdb.connect(user=dbuser, passwd=dbpass, host=dbhost, db=daba)
                        cursor = dbcon.cursor()
                        sql_update_query_info = """UPDATE query_info SET CDSKaKs_flag = 1, uPEPKaKs_flag = 1 WHERE unique_id = "%s";"""
                        cursor.execute(sql_update_query_info % (unique_file_id))
                        dbcon.commit()
                        
                        cursor.close()
                        dbcon.close()

                        uPEPpair = [changedseq[blastmatch[2][0] - 1:blastmatch[2][1]],
                                    alignref[uPEPloc[0] + blastmatch[3][0] - 2: uPEPloc[0] + blastmatch[3][1] - 1]]
                        CDSpair = []
                        if CDSfeature:
                            CDSpair = [changedseq[CDSfeature[0] - 1:CDSfeature[1] - 3], alignref[
                                                                                        details[1][0] - 1:details[1][
                                                                                                              1] - 3]]  ### -1 to correct for python sequences // -3 to remove stop codon
                        if len(uPEPpair[0]) == len(uPEPpair[1]):
                            with open(tempfilename + "yn00", "wt") as yn00file:
                                yn00file.write("  2  %i\r\n" % (len(uPEPpair[0])))
                                yn00file.write("query\r\n" + uPEPpair[0] + "\r\n" + "ref\r\n" + uPEPpair[1] + "\r\n")
                            #yn00file = open(tempfilename + "yn00", "wt")

                            #try:
                                #yn00file.write("  2  %i\r\n" % (len(uPEPpair[0])))
                                #yn00file.write("query\r\n" + uPEPpair[0] + "\r\n" + "ref\r\n" + uPEPpair[1] + "\r\n")
                            #finally:
                                #yn00file.close()
                              
                            #retcode = subprocess.check_call([apps_loc+ "/yn00", tempfilename+ "yn00", tempfilename + "yn00" + "uPEP"])

                            retcode = subprocess.call([apps_loc+ "/yn00", tempfilename+ "yn00", tempfilename + "yn00" + "uPEP"], stdout=sys.stderr, stderr=sys.stderr)
                            
                            os.remove(tempfilename + "yn00")
                            kaksfile = open(tempfilename + "yn00" + "uPEP", "rt")
                            try:
                                try:
                                    for i in range(0, 6):
                                        line = kaksfile.readline().rstrip()
                                    temp = line.split(" ", 1)
                                    kaksuPEP = [temp[0]]
                                    kaksuPEP += temp[1][1:-1].split(" ")
                                except:
                                    sys.stderr.write("Line = {0}".format([temp]))
                                    sys.stderr.flush()
                            finally:
                                kaksfile.close()
                            os.remove(tempfilename + "yn00" + "uPEP")
                            KaKs_success = True
                            try:
                                Ka = abs(float(kaksuPEP[1]))
                                Ks = abs(float(kaksuPEP[2]))
                                KaKs = abs(float(kaksuPEP[0]))
                            except:
                                KaKs_success = False
                            if not KaKs_success:
                                n += 1

                                kkr = "Unable to estimate Ka/Ks ratio of uPEP."
                                dbcon = MySQLdb.connect(user=dbuser, passwd=dbpass, host=dbhost,
                                                                db=daba)
                                cursor = dbcon.cursor()
                                sql_update_query_non_trivial = (
                                """UPDATE query_non_trivial SET upep_kaks = "%i" WHERE (unique_id = "%s" and target_hitdef = "%s");""")
                                cursor.execute(sql_update_query_non_trivial % (n, unique_file_id, hitdef))
                                dbcon.commit()
                                sql_insert_query_uPEP_KaKs = (
                                    """INSERT INTO query_uPEP_KaKs (time_id, unique_id, hitdef, KaKs, upep_kaks_trigger) VALUES (current_timestamp, "%s", "%s", "%s", "%i");""")

                                cursor.execute(sql_insert_query_uPEP_KaKs % (unique_file_id, hitdef, kkr, n))
                                dbcon.commit()
                                cursor.close()
                                dbcon.close()
                            if Ks == 0:
                                n += 1

                                kkr = "Estimated uPEP Ka/Ks ratio: N/A (Ka: %.4f, Ks: %.4f)" % (
                                    Ka, Ks)
                                dbcon = MySQLdb.connect(user=dbuser, passwd=dbpass, host=dbhost,
                                                                db=daba)
                                cursor = dbcon.cursor()
                                sql_update_query_non_trivial = (
                                """UPDATE query_non_trivial SET upep_kaks = "%i" WHERE (unique_id = "%s" and target_hitdef = "%s");""")
                                cursor.execute(sql_update_query_non_trivial % (n, unique_file_id, hitdef))
                                dbcon.commit()
                                sql_insert_query_uPEP_KaKs = (
                                    """INSERT INTO query_uPEP_KaKs (time_id, unique_id, hitdef, KaKs, upep_kaks_trigger) VALUES (current_timestamp, "%s", "%s", "%s", "%i");""")

                                cursor.execute(sql_insert_query_uPEP_KaKs % (unique_file_id, hitdef, kkr, n))
                                dbcon.commit()
                                cursor.close()
                                dbcon.close()
                            else:
                                n += 1

                                kkr = "Estimated uPEP Ka/Ks ratio: %.4f (Ka: %.4f, Ks: %.4f)" % (
                                    KaKs, Ka, Ks)
                                dbcon = MySQLdb.connect(user=dbuser, passwd=dbpass, host=dbhost,
                                                                db=daba)
                                cursor = dbcon.cursor()
                                sql_update_query_non_trivial = (
                                """UPDATE query_non_trivial SET upep_kaks = "%i" WHERE (unique_id = "%s" and target_hitdef = "%s");""")
                                cursor.execute(sql_update_query_non_trivial % (n, unique_file_id, hitdef))
                                dbcon.commit()
                                sql_insert_query_uPEP_KaKs = (
                                    """INSERT INTO query_uPEP_KaKs (time_id, unique_id, hitdef, KaKs, upep_kaks_trigger) VALUES (current_timestamp, "%s", "%s", "%s", "%i");""")

                                cursor.execute(sql_insert_query_uPEP_KaKs % (unique_file_id, hitdef, kkr, n))
                                dbcon.commit()
                                cursor.close()
                                dbcon.close()
                        else:
                            n += 1
                            kkr = "Unable to estimate Ka/Ks ratio of uPEP: Query and Reference uPEPs of unequal size (%i vs. %i)." % (
                                len(uPEPpair[0]), len(uPEPpair[1]))
                            dbcon = MySQLdb.connect(user=dbuser, passwd=dbpass, host=dbhost,
                                                            db=daba)
                            cursor = dbcon.cursor()
                            sql_update_query_non_trivial = (
                            """UPDATE query_non_trivial SET upep_kaks = "%i" WHERE (unique_id = "%s" and target_hitdef = "%s");""")
                            cursor.execute(sql_update_query_non_trivial % (n, unique_file_id, hitdef))
                            dbcon.commit()
                            sql_insert_query_uPEP_KaKs = (
                                """INSERT INTO query_uPEP_KaKs (time_id, unique_id, hitdef, KaKs, upep_kaks_trigger) VALUES (current_timestamp, "%s", "%s", "%s", "%i");""")

                            cursor.execute(sql_insert_query_uPEP_KaKs % (unique_file_id, hitdef, kkr, n))
                            dbcon.commit()
                            cursor.close()
                            dbcon.close()
                        if CDSfeature:
                            self.write_message("- Find CDSfeature.")

                            if len(CDSpair[0]) == len(CDSpair[1]):
                                yn00file = open(tempfilename + "yn00", "wt")
                                try:
                                    yn00file.write("  2  %i\r\n" % (len(CDSpair[0])))
                                    yn00file.write("query\r\n" + CDSpair[0] + "\r\n" + "ref\r\n" + CDSpair[1] + "\r\n")
                                finally:
                                    yn00file.close()
                                retcode = subprocess.call(
                                    [apps_loc+ "/yn00", tempfilename + "yn00", tempfilename + "yn00" + "CDS"])
                                os.remove(tempfilename + "yn00")
                                try:
                                    kaksfile = open(tempfilename + "yn00" + "CDS")
                                    for i in range(0, 6):
                                        line = kaksfile.readline().rstrip()
                                    temp = line.split(" ", 1)
                                    kaksCDS = [temp[0]]
                                    kaksCDS += temp[1][1:-1].split(" ")
                                finally:
                                    kaksfile.close()
                                os.remove(tempfilename + "yn00" + "CDS")
                                KaKs_success = True
                                try:
                                    Ka = abs(float(kaksCDS[1]))
                                    Ks = abs(float(kaksCDS[2]))
                                    KaKs = abs(float(kaksCDS[0]))
                                except:
                                    KaKs_success = False
                                if not KaKs_success:
                                    n += 1
                                    kkr = "Unable to estimate Ka/Ks ratio of CDS."
                                    dbcon = MySQLdb.connect(user=dbuser, passwd=dbpass, host=dbhost,
                                                                    db=daba)
                                    cursor = dbcon.cursor()
                                    sql_update_query_non_trivial = (
                                    """UPDATE query_non_trivial SET cds_kaks = "%i" WHERE (unique_id = "%s" and target_hitdef = "%s");""")
                                    cursor.execute(sql_update_query_non_trivial % (n, unique_file_id, hitdef))
                                    dbcon.commit()
                                    sql_insert_query_CDS_KaKs = (
                                        """INSERT INTO query_CDS_KaKs (time_id, unique_id, hitdef, KaKs, cds_kaks_trigger) VALUES (current_timestamp, "%s", "%s", "%s", "%i");""")
                                    cursor.execute(sql_insert_query_CDS_KaKs % (unique_file_id, hitdef, kkr, n))
                                    dbcon.commit()
                                    cursor.close()
                                    dbcon.close()
                                if Ks == 0:
                                    n += 1

                                    kkr = "Estimated CDS Ka/Ks ratio:&nbsp;&nbsp;N/A (Ka: %.4f, Ks: %.4f)" % (
                                        Ka, Ks)
                                    dbcon = MySQLdb.connect(user=dbuser, passwd=dbpass, host=dbhost,
                                                                    db=daba)
                                    cursor = dbcon.cursor()
                                    sql_update_query_non_trivial = (
                                    """UPDATE query_non_trivial SET cds_kaks = "%i" WHERE (unique_id = "%s" and target_hitdef = "%s");""")
                                    cursor.execute(sql_update_query_non_trivial % (n, unique_file_id, hitdef))
                                    dbcon.commit()
                                    sql_insert_query_CDS_KaKs = (
                                        """INSERT INTO query_CDS_KaKs (time_id, unique_id, hitdef, KaKs, cds_kaks_trigger) VALUES (current_timestamp, "%s", "%s", "%s", "%i");""")
                                    cursor.execute(sql_insert_query_CDS_KaKs % (unique_file_id, hitdef, kkr, n))
                                    dbcon.commit()
                                    cursor.close()
                                    dbcon.close()
                                else:
                                    n += 1

                                    kkr = "Estimated CDS Ka/Ks ratio:&nbsp;&nbsp;%.4f (Ka: %.4f, Ks: %.4f)" % (
                                        KaKs, Ka, Ks)
                                    dbcon = MySQLdb.connect(user=dbuser, passwd=dbpass, host=dbhost,
                                                                    db=daba)
                                    cursor = dbcon.cursor()
                                    sql_update_query_non_trivial = (
                                    """UPDATE query_non_trivial SET cds_kaks = "%i" WHERE (unique_id = "%s" and target_hitdef = "%s");""")
                                    cursor.execute(sql_update_query_non_trivial % (n, unique_file_id, hitdef))
                                    dbcon.commit()
                                    sql_insert_query_CDS_KaKs = (
                                        """INSERT INTO query_CDS_KaKs (time_id, unique_id, hitdef, KaKs, cds_kaks_trigger) VALUES (current_timestamp, "%s", "%s", "%s", "%i");""")
                                    cursor.execute(sql_insert_query_CDS_KaKs % (unique_file_id, hitdef, kkr, n))
                                    dbcon.commit()
                                    cursor.close()
                                    dbcon.close()
                            else:
                                n += 1

                                kkr = "Unable to estimate Ka/Ks ratio of CDS: Query and reference coding sequences of unequal size (%i vs. %i)." % (
                                    len(CDSpair[0]), len(CDSpair[1]))
                                dbcon = MySQLdb.connect(user=dbuser, passwd=dbpass, host=dbhost,
                                                                db=daba)
                                cursor = dbcon.cursor()
                                sql_update_query_non_trivial = (
                                """UPDATE query_non_trivial SET cds_kaks = "%i" WHERE (unique_id = "%s" and target_hitdef = "%s");""")
                                cursor.execute(sql_update_query_non_trivial % (n, unique_file_id, hitdef))
                                dbcon.commit()
                                sql_insert_query_CDS_KaKs = (
                                    """INSERT INTO query_CDS_KaKs (time_id, unique_id, hitdef, KaKs, cds_kaks_trigger) VALUES (current_timestamp, "%s", "%s", "%s", "%i");""")
                                cursor.execute(sql_insert_query_CDS_KaKs % (unique_file_id, hitdef, kkr, n))
                                dbcon.commit()
                                cursor.close()
                                dbcon.close()
                        else:
                            n += 1

                            kkr = "Unable to estimate Ka/Ks ratio of CDS: Unable to define CDS from user entered sequence."
                            dbcon = MySQLdb.connect(user=dbuser, passwd=dbpass, host=dbhost,
                                                            db=daba)
                            cursor = dbcon.cursor()
                            sql_update_query_non_trivial = (
                            """UPDATE query_non_trivial SET cds_kaks = "%i" WHERE (unique_id = "%s" and target_hitdef = "%s");""")
                            cursor.execute(sql_update_query_non_trivial % (n, unique_file_id, hitdef))
                            dbcon.commit()
                            sql_insert_query_CDS_KaKs = (
                                """INSERT INTO query_CDS_KaKs (time_id, unique_id, hitdef, KaKs, cds_kaks_trigger) VALUES (current_timestamp, "%s", "%s", "%s", "%i");""")
                            cursor.execute(sql_insert_query_CDS_KaKs % (unique_file_id, hitdef, kkr, n))
                            dbcon.commit()
                            cursor.close()
                            dbcon.close()
                    features = []
                    reffeatures = []
                    reffeatures.append(details[1] + ["CDS"])
                    for j in details[2]:
                        reffeatures.append(j + ["uPEP"])
                    if CDSfeature:
                        features.append(CDSfeature)
                    features.append(blastmatch[2] + ["uPEP"])
                    ##if True: #i[3][0] == 1:
                    ##    retcode = subprocess.call(["./yn00", "ptp4a1upep.nuc", "ptpuPEP"])
                    self.write_message("- Alignment using LAGAN.")
                    alignfile1 = open(tempfilename + hitdef+"_1", "wt")
                    alignfile2 = open(tempfilename + hitdef+"_2", "wt")
                    try:
                        alignfile1.write(">Query\n" + changedseq + "\n")
                        alignfile2.write(">Reference\n" + alignref + "\n")
                    finally:
                        alignfile1.close()
                        alignfile2.close()
                    os.environ["LAGAN_DIR"] = lagan_loc
                    retcode = subprocess.call(
                        ["perl", lagan_loc+ "/lagan.pl", tempfilename + hitdef+"_1", tempfilename + hitdef+"_2", "-out",
                         tempfilename + hitdef+".mfa", "-mfa"])
                    retcode = subprocess.call(
                        ["perl", lagan_loc+ "/mf_to_align.pl", "-f", tempfilename + hitdef+".mfa", "-out",
                         tempfilename + hitdef+"aligneduf"])
                    window = window_t
                    windowfilename = tempfilename + hitdef+"w"+window
                    with open(tempfilename + hitdef+"aligneduf", "rt") as aligneduffile, open(tempfilename + hitdef+"aligned", "wt") as alignedfile: 
                        query = "".join(aligneduffile.readline().rstrip("\r\n"))
                        match = "".join(aligneduffile.readline().rstrip("\r\n"))
                        ref = "".join(aligneduffile.readline().rstrip("\r\n"))
                        query2 = conserveduPEP.split_sequence(query)
                        match2 = conserveduPEP.split_sequence(match)
                        ref2 = conserveduPEP.split_sequence(ref)
                        for i1, i2, i3 in zip(query2, match2, ref2):
                            alignedfile.write("%s\n%s\n%s\n\n" % (i1, i2, i3))
                        
                        
                        dbcon = MySQLdb.connect(user=dbuser, passwd=dbpass, host=dbhost, db=daba)
                        cursor = dbcon.cursor()
                        sql_update_query_info = """UPDATE query_info SET window = "%i" WHERE unique_id = "%s";"""
                        cursor.execute(sql_update_query_info % (int(window), unique_file_id))
                        dbcon.commit()
                        cursor.close()
                        dbcon.close()
                        
                        
                        conserveduPEP.makeheatmap(windowfilename + ".ppm", conserveduPEP.homologyvector(query, match, int(window)), features, grad, heatmapsize, 25)
                        if refhm == 1:
                            conserveduPEP.makeheatmap(windowfilename + "r.ppm", conserveduPEP.homologyvector(ref, match, int(window)), reffeatures, grad, heatmapsize, 25)
                        
                    
                    self.write_message("- Draw Heatmap.")
                        
                    if refhm == 1:
                        retcode = subprocess.call([apps_loc+ "/png.py", "-c", "9", windowfilename + "r.ppm"])
                        os.remove(windowfilename + "r.ppm")
                    
                    retcode = subprocess.call([apps_loc+ "/png.py", "-c", "9", windowfilename + ".ppm"])
                    os.remove(windowfilename + ".ppm")

                    if refhm == 1:
                        dbcon = MySQLdb.connect(user=dbuser, passwd=dbpass, host=dbhost, db=daba)
                        cursor = dbcon.cursor()
                        sql_update_query_info = """UPDATE query_info SET Refheatmap = 1 WHERE unique_id = "%s";"""
                        cursor.execute(sql_update_query_info % (unique_file_id))
                        dbcon.commit()
                        cursor.close()
                        dbcon.close()
                self.write_message('Complete.')

            else:
                errorstuff = "align error"


class MainHandler(tornado.web.RequestHandler):
    def get(self):
        self.render("index.html", title="Home Page",)
class uPEPHelperForm(tornado.web.RequestHandler):
    def get(self):
        self.render("uhform.html", title="uPEP Database Updater", codons=helpersetting.STARTING_CODONS, databases=helpersetting.REFSEQ_DBS,)

class uPEPHelperProcess(tornado.web.RequestHandler):
    @gen.coroutine
    def post(self):
        bdatabase = self.request.body_arguments["database"][0]
        database = bdatabase.decode("utf-8")
        #databases = helpersetting.REFSEQ_DBS     
        codons = helpersetting.STARTING_CODONS
        input_codons = []
        for codon in codons:
            if codon in self.request.body_arguments:
               input_codons.append(codon)
        override = ""
        if self.request.body_arguments["override"]:
            override = "True"
        #pool = ThreadPoolExecutor(max_workers=settingmain.MAX_WORKERS)
        #yield pool.submit(helpertasks.upephelper_processing, database, input_codons, override)
        print(database)
        print(input_codons)
        print(override)
        q = Queue("upephelper", connection=Redis())
        result = q.enqueue(helpertasks.upephelper_processing, args=(database, input_codons, override,), timeout=343600)
        print(result.id)
        self.render("results_helper.html", title="uPEP Database Updater Result",)
                

class uPEPForm(tornado.web.RequestHandler):
    def get(self):
        self.render("upep.html", title="uPEPperoni - An Online Tool for uPEP Location and Analysis", codons=upepsetting.STARTING_CODONS, databases=upepsetting.DBWEB,)

class uPEPProcess(tornado.web.RequestHandler):
    def post(self):
        parameters = {}
        parameters["timeid"] = time.time()
        parameters["unid"] = "%012x%016x" % (int(parameters["timeid"] * 1000), random.randint(0, 0xFFFFFFFFFFFFFFFF))
        bdb = self.request.body_arguments["database"][0]
        parameters["db"] = bdb.decode("utf-8")
        
        bcodon = self.request.body_arguments["codon"][0]
        parameters["codon"] = bcodon.decode("utf-8")
        bdbv = self.request.body_arguments["database_version"][0]
        parameters["dbv"] = bdbv.decode("utf-8")
        bfilter = self.request.body_arguments["filter"][0]
        parameters["filter"] = bfilter.decode("utf-8")
        parameters["kaks"] = 0
        if "kaks" in self.request.body_arguments:
            parameters["kaks"] = 1
        parameters["rhm"] = 0
        if "refheatmaps" in self.request.body_arguments:
            parameters["rhm"] = 1
        parameters["strict"] = 0
        if "strict" in self.request.body_arguments:
            parameters["strict"] = 1
        parameters["match"] = 0
        if "match" in self.request.body_arguments:
            parameters["match"] = 1
        parameters["mismatch"] = 0
        if "mismatch" in self.request.body_arguments:
            parameters["mismatch"] = 1
        parameters["existence"] = 0
        if "existence" in self.request.body_arguments:
            parameters["existence"] = 1
        parameters["existension"] = 0
        if "extension" in self.request.body_arguments:
            parameters["extension"] = 1
        parameters["min"] = 0
        if "min" in self.request.body_arguments:
            bmin = self.request.body_arguments["min"][0]
            parameters["min"] = bmin.decode("utf-8")
        parameters["max"] = 0
        if "min" in self.request.body_arguments:
            bmax = self.request.body_arguments["max"][0]
            parameters["max"] = bmax.decode("utf-8")
        parameters["grace"] = 0
        if "min" in self.request.body_arguments:
            bgrace = self.request.body_arguments["grace"][0]
            parameters["grace"] = bgrace.decode("utf-8")
        gradi = []
        gradnames = ["black", "blue", "green", "yellow", "red"]
        for i in range(0, 5):
            if "g" + gradnames[i] in self.request.body_arguments:
                if not self.request.body_arguments["g" + gradnames[i]][0] == b"":
                    bcolor = self.request.body_arguments["g" + gradnames[i]][0]
                    color = bcolor.decode("utf-8")
                    gradi.append(float(color) / 100)
                else:
                    gradi.append(None)
        parameters["gradi"] = gradi
        parameters["hms"] = 0
        if "heatmapsize" in self.request.body_arguments:
            bhms = self.request.body_arguments["heatmapsize"][0]
            parameters["hms"] = bhms.decode("utf-8")
        parameters["window"] = 0
        if "window" in self.request.body_arguments:
            bws = self.request.body_arguments["window"][0]
            parameters["window"] = bws.decode("utf-8")
 
        parameters["seqquery"] = ""
        if "seqquery" in self.request.body_arguments:
            dbuser = upepsetting.DATABASES["default"]["USER"]
            dbpass = str(upepsetting.DATABASES["default"]["PASSWORD"])
            dbhost = upepsetting.DATABASES["default"]["HOST"]
            daba = upepsetting.DATABASES["default"]["DB"]
            bseqquery = self.request.body_arguments["seqquery"][0]
            seqquery = bseqquery.decode("utf-8")
            parameters["seqquery"] = seqquery.upper().replace("\n", "").replace("\r", "").replace(" ", "")
            if len(parameters["seqquery"]) == 28:
                dbcon = MySQLdb.connect(user=dbuser, passwd=dbpass, host=dbhost, db=daba)
                cursor = dbcon.cursor()
                sql_retrieve_query = """select * from query_info where unique_id = "%s";"""
                cursor.execute(sql_retrieve_query % parameters["seqquery"])

                if not cursor.rowcount:
                    cursor.close()
                    dbcon.close()
                    
                    cursor.close()
                    dbcon.close()
                    parameters['task'] = 1
                    self.render("upep_process.html", title="uPEP Process", parameters=json.dumps(parameters), filter_org = parameters['filter'], unid = parameters['unid'], strict = str(parameters['strict']),)        
                else:
                    result = cursor.fetchone()
                    cursor.close()
                    dbcon.close()
                    parameters['task'] = 0
                    self.render("upep_process.html", title="uPEP Process", task=0, parameters=json.dumps(parameters), filter_org = parameters['filter'], unid = result[1], strict = str(parameters['strict']),)
            else:
                parameters['task'] = 1
                self.render("upep_process.html", title="uPEP Process", task=1, parameters=json.dumps(parameters), filter_org = parameters['filter'], unid = parameters['unid'], strict = str(parameters['strict']),)
                    
class GetResult(tornado.web.RequestHandler):
    def get(self, job_key, filter_org, strict):
        self.weboutput(job_key, filter_org, int(strict))
        
    def weboutput(self, unique_file_id, filter_org, strict):
        
        dbuser = upepsetting.DATABASES['default']['USER']
        dbpass = str(upepsetting.DATABASES['default']['PASSWORD'])
        dbhost = upepsetting.DATABASES['default']['HOST']
        daba = upepsetting.DATABASES['default']['DB']
        dbcon = MySQLdb.connect(user=dbuser, passwd=dbpass, host=dbhost, db=daba)
        cursor = dbcon.cursor()

        seqq, dbn, prn, dbv, tfn, ini_flag, trv_flag, bm_flag, hm_flag, cdsk_flag, \
        upepk_flag, win_r, rhm, saa = conserveduPEP.parsing_initial(unique_file_id)
        
        mnumber = 0
        webresult = 0
        trivial_row = 0
        ms = ''
        upepkaks_dict = {}
        upepcds_dict = {}
        if ini_flag == 1:
            # print 'Error: Unable to find query mRNA sequence in selected Refseq database.' + seqq + ' not found \n \
            print('No results were recorded. \n \
                    Click <a href="javascript:history.go(-1)"> here </a> to return to the input page.')
            print('</body></html>')
            cursor.close()
            dbcon.close()
            self.render("upep_result.html", title="uPEP Result", page_er = 'Error: Unable to find query mRNA sequence in selected Refseq database.\n'+seqq+' not found\n', seqq=seqq, dbn=dbn, prn=prn, dbv=dbv, tfn=tfn, ini_flag = ini_flag, trv_flag = trv_flag, bm_flag = bm_flag, hm_flag = hm_flag, cdsk_flag = cdsk_flag, upepk_flag = upepk_flag, win_r = win_r, rhm = rhm, saa = saa)
            
        if ini_flag == 2:
            print('<b>QUERY: ' + prn + ' (' + seqq + ')</b>')
        elif ini_flag == 3:
            print('<b>QUERY: User Entered</b>')
        if hm_flag == 0:
            if bm_flag == 1:
                print('<br>')
                # print 'No hits were found for the given input sequence in the \'' + dbn + '\' uPEP database. \n \
                print('No results were recorded. \n \
                Click <a href="javascript:history.go(-1)"> here </a> to return to the input page.')
                print('</body></html>')
                cursor.close()
                dbcon.close()
                self.render("upep_result.html", title="uPEP Result", page_er = 'No hits were found for the given input sequence in the '+dbn+' uPEP database. \n', seqq=seqq, dbn=dbn, prn=prn, dbv=dbv, tfn=tfn, ini_flag = ini_flag, trv_flag = trv_flag, bm_flag = bm_flag, hm_flag = hm_flag, cdsk_flag = cdsk_flag, upepk_flag = upepk_flag, win_r = win_r, rhm = rhm, saa = saa)
                

            elif bm_flag == 2:
                print('<br><b> 1 hit found in raw output.</b><br><br>')
            elif bm_flag == 3:

                sql_bm_retrieve = """select * from query_matchnumber where unique_id = "%s";"""
                cursor.execute(sql_bm_retrieve % unique_file_id)

                for (time_id, unique_id, match_number) in cursor:
                    mnumber = match_number
                    print('<br><b> %i hits found in raw output.</b><br><br>' % mnumber)
            if trv_flag == 1:
                sql_trv_retrieve = """select * from query_trivial where unique_id = "%s";"""
                cursor.execute(sql_trv_retrieve % unique_file_id)

                result = cursor.fetchall()
                trivial_row = cursor.rowcount
                print('<b> %i hit(s) with the same transcript as query sequence </b><br><br>' % trivial_row)
            sql_non_trivial_retrieve, match_string = conserveduPEP.strictmode_and_filter_trigger(unique_file_id, strict, filter_org, saa)

            cursor.execute(sql_non_trivial_retrieve)
            if not cursor.rowcount:
                print('<br>')
                print('Error: Internal uPEP database/Refseq database mismatch. Query not found \n \
                This error has been logged and sent to the webmaster. See help for further details.')
                cursor.close()
                dbcon.close()
                self.render("upep_result.html", title="uPEP Result", page_er = 'Error: Internal uPEP database/Refseq database mismatch. Query not found.\n', seqq=seqq, dbn=dbn, prn=prn, dbv=dbv, tfn=tfn, ini_flag = ini_flag, trv_flag = trv_flag, bm_flag = bm_flag, hm_flag = hm_flag, cdsk_flag = cdsk_flag, upepk_flag = upepk_flag, win_r = win_r, rhm = rhm, saa = saa, mnumber = mnumber, trivial_row = trivial_row)
                return 0

            else:

                result = cursor.fetchall()
                webresult = result
                if cursor.rowcount > 0:
                    ms = match_string % (cursor.rowcount)
                    print(match_string % (cursor.rowcount))

                    for row in result:
                        print('<b><u>HIT/REFERENCE:</u> ' + row[2] + ' (' + row[10] + ')</b><br>')
                        print('<PRE>%s   %s, [%i, %i]\n%s\n%s   %s, [%i, %i]</PRE>' % (
                            row[3], seqq, row[6], row[7], row[8], row[9], row[10],
                            row[11], row[12]))
                        if upepk_flag == 1:
                            sql_upepk_retrieve = """select * from query_uPEP_KaKs where (unique_id = "%s" and hitdef = "%s" and upep_kaks_trigger = "%i");"""
                            cursor.execute(sql_upepk_retrieve % (unique_file_id, row[10], row[13]))
                            upepk_result = cursor.fetchall()
                            for i in upepk_result:
                                upepkaks_dict[row[10]] = i[3]
                                print('<span style="font-size:80%"><b>' + i[3] + '</b></span><br>')
                        if cdsk_flag == 1:
                            sql_cdsk_retrieve = """select * from query_CDS_KaKs where (unique_id = "%s" and hitdef = "%s" and cds_kaks_trigger = "%i");"""
                            cursor.execute(sql_cdsk_retrieve % (unique_file_id, row[10], row[14]))
                            cdsk_result = cursor.fetchall()
                            for i in cdsk_result:
                                upepcds_dict[row[10]] = i[3]
                                print('<span style="font-size:80%"><b>' + i[3] + '</b></span><br>')
                        print('The unformatted aligned sequence can be viewed <a href=' + tfn[8:] + row[
                            10] + 'aligned>here.</a><br>')
                        wfn = tfn + row[10] + 'w' + str(win_r)
                        print('<br><span style="font-size:80%"><b>Heatmap representation of ' + seqq + ':</b></span><br>')
                        print('<br><img src=' + wfn[8:] + '.png' + '><br><br>')
                        if rhm == 1:
                            print('<br><span style="font-size:80%"><b>Heatmap representation of ' + row[10] + ':</b></span><br>')
                            print('<br><img src=' + wfn[8:] + 'r.png' + '><br><br>')
                    print('The heatmap diagrams can be downloaded by right-clicking on the image and \
                                selecting \"Save Picture As...\" (or \"Save Image As...\" in Firefox)')
                    print('</body></html>')
                    cursor.close()
                    dbcon.close()
                else:
                    print('<br><b>No hits according to selected parameters.</b><br>')
            self.render("upep_result.html", title="uPEP Result", seqq=seqq, dbn=dbn, prn=prn, dbv=dbv, tfn=tfn, ini_flag = ini_flag, trv_flag = trv_flag, bm_flag = bm_flag, hm_flag = hm_flag, cdsk_flag = cdsk_flag, upepk_flag = upepk_flag, win_r = win_r, rhm = rhm, saa = saa, mnumber = mnumber, webresult = webresult, trivial_row = trivial_row, ms = ms, upepkaks_dict = upepkaks_dict, upepcds_dict = upepcds_dict, unid = unique_file_id)
        else:
            print('<br>')
            print('No results was recorded. Click <a href="javascript:history.go(-1)"> here </a> to return to the input page.')
            print('</body></html>')
            cursor.close()
            dbcon.close()
            return 0
class WHMForm(tornado.web.RequestHandler):
    def get(self):
        self.render("whm.html", title="uPEP Webheatmap Utility", codons=upepsetting.STARTING_CODONS,)

class WHMProcess(tornado.web.RequestHandler):
    def post(self):
        results = whm.whmstandalone(self.request.body_arguments)
        
        print(results)
        self.render("whm_result.html", title="uPEP Webheatmap Result", results=results,)
        
settings = {
    "autoreload": True,
    "debug": True,
    "static_path": settingmain.APP_STATIC,
    "template_path": settingmain.APP_TEMPLATE,
    }

if __name__ == "__main__":
    application = tornado.web.Application([
        (r"/", MainHandler),
        (r"/ws", wsUI),
        (r"/websocket", WebSocketHandler),
        (r"/upephelper", uPEPHelperForm),
        (r"/upephelper/process", uPEPHelperProcess),
        (r"/upep/result/data/(.*)", tornado.web.StaticFileHandler, {"path": "./upep/data"}),
        (r"/upep/result/(?P<job_key>[^\/]+)/(?P<filter_org>[^\/]+)/(?P<strict>[^\/]+)?/", GetResult),
        (r"/upep", uPEPForm),
        (r"/upep/process", uPEPProcess),
        (r"/whm", WHMForm),
        (r"/whm/process", WHMProcess),
    ], **settings)
    application.listen(9999)
    tornado.ioloop.IOLoop.current().start()
