#!/usr/bin/env python3
"""
fabfile for doing uPEP things

Mitchell Stanton-Cook
m.stantoncook@gmail.com

TODO: No real verification of the successful completion of step
(i.e. verify all files downloaded, verify compaction, compilation, uPEPFinding
& BLAST DB creation). 

"""
import os
import sys
import gzip
import bisect
import glob
import ast
import re
#import configparser
import ftputil
#from fabric.api import task
#import MySQLdb
from . import helpersetting
import tqdm
import MySQLdb

HOST = "ftp.ncbi.nlm.nih.gov"
BASE = "/refseq/release/"
#config = configparser.ConfigParser()
#config.read('config.ini')

database_loc = helpersetting.UPEPHELPER_DATABASE #config.get('Database', 'database_ext')
staging_loc = helpersetting.UPEPHELPER_STAGING #config.get('Database', 'staging')


def setup(outpath, home, db):
    """
    Build outpath and chdir to it

    This is relative to the fabfile path!

    :param outpath: the outpath
    :param home: the directory this fabfile is executed in
    :param db: the databse identifier

    :rtype: the fullpath to the outpath as a string
    """
    if not outpath:
        outpath = os.path.join(home, db)
    else:
        outpath = os.path.join(os.path.expanduser(outpath), db)
    try:
        os.mkdir(outpath)
    except OSError:
        print("Output directory exists")
        pass
    os.chdir(outpath)
    return outpath


def download_db(db):
    """
    Download all RefSeq database files to the current working driectory

    :param db: the database identifier
    """
    name = db.split('-')[-1]
    url = BASE + name
    host = ftputil.FTPHost(HOST, 'anonymous', 'beaton.lab@gmail.com')
    host.chdir(url)
    remote = host.listdir(host.curdir)
    n = 0
    for f in remote:
        if f.endswith("rna.gbff.gz"):
            print('Retrieving %s' % f)
            n += 1
            host.download(f, f)
    return n

#@task
def get_NCBI_RefSeq_release():
    """
    Prints & returns the NCBI RefSeq release number
    
    :rtype: the release number as an integer
    """
    url = BASE + "release-notes/"
    host = ftputil.FTPHost(HOST, 'anonymous', 'beaton.lab@gmail.com')
    host.chdir(url)
    remote = host.listdir(host.curdir)
    number = -1
    for f in remote:
        if f.startswith("RefSeq-"):
            number = f.split('release')[-1][:-4]
            break
    print('Remote Release %s' % (number))
    return int(number)


#@task
def get_uPEP_RefSeq_release(location=database_loc+'.refseq_version'):

    
    """
    Prints & returns the uPEP RefSeq release number
    
    :param location: [default=/var/RefSeq/.refseq_version] the file that the 
                     RefSeq version is stored in

    :rtype: the release number as an integer
    """
    number = -1
    with open(location) as fin:
        number = fin.readlines()[0].strip()
    print('Local Release %s' % (number))
    return int(number)

def upep_mysql_database(unid, dbuser, dbpass, dbhost, db, key, codonstring, override, dbv):
    dbcon = MySQLdb.connect(user=dbuser, passwd=dbpass, host=dbhost, db=db)
    cursor = dbcon.cursor()
    
    #sql_dbu = ("use uPEP;")
    
    #cursor.execute(sql_dbu)
    
    sql_tbc = ("create table if not exists query_info (time_id TIMESTAMP NOT NULL, unique_id VARCHAR(40) NOT NULL,"
               "seqquery VARCHAR(50) NOT NULL, database_name VARCHAR(20) NOT NULL, db_version VARCHAR(4) NOT NULL, protein_name VARCHAR(1000) NOT NULL, tempfilename VARCHAR(50) NOT NULL,"
               "initial_flag TINYINT NOT NULL, trivial_flag TINYINT NOT NULL, blastmacthes_flag TINYINT NOT NULL, CDSKaKs_flag TINYINT NOT NULL, uPEPKaKs_flag TINYINT NOT NULL, heatmap_flag TINYINT NOT NULL, window INT NOT NULL, Refheatmap TINYINT NOT NULL);")
               
    cursor.execute(sql_tbc)
    
    
    sql_tbc2 = ("create table if not exists query_trivial (time_id TIMESTAMP NOT NULL, unique_id VARCHAR(40) NOT NULL,"
                "hitdef VARCHAR(20) NOT NULL);")
    
    cursor.execute(sql_tbc2)
    
    
    sql_tbc3 = ("create table if not exists query_matchnumber (time_id TIMESTAMP NOT NULL, unique_id VARCHAR(40) NOT NULL, match_number INT NOT NULL);")

    cursor.execute(sql_tbc3)
    
    sql_tbc4 = ("create table if not exists query_non_trivial (time_id TIMESTAMP NOT NULL, unique_id VARCHAR(40) NOT NULL,"
                "protein_name VARCHAR(1000) NOT NULL, input_sequence VARCHAR(200) NOT NULL, target_organism VARCHAR(50) NOT NULL, input_hitdef VARCHAR(20) NOT NULL,"
                "input_starting_position INT NOT NULL, input_ending_position INT NOT NULL, alignment VARCHAR(200) NOT NULL, target_sequence VARCHAR(200) NOT NULL,"
                "target_hitdef VARCHAR(20) NOT NULL, target_starting_position INT NOT NULL, target_ending_position INT NOT NULL);")
                
    cursor.execute(sql_tbc4)
    
    
    sql_tbc5 = ("create table if not exists query_uPEP_KaKs (time_id TIMESTAMP NOT NULL, unique_id VARCHAR(40) NOT NULL, hitdef VARCHAR(20) NOT NULL, KaKs VARCHAR(200) NOT NULL, upep_kaks_trigger INT NOT NULL);")
    cursor.execute(sql_tbc5)

    
    sql_tbc6 = ("create table if not exists query_CDS_KaKs (time_id TIMESTAMP NOT NULL, unique_id VARCHAR(40) NOT NULL, hitdef VARCHAR(20) NOT NULL, KaKs VARCHAR(200) NOT NULL, cds_kaks_trigger INT NOT NULL);")
       
    cursor.execute(sql_tbc6)
    
    sql_tbc7 = ("create table if not exists updater_log (time_id_start TIMESTAMP NOT NULL, time_id_finish TIMESTAMP NOT NULL, unique_id VARCHAR(40) NOT NULL, refseq_database VARCHAR(20) NOT NULL, refseq_database_version INT NOT NULL, acc_db TINYINT NOT NULL, gi_db TINYINT NOT NULL, codons VARCHAR(100) NOT NULL, override_refseq_database VARCHAR(5) NOT NULL, success_log TINYINT NOT NULL);")
    cursor.execute(sql_tbc7)
    
    update_log = ("""INSERT INTO updater_log (time_id_start, time_id_finish, unique_id, refseq_database, refseq_database_version, acc_db, gi_db, codons, override_refseq_database, success_log) VALUES (current_timestamp, current_timestamp, "%s", "%s", "%i", 1, 1, "%s", "%s", 0);""")
    cursor.execute(update_log % (unid, key, dbv, codonstring, override))
    dbcon.commit()
    cursor.close()
    dbcon.close()
    
#@task
def build_upep_dbs(outpath=staging_loc, key=None, override=False):

    """
    Update the uPEP databases

    Alternatively if a key is given (one of):
        * RefSeq-complete
        * RefSeq-fungi
        * RefSeq-invertebrate
        * RefSeq-plant
        * RefSeq-vertebrate_mammalian
        * RefSeq-vertebrate_other
    
    the task will only upgrade the givenDB.

    You can set the outpath with outpath='/dump/me/here'

    Setting override to True will update even in local version is the same as
    the NCBI RefSeq version

    This task:
        * checks that you need to update (unless overridden)
        * sets up a download directory
        * fetches all required RefSeq databases from NCBI
        * compacts
        * compiles
        * uPEP finds
        * converts to BLAST db
        * moves the files to the production server (creating correct 
          permissions, symbolic links etc)
        * stores the new local RefSeq version

    :param outpath: the base location to dump the files to (must exist)
                    i.e. /var/RefSeq/staging (full path as a string). If not 
                    given will be dumped in the fabfile directory
    :param key: [def=None] a specific database to grab. If not given will grab 
                all
    :param override: [default = False] upgrade even if local and remote are 
                     the same version
    """
    dbuser = settings.DATABASES['default']['USER']
    dbpass = str(settings.DATABASES['default']['PASSWORD'])
    dbhost = settings.DATABASES['default']['HOST']
    override = ast.literal_eval(str(override))
    remote = get_NCBI_RefSeq_release()
    dbv = str(remote)
    local = get_uPEP_RefSeq_release()
    if remote > local or override == True:
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
                wd = setup(outpath, home, key)
                download_db(key)
                compacted = compact_RefSeq(wd, dbv)
                compile_RefSeq(compacted, dbv, dbuser, dbpass, dbhost)
                os.chdir(home)
            else:
                print("Not a defined db")
                sys.exit(1)
        else:
            for db in dbs:
                print("Working with database" + db)
                wd = setup(outpath, home, db)
                download_db(db)
                compacted = compact_RefSeq(wd, dbv)
                compile_RefSeq(compacted, dbv, dbuser, dbpass, dbhost)
                os.chdir(home)
        os.chdir(outpath)
        try:
            os.mkdir("../tmp/RefSeqdb"+dbv)
        except OSError:
            print("RefSeqdb %s directory exists \n Overriding directory" % dbv)
            pass
        
        try:
            os.system("mv "+staging_loc+"RefSeq* ../tmp/RefSeqdb"+dbv)
        except OSError:
            os.system("cp -r "+staging_loc+"RefSeq* ../tmp/RefSeqdb"+dbv)
            os.system("rm -rf "+staging_loc+"RefSeq*")
        
        for starting_codon in starting_codons:
            _, proc_list = uPEP_finder(codon=starting_codon, db_version=dbv, outpath=outpath)
            build_blast_db(proc_list)
        # Here is where we should go into maintainence mode
        ioloop.IOLoop.current().run_sync(upep_mysql_database(dbuser, dbpass, dbhost)) 
        finalise_update()
        
    else:
        print("No updrade required")


def compact_RefSeq(RefSeq_directory, db_version):
    """
    Removes 'junk' from the RefSeq files

    Provided by Adam Skarshewski

    Originally:

    mkdir {directory}
    cd {directory}
    compaction.py ../{RefSeq_directory}/
    cd ..

    :param RefSeq_directory: the full path as a string to a RefSeq db directory
    """
    stored = RefSeq_directory.split('-')[-1]
    try:
        os.mkdir("../" + db_version + stored)
    except OSError:
        print("Output directory exists")
        pass
    os.chdir("../" + db_version + stored)
    n = 0;

    for root, dirs, files in os.walk(RefSeq_directory):
        for filename in files:
            filepath = root + '/' + filename

            print("Compacting " + filename)

            if filepath[-3:] == '.gz':
                readfile = gzip.open(filepath, 'rt', 9)
                writefile = gzip.open(db_version + filename, 'wt', 9)
                    # writefile = open(filename + '.txt','wb')
                
                skipflag = False
                line = readfile.readline()
                while True:
                    #print(line)
                    
                    if not (line):
                        break
                    if line[0:5] == 'LOCUS':
                        #print(line)
                        identifier = ''
                        n += 1
                        j = 0
                        templine = (line[5:]).lstrip()
                        for i in templine:
                            if i == ' ':
                                identifier = templine[0:j]
                                break
                            j += 1
                        if not (identifier[0:2] == 'XM' or \
                                            identifier[0:2] == 'NM'):
                            #print(identifier)
                            while True:
                                line = readfile.readline()
                                if not line:
                                    break
                                if line.rstrip() == '//':
                                    break
                            continue
                        while True:
                            writefile.write(line)
                            line = readfile.readline()
                            if line[0:9] == 'REFERENCE' or \
                                            line[0:8] == 'FEATURES':
                                break
                        continue
                    if line[0:8] == 'FEATURES':
                        writefile.write(line)
                        line = readfile.readline()
                        while True:
                            if line[5:8] == 'CDS' or line[5:11] == 'source':
                                while True:
                                    writefile.write(line)
                                    line = readfile.readline()
                                    if not (line[0:6] == '      '):
                                        break
                                continue
                            if line[0:6] == 'ORIGIN':
                                while True:
                                    writefile.write(line)
                                    line = readfile.readline()
                                    if line == '//\n':
                                        writefile.write(line)
                                        break
                                break
                            line = readfile.readline()
                    line = readfile.readline()
                readfile.close()
                writefile.close()
    os.chdir("../")
    compacted_directory = db_version + stored
    return compacted_directory


def compile_RefSeq(compacted_dir, db_version, fn, dbuser, dbpass, dbhost, db, timeid):
    """
    Prepares a compacted RefSeq database for uPEP finding

    Provided by Adam Skarshewski

    Originally:
   
    compilation.py {compacted_dir}/

    :param compacted_dir: the full path as a string to the compacted RefSeq db 
                          directory
    """
    dataname = compacted_dir
    acc_table_name = dataname+'_acc'
    gi_table_name = dataname+'_gi'
    dbcon = MySQLdb.connect(user=dbuser, passwd=dbpass, host=dbhost, db=db)
    cursor = dbcon.cursor()
    
    #sql_dbc = ("""create database IF NOT EXISTS uPEP;""")
    
    #cursor.execute(sql_dbc)    
    
    #sql_dbtarget = ("""use uPEP;""")
    
    #cursor.execute(sql_dbtarget)
    
    
    sql_check1 = ("""drop table IF EXISTS %s;""")
    
    cursor.execute(sql_check1 % acc_table_name)
    
    cursor.execute(sql_check1 % gi_table_name)
    
    sql_tbc_ACC = ("""create table %s (accession VARCHAR(20) NOT NULL, organism VARCHAR(50) NOT NULL, position INT NOT NULL, filepath VARCHAR(1000) NOT NULL);""")
    
    cursor.execute(sql_tbc_ACC % acc_table_name)
    
    sql_tbc_GI = ("""create table %s (GI INT NOT NULL, accession VARCHAR(20) NOT NULL);""")

    cursor.execute(sql_tbc_GI % gi_table_name)
    
    n = 0
    k = 0
    #pbar = tqdm.tqdm(total=fn, desc='Total compiling progress', leave=True)
    for root, dirs, files in os.walk(compacted_dir):
        
        for filename in files:
            filepath = root + '/' + filename
            
            accessionlist = []
            GIlist = []
            if filepath[-3:] == '.gz':
                k += 1
                print("\n%i/%i Compiling (compacted) %s" % (k, fn, filename))
                readfile = gzip.open(filepath, 'rt', 9)
                skipflag = False
                pos = readfile.tell()
                while True:
                    line = readfile.readline()
                    if skipflag and not (line == '//\n'):
                        continue
                    # print line.rstrip()
                    if line[0:5] == 'LOCUS':
                        n += 1
                        temp = line[12:]
                        if (temp[:2] == 'NM') or (temp[:2] == 'XM'):
                            accession = temp[:temp.find(' ')]
                        else:
                            skipflag = True
                    elif line[0:7] == 'VERSION':
                        GI = line[line.find(':') + 1:].rstrip()
                    elif line[2:12] == 'ORGANISM  ':
                        organism = line[12:].rstrip()
                        # elif line[5:8] == 'CDS':
                        # temp = line[21:]
                        # print temp[:temp.find(' ')]
                    elif line == '//\n':
                        if not (skipflag):
                            bisect.insort(accessionlist, [accession, organism, pos,
                                                          filepath])
                            bisect.insort(GIlist, [GI, accession])
                        skipflag = False
                        pos = readfile.tell()
                    if not (line):
                        break
                readfile.close()
            # accessiondictionaryfilename = compacted_dir + '-ACCcompletecompact.dict.rna.gbff'
            # openfile = open(accessiondictionaryfilename, 'a')
            # for i in accessionlist:
                # openfile.write('%s\t%s\t%s\n' % (i[0], i[1], i[2]))
            # openfile.close()
            # GIdictionaryfilename = compacted_dir + '-GIcompletecompact.dict.rna.gbff'
            # openfile = open(GIdictionaryfilename, 'a')
            # for i in GIlist:
                # openfile.write('%s\t%s\n' % (i[0], i[1]))
            # openfile.close()
            sql_insert_ACC = ("""INSERT INTO %s (accession, organism, position, filepath) VALUES ("%s", "%s", "%s", "%s");""")
            pbar1 = tqdm.tqdm(total=len(accessionlist), desc='Accession SQL Insert', leave=True)
            for i in accessionlist:
                #print ("""Inserting Accession: %s, Organism: %s, Position: %s, Filepath: %s) VALUES (%s, %s, %s, %s);"""  % (i[0], i[1], i[2], i[3]))
                cursor.execute(sql_insert_ACC % (acc_table_name, i[0], i[1], i[2], i[3]))
                
                dbcon.commit()
                pbar1.update(1)
                
            sql_insert_GI = ("""INSERT INTO %s (GI, accession) VALUES ("%s", "%s");""")
            pbar2 = tqdm.tqdm(total=len(GIlist), desc='GI SQL Insert', leave=True)
            for i in GIlist:
                cursor.execute(sql_insert_GI % (gi_table_name, i[0], i[1]))
                dbcon.commit()
                pbar2.update(1)
            #pbar.update(1)
            
             
    print("Working with %i elements" % (n))
    cursor.close()
    dbcon.close()
    return n
def seqtocaps(seq):
    """
    Helper method for uPEP_finder
    
    Provided by Adam Skarshewski
    """
    temp = ''
    for i in seq:
        if ord(i) > 96:
            temp = temp + chr(ord(i) - 32)
        else:
            temp
    return temp


def returnjoins(string):
    """
    Helper method for uPEP_finder

    Provided by Adam Skarshewski
    """
    temp = []
    for i in range(0, len(string)):
        if string[i] == '(':
            start = i + 1
        if string[i] == ',':
            temp.append(string[start:i])
            start = i + 1
        if string[i] == ')':
            temp.append(string[start:i])
    temp2 = []
    for i in temp:
        for j in range(0, len(i)):
            if i[j] == '.':
                temp2.append([i[0:j], i[j + 2:]])
                break
    return temp2


def converttoint(string):
    """
    Helper method for uPEP_finder

    Provided by Adam Skarshewski
    """
    if string[0] == '<':
        return 0
    elif string[0] == '>':
        return int(string[1:])
    else:
        return int(string)

#@task
def uPEP_finder(codon, db_version, outpath, fn):
    os.chdir(outpath)
    print("\nBuilding uPEP database with starting codon " + codon)
    """
    Find uPEP sequences and generate a DB ready for consumption by BLAST
    
    To finish
    """
    p = 0
    specieslist = []
    proc_list = []
    for top, dirs, files in os.walk('.'):
        for directory in dirs:
            if not directory.startswith('RefSeq'):
                print(directory)
                if top == '.' and directory:
                    savefile = open(directory + codon + '.db', 'wt')
                    proc_list.append(directory + codon + '.db')
                for root, dirs, files in os.walk(directory):
                    for filename in files:
                        if filename.endswith('.gz'):
                            p += 1
                            print("(%i/%i) %s" % (p, fn, filename))
                            try:
                                takejoins = True
                                scanfile = gzip.open(root + '/' + filename, 'rt', '9')
                                while True:
                                    line = scanfile.readline()
                                    #identifier = ''
                                    if not line:
                                        break
                                    if line[0:5] == 'LOCUS':
                                        #n = 0
                                        templine = line[5:].lstrip()
                                        d = re.search('[^\s]+', templine)
                                        identifier = d.group(0)
                                        
                                        #for i in templine:
                                            #if i == ' ':
                                                #identifier = templine[0:n]
                                                #print(identifier)
                                                #break
                                            #n = n + 1
                                    elif line[2:12] == 'ORGANISM  ':
                                        if not (line[12:].rstrip() in specieslist):
                                            specieslist.append(line[12:].rstrip())
                                        if not True:  ##if not(line[12:].rstrip() in specieslist):
                                            while True:
                                                line = scanfile.readline()
                                                if not line:
                                                    break
                                                if line.rstrip() == '//':
                                                    break
                                    elif line[5:8] == 'CDS':
                                        if line[21:25] == 'join':
                                            if takejoins:
                                                more = returnjoins(line[25:].rstrip())
                                                CDSstart = converttoint(more[0][0])
                                                CDSend = converttoint(more[len(more) - 1][1])
                                        else:
                                            i = line[21:]
                                            for j in range(0, len(i)):
                                                if i[j] == '.':
                                                    CDSstart = converttoint(i[0:j])
                                                    CDSend = converttoint(i[j + 2:].rstrip())
                                                    break
                                    
                                    elif line[0:6] == 'ORIGIN':
                                        #print(identifier)
                                        temp = ''
                                        while True:
                                            line = scanfile.readline()
                                            if not line or line[0:2] == '//':
                                                break
                                            for i in line[10:].rstrip():
                                                if not i == ' ':
                                                    temp = temp + i
                                        record = [identifier, temp, CDSstart, CDSend]
                                        uORFs = []
                                        fail = False
                                        if record:
                                            CDSstart = record[2]
                                            for i in range(0, CDSstart):
                                                if not fail:
                                                    if record[1][i] == codon[0] and record[1][i + 1] == codon[1] and record[1][i + 2] == codon[2]:
                                                        j = 0
                                                        while (i + j) < CDSstart + min(100, (
                                                                CDSend - CDSstart - 5)):  # could be CDSstart + gracelength
                                                            j = j + 3
                                                            try:
                                                                if record[1][i + j] == 't' and (
                                                                            record[1][i + j + 1] == 'a' and (
                                                                            record[1][i + j + 2] == 'g' \
                                                                            or record[1][i + j + 2] == 'a') or (
                                                                        record[1][i + j + 1] == 'g' and (
                                                                    record[1][i + j + 2] == 'a'))):
                                                                    uORFs.append([record[1][i:i + j], [i + 1, i + j]])
                                                                    break
                                                            except:
                                                                print(identifier + " failed")
                                                                fail = True
                                                                break
                                        if uORFs:
                                            j = 0
                                            while True:
                                                if len(uORFs[j][0]) < 15 or len(uORFs[j][0]) > 600:
                                                    del uORFs[j]
                                                    j = j - 1
                                                j = j + 1
                                                if j == len(uORFs):
                                                    break
                                            if uORFs:
                                                
                                                for i in range(0, len(uORFs)):
                                                    uPEPloc = uORFs[i][1]
                                                    savefile.write(
                                                        '>gi|' + identifier + '_' + str(i + 1) + '|' + str(uPEPloc) + '\n')
                                                    uPEPtemp = seqtocaps(uORFs[i][0])
                                                    while len(uPEPtemp) > 80:
                                                        savefile.writelines(uPEPtemp[0:80] + '\n')
                                                        uPEPtemp = uPEPtemp[80:]
                                                    savefile.writelines(uPEPtemp + '\n')
                                                    
                                    if not line:
                                        break
                            finally:
                                scanfile.close()
                savefile.close()
    return specieslist, proc_list


def build_blast_db(upep_database_list):
    """
    Using Blast2 build a blastdb from a uPEP db
    """
    #pbar = tqdm.tqdm(total=len(upep_database_list), desc='Building blastdb', leave=True)
    for e in upep_database_list:
        os.system("formatdb -p F -i " + e)
        #pbar.update(1)

#@task
def finalise_update():   #(location=database_loc+'.refseq_version'):

    """
    Updates the RefSeq DB version store, and copies of updated databases

    :param location: the full path as a string to the RefSeq DB version store
    """
    # Update RefSeq Version
    #with open(location, 'w') as fout:
    #    fout.write(str(get_NCBI_RefSeq_release()) + '\n')
    # Create symbolic links (if needed)
    #if os.path.isfile("completed-ACCcompletecompact.dict.rna.gbff"):
    #    os.system("ln -s completed-ACCcompletecompact.dict.rna.gbff ACCcompletecompact.dict.rna.gbff")
    #    os.system("ln -s completed-GIcompletecompact.dict.rna.gbff GIcompletecompact.dict.rna.gbff")
    # Copy over *only* databases
    skip = ['fabfile.py', 'fabfile.pyc', 'README.rst', 'formatdb.log', 'config.ini']
    cur = glob.glob('*')
    for e in cur:
        if e not in skip:
            os.system("chmod 775 " + e)
            os.system("cp -r "+e+" ../")
            os.system("rm -rf "+e)
            #os.system("chmod 775 " + e)
            #os.system("mv " + e + " ../")
