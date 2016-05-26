import gzip
import sys
import os
from upep import upepsetting
import MySQLdb
def seqtocaps(seq):
    temp = ''
    for i in seq:
        if ord(i) > 96:
            temp = temp + chr(ord(i)-32)
        else: temp
    return temp

def ACCloc(accession, dbstuff):
    database_loc = upepsetting.UPEPHELPER_DATABASE
    dbuser = upepsetting.DATABASES['default']['USER']
    dbpass = str(upepsetting.DATABASES['default']['PASSWORD'])
    dbhost = upepsetting.DATABASES['default']['HOST']
    daba = upepsetting.DATABASES['default']['DB']

    dbconnect = MySQLdb.connect(user=dbuser, passwd=dbpass, host=dbhost, db=daba)
    cursor = dbconnect.cursor()
    sql_retrieve = """select * from """+dbstuff+"""_acc where accession = %s"""
    cursor.execute(sql_retrieve, (accession,))
    for (accession, organism, position, filepath) in cursor:
        acc = accession
        org = organism
        pos = position
        filep = filepath
        cursor.close()
        dbconnect.close()
        return acc, org, pos, filep 

def GIloc(accession, dbstuff):
    database_loc = upepsetting.UPEPHELPER_DATABASE
    dbuser = upepsetting.DATABASES['default']['USER']
    dbpass = str(upepsetting.DATABASES['default']['PASSWORD'])
    dbhost = upepsetting.DATABASES['default']['HOST']
    daba = upepsetting.DATABASES['default']['DB']

    dbconnect = MySQLdb.connect(user=dbuser, passwd=dbpass, host=dbhost, db=daba)
    cursor = dbconnect.cursor()
    sql_retrieve = """select * from """+dbstuff+"""_gi where GI = %s"""
    cursor.execute(sql_retrieve, (accession,))
    for (GI, accession) in cursor:
        GI = GI
        acc = accession
        cursor.close()
        dbconnect.close()
        acc2, org2, pos, filep = ACCloc(acc)
        return acc2, org2, pos, filep

   

def getCDS(accession, dbstuff):
    sys.stderr.write("Debug 2: "+str(accession)+'\n')
    sys.stderr.flush()
    if len(accession) > 3:
        if accession[:3] == 'GI:':
            sys.stderr.write("Debug 3: GI \n")
            sys.stderr.flush()
            return GIloc(accession[3:], dbstuff)
        elif (accession[:3] == 'NM_' or accession[:3] == 'XM_'):
            sys.stderr.write("Debug 3: AC \n")
            sys.stderr.flush()
            return ACCloc(accession, dbstuff)
        else:
            return None
    else:
        return None

def getmRNA(accession, dbstuff, codon, minsize, maxsize, gracelength):
    database_loc = upepsetting.UPEPHELPER_DATABASE
    takejoins = True
    details = getCDS(accession, dbstuff)
    sys.stderr.write("Debug 1: "+str(details)+'\n')
    sys.stderr.flush()
    if details:
        handle = gzip.open(database_loc + details[3],'rt','9')
        try:
            print('mrna')
            handle.seek(int(details[2]))
            while True:
                line = handle.readline()
                if line[0:10] == 'DEFINITION':
                    definition = line[12:].rstrip()
                    line = handle.readline()
                    while True:
                        if line[0] == ' ':
                            definition = definition + ' ' + line[12:].rstrip()
                            line = handle.readline()
                        else:
                            break                        
                if line[5:8] == 'CDS':   
                    if line[21:25] == 'join':
                        if takejoins:
                            more = returnjoins(line[25:].rstrip())
                            CDSstart = int(more[0][0])
                            CDSend = int(more[len(more)-1][1])
                    else:
                        i = line[21:]
                        for j in range(0,len(i)):
                            if i[j] == '.':
                                CDSstart = int(i[0:j])
                                CDSend = int(i[j+2:].rstrip())
                                break                 
                elif line[0:6] == 'ORIGIN':
                    temp = ''
                    while True:
                        line = handle.readline()
                        if not line or line[0:2] == '//':
                            break
                        for i in line[10:].rstrip():
                            if not i == ' ':
                                temp = temp + i                         
                    uORFs = []
                    if maxsize:
                        for i in range(0, CDSstart):
                            if temp[i] == codon[0] and temp[i+1] == codon[1] and temp[i+2] == codon[2]:
                                j = 0
                                while (i + j) < CDSstart + gracelength:
                                    j = j + 3
                                    if temp[i+j] == 't' and (temp[i+j+1] == 'a' and (temp[i+j+2] == 'g' \
                                        or temp[i+j+2] == 'a') or (temp[i+j+1] == 'g' and (temp[i+j+2] == 'a'))):
                                        uORFs.append([i,i+j])
                                        break    
                    if uORFs:
                        j = 0
                        while True:
                            if uORFs[j][1]-uORFs[j][0] < minsize or uORFs[j][1]-uORFs[j][0] > maxsize:
                                del uORFs[j]
                                j = j - 1
                            j = j + 1
                            if j == len(uORFs):
                                break
                    temp = seqtocaps(temp)               
                    return [temp, [CDSstart, CDSend], uORFs, definition, details[1]]       
        finally:
            handle.close()                           
    return None
