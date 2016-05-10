import os, subprocess, time, random, gzip, ConfigParser
import accessions

def addnode(openfile, line):
    if not line:
        line = openfile.readline().rstrip().lstrip()
    if line:
        temp = []
        value = ''
        start = line.find('<')
        if not (start == -1):
            end = line.index('>')
            if not(end == -1):
                name = line[start+1:end]
                line = line[end+1:]
                while True:
                    start = line.find('<')
                    if not(start == -1):
                        if line[start+1] == '/':
                            end = line.index('>')
                            if not(end == -1):
                                if temp:
                                    return [name, temp]
                                else:
                                    value += line[:start]
                                    return [name, value]
                        else:
                            temp.append(addnode(openfile, line))
                    value += line
                    line = openfile.readline().rstrip().lstrip()
                    if not line:
                        break
    return []
    
                       
def parseXMLfile(XMLfile):
    openfile = open(XMLfile, 'rb')
    parselist = []
    try:        
        line = openfile.readline().rstrip().lstrip()
        while not (line == '<BlastOutput>') and line:
            line = openfile.readline().rstrip().lstrip()
        parselist = addnode(openfile, line)
    finally:
        openfile.close()
        os.remove(XMLfile)
    return parselist

def blastit(query, databases, output_dir, SEG_filter=False):
    config = ConfigParser.ConfigParser()
    config.read('config.ini')
    database_loc = config.get('Database', 'database')
    output = []
    for database in databases[0]:
        databaseoutput = []
        tempfilename =  '%s/%012x%016x' %(output_dir, int(time.time()*1000), random.randint(0,0xFFFFFFFFFFFFFFFF))
        queryfilename = tempfilename + '.fa'
        queryfile = open(queryfilename, 'wb')
        try:
            queryfile.write('>Query\n')
            queryfile.write(query)
        finally:
            queryfile.close()
        XMLfile = tempfilename + ".xml"
        if SEG_filter:
            retcode = subprocess.call(["blastall", "-p", "tblastx",  "-d", database, "-i", queryfilename, "-F", "-e", "1e-5", "-m", "7", "-o",  XMLfile ])
        else:
            retcode = subprocess.call(["blastall", "-p", "tblastx",  "-d", database, "-i", queryfilename,  "-e", "1e-5", "-m", "7", "-o",  XMLfile ])
        os.remove(queryfilename)
        if not os.path.exists(XMLfile):
            return []
        parselist = parseXMLfile(XMLfile)
        if parselist == []:
            return []
        iters = []
        for i in range(0, len(parselist[1])):
            if parselist[1][i][0] == 'BlastOutput_iterations':
                iters = parselist[1][i][1]
                break
        hits = []
        for i in range(0, len(iters[0][1])):
            if iters[0][1][i][0] == 'Iteration_hits':
                hits = iters[0][1][i][1]
        for i in range(0, len(hits)):
            for j in range(0, len(hits[i][1])):
                if hits[i][1][j][0] == 'Hit_def':
                    hit_def = hits[i][1][j][1]
                if hits[i][1][j][0] == 'Hit_hsps':
                    hsps = hits[i][1][j][1]
                    for k in hsps:
                        score = -1
                        queryrange = [-1, -1]
                        hitrange = [-1, -1]
                        queryframe = 0
                        hitframe = 0
                        hitseq = ''
                        queryseq = ''
                        match = ''
                        for m in k[1]:
                            if m[0] == 'Hsp_score':
                                score = int(m[1])
                            elif m[0] == 'Hsp_query-from':
                                queryrange[0] = int(m[1])
                            elif m[0] == 'Hsp_query-to':
                                queryrange[1] = int(m[1])                        
                            elif m[0] == 'Hsp_hit-from':
                                hitrange[0] = int(m[1])                   
                            elif m[0] == 'Hsp_hit-to':
                                hitrange[1] = int(m[1])                
                            elif m[0] == 'Hsp_query-frame':
                                queryframe = int(m[1])                
                            elif m[0] == 'Hsp_hit-frame':
                                hitframe = int(m[1])
                            elif m[0] == 'Hsp_qseq':
                                queryseq = m[1]
                            elif m[0] == 'Hsp_hseq':
                                hitseq = m[1]
                            elif m[0] == 'Hsp_midline':
                                match = m[1]
                        if hitframe == 1:
                            if databases[1]: ## Organism Check
                                defn = hit_def.split('|')[1].rsplit('_',1)[0]
                                details = accessions.getCDS(defn)
                                if details:
                                    handle = gzip.open(database_loc + details[2],'rb','9')
                                    try:
                                        while True:
                                            line = handle.readline()
                                            if line[2:10] == 'ORGANISM':
                                                organism = line[12:].rstrip()
                                                if organism == databases[1]:
                                                     databaseoutput.append([hit_def, score, queryrange, hitrange, queryframe, queryseq, hitseq, match])
                                                break
                                            if not line or line[0:2] == '//':
                                                break
                                    finally:
                                        handle.close()
                            else:
                                databaseoutput.append([hit_def, score, queryrange, hitrange, queryframe, queryseq, hitseq, match])
        i = 0
        while i + 1 < len(databaseoutput):
            j = i + 1
            endi = databaseoutput[i][0].find('_',6)
            comp = databaseoutput[i][0][3:endi]
            while j < len(databaseoutput):
                endj = databaseoutput[j][0].find('_',6)
                if (comp == databaseoutput[j][0][3:endj]) and databaseoutput[i][4] == databaseoutput[j][4]:
                    isize = eval(databaseoutput[i][0][databaseoutput[i][0].rfind('|')+1:])
                    jsize = eval(databaseoutput[j][0][databaseoutput[j][0].rfind('|')+1:])
                    if (isize[0] <= jsize[0]) and (isize[1] >= jsize[1]):
                        del databaseoutput[j]
                        continue                
                    elif (jsize[0] <= isize[0]) and (jsize[1] >= isize[1]):
                        del databaseoutput[i]
                        i -= 1
                        break            
                j += 1
            i += 1
        output += databaseoutput  
    i = 0
    while i < len(output):
        if output[i][0][3:5] == 'XM':
            del output[i]
        else:
            i = i+1         
    return output                        

#blastit('ACCTCGTTCATCTCCTTCACCTCCGAAATGATCTCGCTTTTGGGTTTACGGCCGGTCTCTTCACCTGGAGCATCAGCCGGGAAGGTCAGGGTCGCCCTGGCTCGGGCCTGTTCACATTGGGGTCAAAGGCACACATTGGGGGCTCAACCAAGGCGAGCTGCGTTCGCGGGGCCGGGTCTTTCCGCACAGGCGGAGGGCGGTGGCGGGCGCGGAGGCGTCGCGCGAGCCAGGGGGCAGCCACGGGCCGGGGGTACCTAGCGCCACCCGCTTCGCTTGCATCAGCTGCGCGCCCCATCCCGAGGAATGGTAGAGGCAGCCCCGCCCCCGGCCCGCCCCCGCCTTTCCATTGGCTGCCGCGCGGGGCGGGGAGCGGGGTCGGCTCAGTGGCCCTGAGACCCTAGCTCTGCTCTCGGTCCGCTCGCTGTCCGCTAGCCCGCTGCGATGTTGCGCGCTGCCGCCCGC',
#[['../../RefSeq/vertebrate_mammalian.db', '../../RefSeq/vertebrate_other.db'],None])

