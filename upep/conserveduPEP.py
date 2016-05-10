#!/usr/bin/env python
import os, sys, subprocess, time, random, re, MySQLdb
from upep import images, tblastx, accessions
from upep import upepsetting
### TODO ####
'''
Internal upep/refseq database mismatch weblink

Implement Min/Max of uPEP sizes (grace length partially works)

Un-hardcode several directories
sudo chown -R www-data:www-data /var/www

mysql
delete from query_info where time_id < UNIX_TIMESTAMP(DATE_SUB(NOW(), interval 7 days));
'''


def percIdent(alignstring):
    return alignstring.count('|') / float(len(alignstring))


def homologyvector(query, align, calcwindow):
    returnlist = []
    for i in range(0, len(query)):
        if query[i] == '-':
            continue  ##Only count those that are part of the query sequence (resulting vector same size as query sequence)
        else:
            if i < calcwindow / 2:
                returnlist.append(percIdent(align[:(i + calcwindow / 2)]))
            elif (i > (len(query) - calcwindow / 2)):
                returnlist.append(percIdent(align[(i - calcwindow / 2):]))
            else:
                returnlist.append(percIdent(align[(i - calcwindow / 2):(1 + i + calcwindow / 2)]))
    return returnlist


def char(value, bytes, littleendian=True):
    ans = ''
    for i in range(0, bytes):
        a = divmod(value, 256)
        ans = ans + (chr(a[1]))
        value = a[0]
    if littleendian:
        return ans
    else:
        return ans[::-1]


def ordinal(string, littleendian=True):
    value = 0
    if not (littleendian):
        temp = string[::-1]
    else:
        temp = string
    for i in range(0, len(temp)):
        value = value + pow(256, i) * ord(temp[i])
    return value


def valuetocolour(value):
    if value < 0.25:
        return (char(0, 1) + char(0, 1) + char(int(1020 * value), 1))
    elif value < 0.5:
        return (char(0, 1) + char(int(1020 * value - 255), 1) + char(int(510 - 1020 * value), 1))
    elif value < 0.75:
        return (char(int(1020 * value - 510), 1) + char(255, 1) + char(0, 1))
    else:
        return (char(255, 1) + char(int(1020 - 1020 * value), 1) + char(0, 1))


def homologytocolour(percentidentity, gradvector):
    minpos = None
    maxpos = None
    for i in range(0, 5):
        if not (gradvector[i] == None):
            if percentidentity >= gradvector[i]:
                minpos = i
            else:
                maxpos = i
                if (minpos == None):
                    return valuetocolour(maxpos * 0.25)
                else:
                    value = ((percentidentity - gradvector[minpos]) / (gradvector[maxpos] - gradvector[minpos])) * (
                        maxpos * 0.25 - minpos * 0.25) + minpos * 0.25
                    return valuetocolour(value)
    if maxpos == None:
        return valuetocolour(minpos * 0.25)
    else:
        raise ArithmeticError


def vectorResize(vector, newsize):
    if len(vector) == newsize:
        return vector
    returnlist = []
    oldsize = len(vector)
    dist = (oldsize - 1) / float(newsize - 1)
    for i in range(0, newsize):
        loc = i * dist
        frac = loc - int(loc)
        if frac < float(1) / (newsize + 1):
            returnlist.append(vector[int(loc)])
        else:
            # print (len(vector), loc, frac)
            returnlist.append((1 - frac) * vector[int(loc)] + frac * vector[int(loc) + 1])
    return returnlist


def testbmp(filename='legend.bmp'):
    width = 500
    height = 10
    openhandle = open(filename, 'wb')
    openhandle.write('BM')
    openhandle.write(char(54 + height * width * 3, 4))  # FileSize
    openhandle.write(char(0, 4))
    openhandle.write(char(54, 4))
    openhandle.write(char(40, 4))  # headersize
    openhandle.write(char(width, 4))
    openhandle.write(char(height, 4))
    openhandle.write(char(1, 2))
    openhandle.write(char(24, 2))
    openhandle.write(char(0, 4))  # Compression
    openhandle.write(char(height * width * 3, 4))
    openhandle.write(char(0, 16))
    vector = [0, 1]
    colourlist = []
    resizedvector = vectorResize(vector, width)
    for i in resizedvector:
        colourlist.append(homologytocolour(i))
    for i in range(0, height):
        for j in colourlist:
            openhandle.write(j)
    openhandle.close()


def makeheatmap(filename, vector, features, gradient, width, height):
    scale = float(width) / len(vector)
    newfeatures = []
    maximum = 1
    for i in features:
        if (i[0] > len(vector)) or (i[1] > len(vector)):
            print('<b>The parameters for %s, [%d, %d], are outside the length of the query sequence (%d nucleotides)</b></br>' % (
                i[2], i[0], i[1], len(vector)))
            continue
        temp = [int((i[0] - 1) * scale), int(i[1] * scale), i[2]]
        newset = list(range(1, maximum + 2))
        for k in newfeatures:
            if (((temp[0] > k[0] - 3) and (temp[0] < k[1] + 3)) or ((temp[1] > k[0] - 3) and (temp[1] < k[1] + 3))):
                if k[3] in newset:
                    newset.remove(k[3])
        temp.append(min(newset))
        newfeatures.append(temp)
        maximum = max(maximum, min(newset))
    PPMstring = ''
    PPMstring += 'P6\n#temp\n'
    PPMstring += '%d %d\n' % (width, height + 80 + 10 * maximum)
    PPMstring += '255\n'  # Max Colour
    
    tempstring = char(255, 1) * 30 * width * maximum
    for k in newfeatures:
        start = ((k[0] + width * 10 * (maximum - k[3])) * 3)
        end = ((k[1] + width * 10 * (maximum - k[3])) * 3)
        for i in range(0, 8):
            tempstring = tempstring[:start + width * 3 * i] + char(0, 1) * (end - start) + tempstring[end + width * 3 * i:]
    PPMstring += tempstring
    colourlist = []
    resizedvector = vectorResize(vector, width)
    for i in resizedvector:
        colourlist.append(homologytocolour(i, gradient))
    
    for i in range(0, height):
        for j in colourlist:
            PPMstring += j
    new = []
    for i in gradient:
        if i is not None:
            new.append(i)
    lower = min(new)
    upper = max(new)
    if lower > 0:
        lowertext = '<%i%%' % (int(lower * 100))
    else:
        lowertext = '0%'
    if upper < 1:
        uppertext = '>%i%%' % (int(upper * 100))
    else:
        uppertext = '100%'
    midtext = '%i%%' % (int((lower + upper) * 50))
    resizedvector = vectorResize([lower, upper], 290)
    readfile = open(upepsetting.APPS+'/legend.ppms', 'rb')
    openfile = open(filename, 'wb')
    
    openfile.write(PPMstring)
    
    openfile.write((char(255, 1) * 3 * width * 23))
    temp = images.addtext(lowertext, char(255, 1) * 3 * width * 12, width, width - 309)
    temp = images.addtext(midtext, temp, width, width - 162)
    
    openfile.write((images.addtext(uppertext, temp, width, width - 19)))
    openfile.write((char(255, 1) * 3 * width * 2))
    for i in range(0, 6):
        readfile.seek(3 * (800 * (38 + i) + 494))
        temp = (char(255, 1) * 3 * (width - 306)) + readfile.read(3 * 306)
        openfile.write(temp)
    line = char(255, 1) * 3 * (width - 306)
    for i in resizedvector:
        line += (homologytocolour(i, gradient))
    line += char(255, 1) * 3 * 16
    openfile.write((line * 16))
    for i in range(0, 21):
        readfile.seek(3 * (800 * (59 + i) + 494))
        temp = (char(255, 1) * 3 * (width - 306)) + readfile.read(3 * 306)
        openfile.write(temp)
    temp = readfile.read()
    openfile.write(temp)
    openfile.close()
    readfile.close()


def seqcheck(seq, allowcharacters):
    newseq = []
    for i in seq:
        if i in allowcharacters:
            newseq.append(i)
        if not (i in allowcharacters):
            newseq.append('N')
    checkedseq = ''.join(newseq)
    return checkedseq


def cleanDataDir():
    data_loc = upepsetting.DATA

    for root, dirs, files in os.walk(data_loc):
        for filename in files:
            filepath = root + "/" + filename
            fstat = os.stat(filepath)
            # Check for files 1 week old or more and remove
            if ((time.time() - fstat.st_mtime) > (60 * 60 * 24 * 7)):
                os.remove(filepath)

def split_sequence(input_stuff):
    n = "."*80+"?"
    wrapped = re.findall(n, input_stuff)
    return wrapped

def weboutput(unique_file_id, tid, cod, db_version, seqqu, database_n, minm,
              maxm, gracem, heatms, KaKs_t, window_t, grad, refhm):
    upep_loc = upepsetting.APP_ROOT
    apps_loc = upepsetting.APPS
    lagan_loc = upepsetting.LAGAN
    
    upeplib_loc = upepsetting.UPEPLIB
    data_loc = upepsetting.DATA
    database = upepsetting.UPEPHELPER_DATABASE

    dbuser = upepsetting.DATABASES['default']['USER']
    dbpass = str(upepsetting.DATABASES['default']['PASSWORD'])
    dbhost = upepsetting.DATABASES['default']['HOST']
    daba = upepsetting.DATABASES['default']['DB']

    tempfilename = os.path.join(data_loc, unique_file_id)

    codon = cod.replace("U", "t").replace("A", "a").replace("G", "g").replace("C", "c")

    mammaldb = os.path.join(database, db_version + upepsetting.MAMMAL + codon + '.db')
    non_mammalvertdb = os.path.join(database, db_version + upepsetting.NON_MAMMALIAN_VERTEBRATES + codon + '.db')
    invertdb = os.path.join(database, db_version + upepsetting.INVERTEBRATES + codon + '.db')
    plantdb = os.path.join(database, db_version + upepsetting.PLANTS + codon + '.db')
    fungidb = os.path.join(database, db_version + upepsetting.FUNGI + codon + '.db')
    completedb = os.path.join(database + db_version + upepsetting.COMPLETE + codon + '.db')
    dbcon = MySQLdb.connect(user=dbuser, passwd=dbpass, host=dbhost, db=daba)
    cursor = dbcon.cursor()
    codon_aminoa = {'AUG': 'M', 'AUA': 'I', 'AUC': 'I', 'ACG': 'T', 'AUU': 'I', 'AAG': 'K', 'AGG': 'R', 'CUG': 'L',
                    'UUG': 'L', 'GUG': 'V'}
    sql_insert_query_info = (
        """INSERT INTO query_info (time_id, unique_id , seqquery, database_name, db_version, protein_name, tempfilename, """
        """initial_flag, trivial_flag, blastmatches_flag, CDSKaKs_flag, uPEPKaKs_flag, heatmap_flag, window, Refheatmap, """
        """starting_aa) """
        """VALUES ("%f", "%s", 0, 0, 0, 0, "%s", 0, 0, 0, 0, 0, 0, 0, 0, "%s");""")
    cursor.execute(sql_insert_query_info % (tid, unique_file_id, tempfilename, codon_aminoa[cod]))
    dbcon.commit()
    cursor.close()
    dbcon.close()
    CDSfeature = []
    calcKaKs = True
    dictionaries = dict([('Human', [[mammaldb], 'Homo sapiens']),
                         ('Mouse', [[mammaldb], 'Mus musculus']),
                         ('Mammals', [[mammaldb], None]),
                         ('Non-mammalian vertebrates', [[non_mammalvertdb], None]),
                         ('All vertebrates', [[mammaldb, non_mammalvertdb], None]),
                         ('Invertebrates', [[invertdb], None]),
                         ('Plants', [[plantdb], None]),
                         ('Fungi', [[fungidb], None]),
                         ('Complete', [[completedb], None])])
    dbstuff = db_version + upepsetting.DBRE[database_n]

    if seqqu:
        seqquery = seqqu.upper().replace('\n', '').replace('\r', '').replace(' ', '')
        if (seqquery[:3] == "NM_") or (seqquery[:3] == "XM_") or (seqquery[:3] == "GI:"):
            calcKaKs = True
            mRNAparams = [0, 0, 20]
            if minm:
                mRNAparams[0] = int(minm)
            if maxm:
                mRNAparams[1] = int(maxm)
            if gracem:
                mRNAparams[2] = int(gracem)
            details = accessions.getmRNA(seqquery, dbstuff, codon, mRNAparams[0] * 3, mRNAparams[1] * 3,
                                         mRNAparams[2])
            if details:
                alignquery = details[0]
                blastquery = alignquery[
                             :details[1][0] + mRNAparams[2]]  # send only the 5'UTR + grace of alignquery to tblastx.py
                CDSfeature = details[1] + ['CDS']

                dbcon = MySQLdb.connect(user=dbuser, passwd=dbpass, host=dbhost, db=daba)
                cursor = dbcon.cursor()
                sql_update_query_info = (
                    """UPDATE query_info SET seqquery = "%s", database_name = "%s", db_version = "%s", protein_name = "%s", initial_flag = 2 WHERE unique_id = "%s";""")
                cursor.execute(sql_update_query_info % (seqquery, database_n, db_version, details[3], unique_file_id))
                dbcon.commit()
                cursor.close()
                dbcon.close()
            else:
                dbcon = MySQLdb.connect(user=dbuser, passwd=dbpass, host=dbhost, db=daba)
                cursor = dbcon.cursor()
                sql_update_query_info = (
                    """UPDATE query_info SET seqquery = "%s", database_name = "%s", db_version = "%s", initial_flag = 1 WHERE unique_id = "%s";""")
                cursor.execute(sql_update_query_info % (seqquery, database_n, db_version, unique_file_id))
                dbcon.commit()
                cursor.close()
                dbcon.close()
                return 0
        else:
            print('<b>QUERY: User Entered</b>')
            queryname = 'Query'
            dbcon = MySQLdb.connect(user=dbuser, passwd=dbpass, host=dbhost, db=daba)
            cursor = dbcon.cursor()
            sql_update_query_info = (
                """UPDATE query_info SET seqquery = "Query", database_name = "%s", db_version = "%s", initial_flag = 3 WHERE unique_id = "%s";""")
            cursor.execute(sql_update_query_info % (database_n, db_version, unique_file_id))
            dbcon.commit()
            cursor.close()
            dbcon.close()

            alignquery = seqquery
            blastquery = alignquery
        ## TODO Un-hardcore these directories to make this code less terrible....
        if seqcheck(alignquery, ('A', 'T', 'C', 'G')) is not None:
            changedseq = seqcheck(alignquery, ('A', 'T', 'C', 'G'))
            queryname = seqquery

            if "heatmapsize" not in form:
                dbcon = MySQLdb.connect(user=dbuser, passwd=dbpass, host=dbhost, db=daba)
                cursor = dbcon.cursor()
                sql_update_query_info = ("""UPDATE query_info SET heatmap_flag = 1 WHERE unique_id = "%s";""")
                cursor.execute(sql_update_query_info % (unique_file_id))
                dbcon.commit()
                cursor.close()
                dbcon.close()

                return 0
            heatmapsize = int(heatms)
            if (heatmapsize < 400) or (heatmapsize > 10000):
                dbcon = MySQLdb.connect(user=dbuser, passwd=dbpass, host=dbhost, db=daba)
                cursor = dbcon.cursor()
                sql_update_query_info = ("""UPDATE query_info SET heatmap_flag = 2 WHERE unique_id = "%s";""")
                cursor.execute(sql_update_query_info % (unique_file_id))
                dbcon.commit()
                cursor.close()
                dbcon.close()

                return 0
            blastmatches = tblastx.blastit(blastquery, dictionaries[database_n], data_loc)
            # print blastmatches
            if len(blastmatches) == 0:
                dbcon = MySQLdb.connect(user=dbuser, passwd=dbpass, host=dbhost, db=daba)
                cursor = dbcon.cursor()
                sql_update_query_info = ("""UPDATE query_info SET blastmatches_flag = 1 WHERE unique_id = "%s";""")
                cursor.execute(sql_update_query_info % (unique_file_id))
                dbcon.commit()
                cursor.close()
                dbcon.close()

                return 0
            elif len(blastmatches) == 1:

                dbcon = MySQLdb.connect(user=dbuser, passwd=dbpass, host=dbhost, db=daba)
                cursor = dbcon.cursor()
                sql_update_query_info = ("""UPDATE query_info SET blastmatches_flag = 2 WHERE unique_id = "%s";""")

                cursor.execute(sql_update_query_info % (unique_file_id))
                dbcon.commit()
                sql_insert_query_matchnumber = (
                    """INSERT INTO query_matchnumber (time_id, unique_id , match_number) VALUES (current_timestamp, "%s", "%i");""")
                cursor.execute(sql_insert_query_matchnumber % (unique_file_id, len(blastmatches)))
                dbcon.commit()
                cursor.close()
                dbcon.close()
            else:

                dbcon = MySQLdb.connect(user=dbuser, passwd=dbpass, host=dbhost, db=daba)
                cursor = dbcon.cursor()
                sql_update_query_info = ("""UPDATE query_info SET blastmatches_flag = 3 WHERE unique_id = "%s";""")
                cursor.execute(sql_update_query_info % (unique_file_id))
                dbcon.commit()
                sql_insert_query_matchnumber = (
                    """INSERT INTO query_matchnumber (time_id, unique_id , match_number) VALUES (current_timestamp, "%s", "%i");""")
                cursor.execute(sql_insert_query_matchnumber % (unique_file_id, len(blastmatches)))
                dbcon.commit()
                cursor.close()
                dbcon.close()
            n = 0
            for blastmatch in blastmatches:

                hitdef = blastmatch[0][3:blastmatch[0].find('_', 6)]  ##Extract RefSeq Accession
                uPEPloc = eval(blastmatch[0][blastmatch[0].find('['):])  ##Extract uPEP location
                if hitdef == seqquery:
                    dbcon = MySQLdb.connect(user=dbuser, passwd=dbpass, host=dbhost, db=daba)
                    cursor = dbcon.cursor()
                    sql_update_query_info = ("""UPDATE query_info SET trivial_flag = 1 WHERE unique_id = "%s";""")
                    cursor.execute(sql_update_query_info % (unique_file_id))
                    dbcon.commit()
                    sql_insert_query_trivial = (
                        """INSERT INTO query_trivial (time_id, unique_id, hitdef) VALUES (current_timestamp, "%s", "%s");""")
                    cursor.execute(sql_insert_query_trivial % (unique_file_id, hitdef))
                    dbcon.commit()
                    cursor.close()
                    dbcon.close()
                    continue
                details = accessions.getmRNA(hitdef, dbstuff, codon, 60, 300, 20)
                if details:
                    alignref = details[0]

                else:

                    sys.stderr.write('Error: Internal uPEP database/Refseq database mismatch. ' + hitdef)
                    sys.stderr.flush()
                    continue

                dbcon = MySQLdb.connect(user=dbuser, passwd=dbpass, host=dbhost, db=daba)
                cursor = dbcon.cursor()
                sql_insert_query_non_trvial = (
                    """INSERT INTO query_non_trivial (time_id, unique_id, protein_name, input_sequence, target_organism, input_hitdef,"""
                    """input_starting_position, input_ending_position, alignment, target_sequence, target_hitdef, target_starting_position,"""
                    """target_ending_position, upep_kaks, cds_kaks) VALUES (current_timestamp, "%s", "%s", "%s", "%s", "%s", "%i", "%i", "%s", "%s", "%s", "%i", "%i", 0, 0);""")
                cursor.execute(sql_insert_query_non_trvial % (
                    unique_file_id, details[3], blastmatch[5], details[4], queryname, blastmatch[2][0],
                    blastmatch[2][1],
                    blastmatch[7], blastmatch[6], hitdef, uPEPloc[0] + blastmatch[3][0] - 1,
                    uPEPloc[0] + blastmatch[3][1] - 1))

                dbcon.commit()
                cursor.close()
                dbcon.close()
                ### GET DETAILS FOR KaKs -- yn00 ####
                if KaKs_t == 1 and calcKaKs:
                    dbcon = MySQLdb.connect(user=dbuser, passwd=dbpass, host=dbhost, db=daba)
                    cursor = dbcon.cursor()
                    sql_update_query_info = (
                    """UPDATE query_info SET CDSKaKs_flag = 1, uPEPKaKs_flag = 1 WHERE unique_id = "%s";""")
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
                        yn00file = open(tempfilename + 'yn00', 'wt')

                        try:
                            yn00file.write('  2  %i\r\n' % (len(uPEPpair[0])))
                            yn00file.write('query\r\n' + uPEPpair[0] + '\r\n' + 'ref\r\n' + uPEPpair[1] + '\r\n')
                        finally:
                            yn00file.close()
                        retcode = subprocess.call(
                            [os.path.join(apps_loc, "yn00"), tempfilename + 'yn00', tempfilename + 'yn00' + 'uPEP'],
                            stdout=sys.stderr, stderr=sys.stderr)
                        os.remove(tempfilename + 'yn00')
                        kaksfile = open(tempfilename + 'yn00' + 'uPEP', 'rt')
                        try:
                            try:
                                for i in range(0, 6):
                                    line = kaksfile.readline().rstrip()
                                temp = line.split(' ', 1)
                                kaksuPEP = [temp[0]]
                                kaksuPEP += temp[1][1:-1].split(' ')
                            except:
                                sys.stderr.write("Line = {0}".format([temp]))
                                sys.stderr.flush()
                        finally:
                            kaksfile.close()
                        os.remove(tempfilename + 'yn00' + 'uPEP')
                        KaKs_success = True
                        try:
                            Ka = abs(float(kaksuPEP[1]))
                            Ks = abs(float(kaksuPEP[2]))
                            KaKs = abs(float(kaksuPEP[0]))
                        except:
                            KaKs_success = False
                        if not KaKs_success:
                            n += 1

                            kkr = 'Unable to estimate Ka/Ks ratio of uPEP.'
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

                            kkr = 'Estimated uPEP Ka/Ks ratio: N/A (Ka: %.4f, Ks: %.4f)' % (
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

                            kkr = 'Estimated uPEP Ka/Ks ratio: %.4f (Ka: %.4f, Ks: %.4f)' % (
                                KaKs, Ka, Ks)
                            dbcon = mysql.connector.connect(user=dbuser, passwd=dbpass, host=dbhost,
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
                        kkr = 'Unable to estimate Ka/Ks ratio of uPEP: Query and Reference uPEPs of unequal size (%i vs. %i).' % (
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
                        if len(CDSpair[0]) == len(CDSpair[1]):
                            yn00file = open(tempfilename + 'yn00', 'wt')
                            try:
                                yn00file.write('  2  %i\r\n' % (len(CDSpair[0])))
                                yn00file.write('query\r\n' + CDSpair[0] + '\r\n' + 'ref\r\n' + CDSpair[1] + '\r\n')
                            finally:
                                yn00file.close()
                            retcode = subprocess.call(
                                [os.path.join(apps_loc, "yn00"), tempfilename + 'yn00', tempfilename + 'yn00' + 'CDS'])
                            os.remove(tempfilename + 'yn00')
                            try:
                                kaksfile = open(tempfilename + 'yn00' + 'CDS', 'rt')
                                for i in range(0, 6):
                                    line = kaksfile.readline().rstrip()
                                temp = line.split(' ', 1)
                                kaksCDS = [temp[0]]
                                kaksCDS += temp[1][1:-1].split(' ')
                            finally:
                                kaksfile.close()
                            os.remove(tempfilename + 'yn00' + 'CDS')
                            KaKs_success = True
                            try:
                                Ka = abs(float(kaksCDS[1]))
                                Ks = abs(float(kaksCDS[2]))
                                KaKs = abs(float(kaksCDS[0]))
                            except:
                                KaKs_success = False
                            if not KaKs_success:
                                n += 1
                                kkr = 'Unable to estimate Ka/Ks ratio of CDS.'
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

                                kkr = 'Estimated CDS Ka/Ks ratio:&nbsp;&nbsp;N/A (Ka: %.4f, Ks: %.4f)' % (
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

                                kkr = 'Estimated CDS Ka/Ks ratio:&nbsp;&nbsp;%.4f (Ka: %.4f, Ks: %.4f)' % (
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

                            kkr = 'Unable to estimate Ka/Ks ratio of CDS: Query and reference coding sequences of unequal size (%i vs. %i).' % (
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

                        kkr = 'Unable to estimate Ka/Ks ratio of CDS: Unable to define CDS from user entered sequence.'
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
                reffeatures.append(details[1] + ['CDS'])
                for j in details[2]:
                    reffeatures.append(j + ['uPEP'])
                if CDSfeature:
                    features.append(CDSfeature)
                features.append(blastmatch[2] + ['uPEP'])
                ##if True: #i[3][0] == 1:
                ##    retcode = subprocess.call(["./yn00", 'ptp4a1upep.nuc', 'ptpuPEP'])
                alignfile1 = open(tempfilename + hitdef + "_1", 'wt')
                alignfile2 = open(tempfilename + hitdef + "_2", 'wt')
                try:
                    alignfile1.write(">Query\n" + changedseq + "\n")
                    alignfile2.write(">Reference\n" + alignref + "\n")
                finally:
                    alignfile1.close()
                    alignfile2.close()
                os.environ["LAGAN_DIR"] = lagan_loc
                retcode = subprocess.call(
                    ["perl", os.path.join(lagan_loc, "lagan.pl"), tempfilename + hitdef + "_1", tempfilename + hitdef + "_2", "-out",
                     tempfilename + hitdef + ".mfa", "-mfa"])
                retcode = subprocess.call(
                    ["perl", os.path.join(lagan_loc, "mf_to_align.pl"), "-f", tempfilename + hitdef + ".mfa", "-out",
                     tempfilename + hitdef + "aligneduf"])
                aligneduffile = open(tempfilename + hitdef + "aligneduf", 'rt')
                try:
                    query = "".join(aligneduffile.readline().rstrip('\r\n'))
                    match = "".join(aligneduffile.readline().rstrip('\r\n'))
                    ref = "".join(aligneduffile.readline().rstrip('\r\n'))
                    alignedfile = open(tempfilename + hitdef + "aligned", 'wt')
                    query2 = split_sequence(query)
                    match2 = split_sequence(match)
                    ref2 = split_sequence(ref)
                    for i1, i2, i3 in zip(query2, match2, ref2):
                        alignedfile.write('%s\n%s\n%s\n\n' % (i1, i2, i3))
                    alignedfile.close()
                    window = window_t
                    dbcon = MySQLdb.connect(user=dbuser, passwd=dbpass, host=dbhost, db=daba)
                    cursor = dbcon.cursor()
                    sql_update_query_info = ("""UPDATE query_info SET window = "%i" WHERE unique_id = "%s";""")
                    cursor.execute(sql_update_query_info % (int(window), unique_file_id))
                    dbcon.commit()
                    cursor.close()
                    dbcon.close()
                    windowfilename = tempfilename + hitdef + 'w' + window

                    if refhm == 1:
                        makeheatmap(windowfilename + 'r.ppm', homologyvector(ref, match, int(window)), reffeatures,
                                    grad, heatmapsize, 25)
                    makeheatmap(windowfilename + '.ppm', homologyvector(query, match, int(window)), features, grad,
                                heatmapsize, 25)
                finally:
                    aligneduffile.close()
                if refhm == 1:
                    retcode = subprocess.call([os.path.join(apps_loc, "png.py"), "-c", "9", windowfilename + 'r.ppm'])
                    os.remove(windowfilename + 'r.ppm')
                retcode = subprocess.call([os.path.join(apps_loc, "png.py"), "-c", "9", windowfilename + '.ppm'])
                os.remove(windowfilename + '.ppm')

                if refhm == 1:
                    dbcon = MySQLdb.connect(user=dbuser, passwd=dbpass, host=dbhost, database='uPEP_result')
                    cursor = dbcon.cursor()
                    sql_update_query_info = ("""UPDATE query_info SET Refheatmap = 1 WHERE unique_id = "%s";""")
                    cursor.execute(sql_update_query_info % (unique_file_id))
                    dbcon.commit()
                    cursor.close()
                    dbcon.close()

        else:
            errorstuff = 'align error'


def parsing_initial(unid):
    dbuser = upepsetting.DATABASES['default']['USER']
    dbpass = str(upepsetting.DATABASES['default']['PASSWORD'])
    dbhost = upepsetting.DATABASES['default']['HOST']
    daba = upepsetting.DATABASES['default']['DB']
    dbcon = MySQLdb.connect(user=dbuser, passwd=dbpass, host=dbhost, db=daba)
    cursor = dbcon.cursor()
    sql_retrieve_query = """select * from query_info where unique_id = "%s";"""
    cursor.execute(sql_retrieve_query % unid)
    for (
            time_id, unique_id, seqquery, database_name, db_version, protein_name, tempfilename, initial_flag,
            trivial_flag,
            blastmatches_flag, CDSKaKs_flag, uPEPKaKs_flag, heatmap_flag, window, Refheatmap, starting_aa) in cursor:
        seqq = seqquery
        dbn = database_name
        prn = protein_name
        dbv = db_version
        tfn = tempfilename
        ini_flag = initial_flag
        trv_flag = trivial_flag
        bm_flag = blastmatches_flag
        hm_flag = heatmap_flag
        cdsk_flag = CDSKaKs_flag
        upepk_flag = uPEPKaKs_flag
        win_r = window
        rhm = Refheatmap
        saa = starting_aa

        return seqq, dbn, prn, dbv, tfn, ini_flag, trv_flag, bm_flag, hm_flag, cdsk_flag, upepk_flag, win_r, rhm, saa


def strictmode_and_filter_trigger(unique_file_id, strict, filter_org, starting_amino):

    if strict == 1:
        if filter_org == 'None':
            sql_non_trivial_retrieve = """select * from query_non_trivial where (unique_id = '%s' and target_sequence like "%s");""" % (
            unique_file_id, (starting_amino+'%'))
            match_string = 'Non-strict mode: %i hits found.'
            return sql_non_trivial_retrieve, match_string

        else:
            sql_non_trivial_retrieve = """select * from query_non_trivial where (unique_id = '%s' and target_organism = '%s' and target_sequence like "%s");""" % (
            unique_file_id, filter_org, (starting_amino+'%'))
            match_string = 'Strict mode: %i filtered hits found.'
            return sql_non_trivial_retrieve, match_string
    else:
        if filter_org == 'None':
            sql_non_trivial_retrieve = """select * from query_non_trivial where unique_id = "%s";""" % unique_file_id
            match_string = 'Non-strict mode: %i hits found.'
            return sql_non_trivial_retrieve, match_string
        else:
            sql_non_trivial_retrieve = """select * from query_non_trivial where (unique_id = '%s' and target_organism = "%s");""" % (
            unique_file_id, filter_org)
            match_string = 'Strict mode: %i filtered hits found.'
            return sql_non_trivial_retrieve, match_string


def weboutput2(unique_file_id, filter_org, strict):
    dbuser = upepsetting.DATABASES['default']['USER']
    dbpass = str(upepsetting.DATABASES['default']['PASSWORD'])
    dbhost = upepsetting.DATABASES['default']['HOST']
    daba = upepsetting.DATABASES['default']['DB']
    dbcon = MySQLdb.connect(user=dbuser, passwd=dbpass, host=dbhost, db=daba)
    cursor = dbcon.cursor()

    seqq, dbn, prn, dbv, tfn, ini_flag, trv_flag, bm_flag, hm_flag, cdsk_flag, \
    upepk_flag, win_r, rhm, saa = parsing_initial(unique_file_id)

    print('<b>Using RefSeq Release: <b>' + dbv + '</b><br>')
    print('<b>Database:' + dbn + '</b><br>')
    print('<b>Job ID: </b>' + unique_file_id + '<br>')
    print()
    if ini_flag == 1:
        # print 'Error: Unable to find query mRNA sequence in selected Refseq database.' + seqq + ' not found \n \
        print('No results were recorded. \n \
                Click <a href="javascript:history.go(-1)"> here </a> to return to the input page.')
        print('</body></html>')
        cursor.close()
        dbcon.close()
        return 0
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
            return 0

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
            print('<b> %i hit(s) with the same transcript as query sequence </b><br><br>' % cursor.rowcount)
        sql_non_trivial_retrieve, match_string = strictmode_and_filter_trigger(unique_file_id, strict, filter_org, saa)

        cursor.execute(sql_non_trivial_retrieve)
        if not cursor.rowcount:
            print('<br>')
            print('Error: Internal uPEP database/Refseq database mismatch. Query not found \n \
            This error has been logged and sent to the webmaster. See help for further details.')
            cursor.close()
            dbcon.close()
            return 0

        else:

            result = cursor.fetchall()
            if cursor.rowcount > 0:
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
                            print('<span style="font-size:80%"><b>' + i[3] + '</b></span><br>')
                    if cdsk_flag == 1:
                        sql_cdsk_retrieve = """select * from query_CDS_KaKs where (unique_id = "%s" and hitdef = "%s" and cds_kaks_trigger = "%i");"""
                        cursor.execute(sql_cdsk_retrieve % (unique_file_id, row[10], row[14]))
                        cdsk_result = cursor.fetchall()
                        for i in cdsk_result:
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

    else:
        print('<br>')
        print('No results was recorded. Click <a href="javascript:history.go(-1)"> here </a> to return to the input page.')
        print('</body></html>')
        cursor.close()
        dbcon.close()
        return 0

        # print 'Error: The width of the heatmap was not specified.  \
        # Click <a href="javascript:history.go(-1)"> here </a> to return to the input page.'
        # print '</body></html>'
        # cursor.close()
        # dbcon.close()
        # return 0

        # elif hm_flag == 2:
        # print '<br>'
        # print 'Error: The heapmap width must be between 400 and 10000 pixels.  \
        # Click <a href="javascript:history.go(-1)"> here </a> to return to the input page.'
        # print '</body></html>'
        # cursor.close()
        # dbcon.close()
        # return 0
