#!/usr/bin/env python
import os, sys, subprocess, time, random, cgi
from upeplib import images, tblastx, accessions

### TODO ####
'''
Internal upep/refseq database mismatch weblink

Implement Min/Max of uPEP sizes (grace length partially works)

Un-hardcode several directories
'''


def percIdent(alignstring):
    return alignstring.count('|')/float(len(alignstring))


def homologyvector(query, align, calcwindow):
    returnlist = []
    for i in range(0, len(query)):
        if query[i] == '-':
            continue   ##Only count those that are part of the query sequence (resulting vector same size as query sequence)
        else:
            if i < calcwindow/2:
                returnlist.append(percIdent(align[:(i+calcwindow/2)]))
            elif (i > (len(query) - calcwindow/2)):
                returnlist.append(percIdent(align[(i-calcwindow/2):]))
            else:
                returnlist.append(percIdent(align[(i-calcwindow/2):(1+i+calcwindow/2)]))
    return returnlist

def char(value, bytes, littleendian = True):
    ans = ''
    for i in range(0, bytes):
        a = divmod(value, 256)
        ans = ans + (chr(a[1]))
        value = a[0]
    if littleendian:
        return ans
    else:
        return ans[::-1]


def ordinal(string, littleendian = True):
    value = 0
    if not(littleendian):
        temp = string[::-1]
    else:
        temp = string
    for i in range(0, len(temp)):
        value = value + pow(256,i)*ord(temp[i])
    return value

def valuetocolour(value):
    if value < 0.25:
        return (char(0, 1) + char(int(1020*value), 1) + char(0, 1))
    elif value < 0.5:
        return (char(int(1020*value-255), 1) + char(255, 1) + char(0, 1))
    else:
        return (char(255, 1) + char(int(510-510*value), 1) + char(0, 1))

def homologytocolour(percentidentity, gradvector):
    minpos = None
    maxpos = None
    for i in range(0,5):
        if not(gradvector[i] == None):
            if percentidentity >= gradvector[i]:
                minpos = i
            else:
                maxpos = i
                if (minpos == None):
                    return valuetocolour(maxpos*0.25)
                else:
                    value = ((percentidentity-gradvector[minpos])/(gradvector[maxpos]-gradvector[minpos]))*(maxpos*0.25-minpos*0.25)+minpos*0.25
                    return valuetocolour(value)
    if maxpos == None:
        return valuetocolour(minpos*0.25)
    else:
        raise ArithmeticError

def vectorResize(vector, newsize):
    if len(vector) == newsize:
        return vector
    returnlist = []
    oldsize = len(vector)
    dist = (oldsize - 1)/float(newsize-1)
    for i in range(0, newsize):
        loc = i*dist
        frac =  loc - int(loc)
        if frac < float(1)/(newsize+1):
            returnlist.append(vector[int(loc)])
        else:
            #print (len(vector), loc, frac)
            returnlist.append((1-frac)*vector[int(loc)]+frac*vector[int(loc)+1])
    return returnlist

def testbmp(filename='legend.bmp'):
    width = 500
    height = 10
    openhandle = open(filename,'wb')
    openhandle.write('BM')
    openhandle.write(char(54 + height * width * 3, 4)) #FileSize
    openhandle.write(char(0, 4))
    openhandle.write(char(54, 4))
    openhandle.write(char(40, 4)) #headersize
    openhandle.write(char(width, 4))
    openhandle.write(char(height, 4))
    openhandle.write(char(1, 2))
    openhandle.write(char(24, 2))
    openhandle.write(char(0, 4)) #Compression
    openhandle.write(char(height * width * 3, 4))
    openhandle.write(char(0, 16))
    vector = [0, 1]
    colourlist = []
    resizedvector = vectorResize(vector, width)
    for i in resizedvector:
        colourlist.append(homologytocolour(i))
    for i in range(0,height):
        for j in colourlist:
            openhandle.write(j)
    openhandle.close()


def makeheatmap(filename, vector, features, gradient, width, height):
    scale = float(width)/len(vector)
    newfeatures = []
    maximum = 1
    for i in features:
        if (i[0] > len(vector)) or (i[1] > len(vector)):
            print '<b>The parameters for %s, [%d, %d], are outside the length of the query sequence (%d nucleotides)</b></br>' %(i[2], i[0], i[1], len(vector))
            continue
        temp = [int((i[0]-1)*scale), int(i[1]*scale), i[2]]
        newset = range(1, maximum+2)
        for k in newfeatures:
            if (((temp[0] > k[0]-3) and (temp[0] < k[1]+3)) or ((temp[1] > k[0]-3) and (temp[1] < k[1]+3))):
                if k[3] in newset:
                    newset.remove(k[3])
        temp.append(min(newset))
        newfeatures.append(temp)
        maximum = max(maximum, min(newset))
    PPMstring = ''
    PPMstring += 'P6\n#temp\n'
    PPMstring += '%d %d\n' % (width, height+80+10*maximum)
    PPMstring += '255\n' #Max Colour
    tempstring = char(255,1)*30*width*maximum
    for k in newfeatures:
        start = (k[0]+width*10*(maximum-k[3]))*3
        end = (k[1]+width*10*(maximum-k[3]))*3
        for i in range(0, 8):
            tempstring = tempstring[:start+width*3*i] + char(0,1)*(end-start) + tempstring[end+width*3*i:]
    PPMstring += tempstring
    colourlist = []
    resizedvector = vectorResize(vector, width)
    for i in resizedvector:
        colourlist.append(homologytocolour(i, gradient))
    for i in range(0,height):
        for j in colourlist:
            PPMstring += j
    new = []
    for i in gradient:
        if i is not None:
            new.append(i)
    lower = min(new)
    upper = max(new)
    if lower > 0:
        lowertext = '<%i%%' % (int(lower*100))
    else:
        lowertext = '0%'
    if upper < 1:
        uppertext = '>%i%%' % (int(upper*100))
    else:
        uppertext = '100%'
    midtext = '%i%%' % (int((lower+upper)*50))
    resizedvector = vectorResize([lower, upper], 290)
    readfile = open('legend.ppms', 'rb')
    openfile = open(filename, 'wb')
    openfile.write(PPMstring)
    openfile.write(char(255,1)*3*width*23)
    temp = images.addtext(lowertext, char(255,1)*3*width*12, width, width-309)
    temp = images.addtext(midtext, temp, width, width-162)
    openfile.write(images.addtext(uppertext, temp, width, width-19))
    openfile.write(char(255,1)*3*width*2)
    for i in range(0, 6):
        readfile.seek(3*(800*(38+i)+494))
        temp = char(255,1)*3*(width-306) + readfile.read(3*306)
        openfile.write(temp)
    line = char(255,1)*3*(width-306)
    for i in resizedvector:
        line += (homologytocolour(i, gradient))
    line += char(255,1)*3*16
    openfile.write(line*16)
    for i in range(0, 21):
        readfile.seek(3*(800*(59+i)+494))
        temp = char(255,1)*3*(width-306) + readfile.read(3*306)
        openfile.write(temp)
    temp = readfile.read()
    openfile.write(temp)
    openfile.close()
    readfile.close()

def seqcheck(seq, set):
    for i in seq:
        if not(i in set):
            return False
    return True

def cleanDataDir():
    for root, dirs, files in os.walk("../data"):
        for filename in files:
            filepath = root + "/" + filename
            fstat = os.stat(filepath)
            # Check for files 1 week old or more and remove
            if ((time.time() - fstat.st_mtime) > (60 * 60 * 24 * 7)):
                os.remove(filepath)


def weboutput():
    tempfilename =  '../data/%012x%016x' %(int(time.time()*1000), random.randint(0,0xFFFFFFFFFFFFFFFF))
    form = cgi.FieldStorage()
    if form.keys() == []:
        print "Refresh: 5; url=http://upep-scmb.biosci.uq.edu.au/"
        print "Content-type: text/html"
        print ## <----- VERY IMPORTANT!!!!!! Blank "print" - Separate Headers
        print 'You will be redirected to uPEPperoni in 5 seconds. \n\
        Please click <a href="../../">here</a> if not automatically redirected.'
        return 0
    print "Content-type: text/html"
    print ## <----- VERY IMPORTANT!!!!!! Blank "print" - Separate Headers
    print '<html><head>'
    print '<title>'
    print 'uPEPperoni - Results'
    print '</title>'
    print '</head><body>'
    with open("../RefSeq/.refseq_version") as refseq_rel_in:
        rsr = refseq_rel_in.readlines()[0].strip()
    print '<b>Using RefSeq Release: <b>'+rsr+'</b><br>'
    cleanDataDir()
    CDSfeature = []
    calcKaKs = True
    if form.has_key("seqquery"):
        seqquery = form["seqquery"].value.upper().replace('\n','').replace('\r','').replace(' ','')
        if (seqquery[:3] == "NM_") or (seqquery[:3] == "XM_") or (seqquery[:3] == "GI:"):
            calcKaKs = True
            mRNAparams = [0, 0, 20]
            if form.has_key("min"):
                mRNAparams[0] = int(form["min"].value)
            if form.has_key("max"):
                mRNAparams[1] = int(form["max"].value)
            if form.has_key("grace"):
                mRNAparams[2] = int(form["grace"].value)
            details = accessions.getmRNA(seqquery, mRNAparams[0]*3, mRNAparams[1]*3, mRNAparams[2])
            if details:
                alignquery = details[0]
                blastquery = alignquery[:details[1][0]+mRNAparams[2]] #send only the 5'UTR + grace of alignquery to tblastx.py
                CDSfeature = details[1] + ['CDS']
                print '<b>QUERY: ' + details[3] + ' (' + seqquery + ')</b>'
                queryname = seqquery
            else:
                print '<br>'
                print 'Error: Unable to find query mRNA sequence in selected Refseq database. ' + seqquery + ' not found \n \
                Click <a href="javascript:history.go(-1)"> here </a> to return to the input page.'
                print '</body></html>'
                return 0
        else:
            print '<b>QUERY: User Entered</b>'
            queryname = 'Query'
            alignquery = seqquery
            blastquery = alignquery
        ## TODO Un-hardcore these directories to make this code less terrible....
        if seqcheck(alignquery, ('A', 'T', 'C', 'G')):
            dictionaries = dict([('Human', [['../RefSeq/vertebrate_mammalian.db'],'Homo sapiens']),
                                        ('Mouse', [['../RefSeq/vertebrate_mammalian.db'],'Mus musculus']),
                                        ('Mammals', [['../RefSeq/vertebrate_mammalian.db'], None]),
                                        ('Non-mammalian vertebrates',[['../RefSeq/vertebrate_other.db'], None]),
                                        ('All vertebrates', [['../RefSeq/vertebrate_mammalian.db', '../../RefSeq/vertebrate_other.db'],None]),
                                        ('Invertebrates', [['../RefSeq/invertebrate.db'], None]),
                                        ('Plants', [['../RefSeq/plant.db'], None]),
                                        ('Fungi', [['../RefSeq/fungi.db'], None]),
                                        ('Complete',[['../RefSeq/complete.db'], None])])
            if not form.has_key("heatmapsize"):
                print '<br>'
                print 'Error: The width of the heatmap was not specified.  \
                Click <a href="javascript:history.go(-1)"> here </a> to return to the input page.'
                print '</body></html>'
                return 0
            heatmapsize = int(form["heatmapsize"].value)
            if (heatmapsize < 400) or (heatmapsize > 10000):
                print '<br>'
                print 'Error: The heapmap width must be between 400 and 10000 pixels.  \
                Click <a href="javascript:history.go(-1)"> here </a> to return to the input page.'
                print '</body></html>'
                return 0
            blastmatches = tblastx.blastit(blastquery, dictionaries[form["database"].value], "../data")
            #print blastmatches
            if len(blastmatches) == 0:
                print '<br>'
                print 'No hits were found for the given input sequence in the \'' + form["database"].value + '\' uPEP database. \n \
                Click <a href="javascript:history.go(-1)"> here </a> to return to the input page.'
                print '</body></html>'
                return 0
            elif len(blastmatches) == 1:
                print '<br><b> 1 hit found.</b><br><br>'
            else:
                print '<br><b> %i hits found.</b><br><br>' %len(blastmatches)
            for blastmatch in blastmatches:
                hitdef = blastmatch[0][3:blastmatch[0].find('_',6)] ##Extract RefSeq Accession
                uPEPloc =  eval(blastmatch[0][blastmatch[0].find('['):]) ##Extract uPEP location
                if hitdef == seqquery:
                    print '<b><u>HIT/REFERENCE:</u> Trivial hit. Same transcript as query sequence (' + hitdef + ')</b><br><br>'
                    continue
                details = accessions.getmRNA(hitdef, 60, 300, 20)
                if details:
                    alignref = details[0]
                    print '<b><u>HIT/REFERENCE:</u> ' + details[3] + ' ('+ hitdef+ ')</b><br>'
                else:
                    print '<br>'
                    print 'Error: Internal uPEP database/Refseq database mismatch. ' + hitdef + ' not found \n \
                    This error has been logged and sent to the webmaster. See help for further details.'
                    sys.stderr.write('Error: Internal uPEP database/Refseq database mismatch. ' + hitdef)
                    sys.stderr.flush()
                    continue
                print '<PRE>%s   %s, [%i, %i]\n%s\n%s   %s, [%i, %i]</PRE>' %(blastmatch[5], queryname, blastmatch[2][0], blastmatch[2][1], blastmatch[7], blastmatch[6], hitdef, uPEPloc[0] + blastmatch[3][0] - 1, uPEPloc[0] + blastmatch[3][1] - 1)
                ### GET DETAILS FOR KaKs -- yn00 ####
                if form.has_key("kaks") and calcKaKs:
                    uPEPpair = [alignquery[blastmatch[2][0]-1:blastmatch[2][1]],alignref[uPEPloc[0] + blastmatch[3][0] - 2: uPEPloc[0] + blastmatch[3][1] - 1]]
                    CDSpair = []
                    if CDSfeature:
                        CDSpair = [alignquery[CDSfeature[0]-1:CDSfeature[1]-3],alignref[details[1][0]-1:details[1][1]-3]] ### -1 to correct for python sequences // -3 to remove stop codon
                    if len(uPEPpair[0]) == len(uPEPpair[1]):
                        yn00file = open(tempfilename + 'yn00','wb')
                        try:
                            yn00file.write('  2  %i\r\n' %(len(uPEPpair[0])))
                            yn00file.write('query\r\n' + uPEPpair[0] + '\r\n' + 'ref\r\n' +uPEPpair[1] + '\r\n')
                        finally:
                            yn00file.close()
                        retcode = subprocess.call(["../apps/yn00", tempfilename + 'yn00', tempfilename + 'yn00'+ 'uPEP'], stdout=sys.stderr, stderr=sys.stderr)
                        os.remove(tempfilename + 'yn00')
                        kaksfile = open(tempfilename + 'yn00'+ 'uPEP')
                        try:
                            try:
                                for i in range(0, 6):
                                    line = kaksfile.readline().rstrip()
                                temp = line.split(' ',1)
                                kaksuPEP = [temp[0]]
                                kaksuPEP += temp[1][1:-1].split(' ')
                            except:
                                sys.stderr.write("Line = {0}".format([temp]))
                                sys.stderr.flush()
                        finally:
                            kaksfile.close()
                        os.remove(tempfilename + 'yn00'+ 'uPEP')
                        KaKs_success = True
                        try:
                            Ka = abs(float(kaksuPEP[1]))
                            Ks = abs(float(kaksuPEP[2]))
                            KaKs = abs(float(kaksuPEP[0]))
                        except:
                            KaKs_success = False
                        if (not KaKs_success):
                            print '<span style="font-size:80%"><b>Unable to estimate Ka/Ks ratio of uPEP.</b></span><br>'
                        elif (Ks == 0):
                            print '<span style="font-size:80%%"><b>Estimated uPEP Ka/Ks ratio: N/A (Ka: %.4f, Ks: %.4f)</b></span><br>' %(Ka, Ks)
                        else:
                            print '<span style="font-size:80%%"><b>Estimated uPEP Ka/Ks ratio: %.4f (Ka: %.4f, Ks: %.4f)</b></span><br>' %(KaKs, Ka, Ks)
                    else:
                        print '<span style="font-size:80%"><b>Unable to estimate Ka/Ks ratio of uPEP: Query and Reference uPEPs of unequal size (%i vs. %i).</b></span><br>' %(len(uPEPpair[0]), len(uPEPpair[1]))
                    if CDSfeature:
                        if len(CDSpair[0]) == len(CDSpair[1]):
                            yn00file = open(tempfilename + 'yn00','wb')
                            try:
                                yn00file.write('  2  %i\r\n' %(len(CDSpair[0])))
                                yn00file.write('query\r\n' + CDSpair[0] + '\r\n' + 'ref\r\n' +CDSpair[1] + '\r\n')
                            finally:
                                yn00file.close()
                            retcode = subprocess.call(["../apps/yn00", tempfilename + 'yn00', tempfilename + 'yn00'+ 'CDS'])
                            os.remove(tempfilename + 'yn00')
                            try:
                                kaksfile = open(tempfilename + 'yn00'+ 'CDS')
                                for i in range(0, 6):
                                    line = kaksfile.readline().rstrip()
                                temp = line.split(' ',1)
                                kaksCDS = [temp[0]]
                                kaksCDS += temp[1][1:-1].split(' ')
                            finally:
                                kaksfile.close()
                            os.remove(tempfilename + 'yn00'+ 'CDS')
                            KaKs_success = True
                            try:
                                Ka = abs(float(kaksCDS[1]))
                                Ks = abs(float(kaksCDS[2]))
                                KaKs = abs(float(kaksCDS[0]))
                            except:
                                KaKs_success = False
                            if (not KaKs_success):
                                print '<span style="font-size:80%"><b>Unable to estimate Ka/Ks ratio of CDS.</b></span><br>'
                            elif (Ks == 0):
                                print '<span style="font-size:80%%"><b>Estimated CDS Ka/Ks ratio:&nbsp;&nbsp;N/A (Ka: %.4f, Ks: %.4f)</b></span><br><br>' %(Ka, Ks)
                            else:
                                print '<span style="font-size:80%%"><b>Estimated CDS Ka/Ks ratio:&nbsp;&nbsp;%.4f (Ka: %.4f, Ks: %.4f)</b></span><br><br>' %(KaKs, Ka, Ks)
                        else:
                            print '<span style="font-size:80%%"><b>Unable to estimate Ka/Ks ratio of CDS: Query and reference coding sequences of unequal size (%i vs. %i).</b></span><br><br>' %(len(CDSpair[0]), len(CDSpair[1]))
                    else:
                        print '<span style="font-size:80%"><b>Unable to estimate Ka/Ks ratio of CDS: Unable to define CDS from user entered sequence.</b></span><br><br>'
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
                alignfile1 = open(tempfilename + hitdef + "_1",'wb')
                alignfile2 = open(tempfilename + hitdef +"_2",'wb')
                try:
                    alignfile1.write(">Query\n" + alignquery + "\n")
                    alignfile2.write(">Reference\n" + alignref + "\n")
                finally:
                    alignfile1.close()
                    alignfile2.close()
                os.environ["LAGAN_DIR"] = "../apps/lagan/"
                retcode = subprocess.call(["perl", "../apps/lagan/lagan.pl", tempfilename + hitdef + "_1", tempfilename + hitdef + "_2", "-out", tempfilename + hitdef + ".mfa" ,"-mfa"])
                retcode = subprocess.call(["perl", "../apps/lagan/mf_to_align.pl", "-f", tempfilename + hitdef + ".mfa", "-out", tempfilename + hitdef + "aligned"])
                alignedfile = open(tempfilename + hitdef + "aligned",'rb')
                try:
                    query = alignedfile.readline().rstrip('\r\n')
                    match = alignedfile.readline().rstrip('\r\n')
                    ref = alignedfile.readline().rstrip('\r\n')
                    window = form["window"].value
                    windowfilename = tempfilename + hitdef+ 'w'+ window
                    grad = []
                    gradnames = ["black", "green", "yellow", "orange","red"]
                    for i in range(0, 5):
                        if form.has_key("g"+gradnames[i]):
                            grad.append(float(form["g"+gradnames[i]].value)/100)
                        else:
                            grad.append(None)
                    if form.has_key("refheatmaps"):
                        makeheatmap(windowfilename + 'r.ppm', homologyvector(ref, match, int(window)), reffeatures, grad, heatmapsize, 25)
                    makeheatmap(windowfilename + '.ppm', homologyvector(query, match, int(window)), features, grad, heatmapsize, 25)
                finally:
                    alignedfile.close()
                if form.has_key("refheatmaps"):
                    retcode = subprocess.call(["../apps/png.py", "-c", "9", windowfilename+'r.ppm'])
                    os.remove(windowfilename+'r.ppm')
                retcode = subprocess.call(["../apps/png.py", "-c", "9", windowfilename+'.ppm'])
                os.remove(windowfilename+'.ppm')
                print 'The unformatted aligned sequence can be viewed <a href=../data/' + tempfilename[8:] + hitdef + 'aligned>here.</a><br>'
                print '<br><span style="font-size:80%"><b>Heatmap representation of ' + queryname + ':</b></span><br>'
                print '<br><img src=../data/' + windowfilename[8:] + '.png' + '><br><br>'
                if form.has_key("refheatmaps"):
                    print '<br><span style="font-size:80%"><b>Heatmap representation of ' + hitdef + ':</b></span><br>'
                    print '<br><img src=../data/' + windowfilename[8:] + 'r.png' + '><br><br>'
            print 'The heatmap diagrams can be downloaded by right-clicking on the image and \
                    selecting \"Save Picture As...\" (or \"Save Image As...\" in Firefox)'
            print '</body></html>'
        else:
            print [alignquery]
            print '<br>'
            print 'Error: Unable to parse query sequence. Invalid characters detected. \n \
                Click <a href="javascript:history.go(-1)"> here </a> to return to the input page.'
            print '</body></html>'


weboutput()
