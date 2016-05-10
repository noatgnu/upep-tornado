#!/usr/bin/env python
import cgi, subprocess, time, random, os, sys
from upeplib import images, accessions


def percIdent(alignstring):
    return alignstring.count('|')/float(len(alignstring))


def homologyvector(query, align, calcwindow):
    returnlist = []
    for i in range(0, len(query)):
        if query[i] == '-':
            continue
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
    for i in range(0,6):
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

def weboutput():
    tempfilename =  '../data/%012x%016x' %(int(time.time()*1000), random.randint(0,0xFFFFFFFFFFFFFFFF))
    form = cgi.FieldStorage()
    if form.keys() == []:
        print "Refresh: 5; url=../../"
        print "Content-type: text/html"
        print 
        print 'You will be redirected to uPEPperoni in 5 seconds. \n\
        Please click <a href="../../">here</a> if not automatically redirected.'
        return 0
    print "Content-type: text/html"
    ##print "Content-disposition: attachment; filename=fname.ext"
    ##print "Set-Cookie: " +tempfilename[8:]+ "=0\n"
    print ## <----- VERY IMPORTANT!!!!!! Blank "print" - Separate Headers
    print '<html><head>'
    print '<title>'
    print 'uPEPperoni - Results'
    print '</title>'
    print '</head><body>'
    ##print form
    ##print '<br><br>'
    features = []
    if form.has_key("alignquery") and form.has_key("alignref"):
        if (form["alignquery"].value[:3] == "NM_") or (form["alignquery"].value[:3] == "XM_") or (form["alignquery"].value[:3] == "GI:"):
            mRNAparams = [0, 0, 0]
            if form.has_key("min"):
                mRNAparams[0] = int(form["min"].value)
            if form.has_key("max"):
                mRNAparams[1] = int(form["max"].value)
            if form.has_key("grace"):
                mRNAparams[2] = int(form["grace"].value)
            details = accessions.getmRNA(form["alignquery"].value, mRNAparams[0]*3, mRNAparams[1]*3, mRNAparams[2])
            if details:
                alignquery = details[0]
                if form.has_key("drawcds"):
                    features.append(details[1] + ['CDS'])
                for i in details[2]:
                    features.append(i + ['uORF'])
                print '<b>QUERY: ' + details[3] + '</b><br>'
            else:
                print '<br>' 
                print 'Error: Unable to find query mRNA sequence in selected Refseq database. ' + form["alignquery"].value + ' not found \n \
                Click <a href="javascript:history.go(-1)"> here </a> to return to the input page.'
                print '</body></html>'
                return 0
        else:
            print '<b>QUERY: User Entered</b><br>'
            alignquery = form["alignquery"].value.upper().replace('\n','').replace('\r','').replace(' ','')
        if (form["alignref"].value[:3] == "NM_") or (form["alignref"].value[:3] == "XM_") or (form["alignref"].value[:3] == "GI:"):
            details = accessions.getmRNA(form["alignref"].value, 0, 0, 0)
            if details:
                alignref = details[0]
                print '<b>REFERENCE: ' + details[3] + '</b><br>'
            else:
                print '<br>' 
                print 'Error: Unable to find reference mRNA sequence in selected Refseq database. ' + form["alignref"].value + ' not found \n \
                Click <a href="javascript:history.go(-1)"> here </a> to return to the input page.'
                print '</body></html>'
                return 0
        else:
            print '<b>REFERENCE: User Entered</b><br>'
            alignref = form["alignref"].value.upper().replace('\n','').replace('\r','').replace(' ','')
        if seqcheck(alignquery, ('A', 'T', 'C', 'G', 'N')):
            if seqcheck(alignref, ('A', 'T', 'C', 'G', 'N')):
                alignfile1 = open(tempfilename+"_1",'wb')
                alignfile2 = open(tempfilename+"_2",'wb')
                try:
                    alignfile1.write(">Query\n" + alignquery + "\n")
                    alignfile2.write(">Reference\n" + alignref + "\n")
                finally:
                    alignfile1.close()
                    alignfile2.close()
                os.environ["LAGAN_DIR"] = "../apps/lagan/"
                retcode = subprocess.call(["perl", "../apps/lagan/lagan.pl", tempfilename+"_1", tempfilename+"_2", "-out", tempfilename+".mfa" ,"-mfa"])
                ##print [form["match"].value, form["mismatch"].value, form["existence"].value, form["extension"].value]
                retcode = subprocess.call(["perl", "../apps/lagan/mf_to_align.pl", "-f", tempfilename+".mfa", "-out", tempfilename+".aligned"])
                alignedfile = open(tempfilename+".aligned",'rb')
                try:
                    query = alignedfile.readline().rstrip()
                    match = alignedfile.readline().rstrip()
                    ref = alignedfile.readline().rstrip()
                    window = form["window"].value
                    windowfilename = tempfilename+ 'w'+ window
                    if form.has_key("markers"):
                        if type(form["markers"]) == list:
                            for i in form["markers"]:
                                firstsplit = i.value.index(',')
                                secondsplit = i.value.index(']')
                                temp = [int(i.value[1:firstsplit]), int(i.value[firstsplit+2:secondsplit]), i.value[secondsplit+2:]]
                                features.append(temp)
                        else:
                            firstsplit = form["markers"].value.index(',')
                            secondsplit = form["markers"].value.index(']')
                            features.append([int(form["markers"].value[1:firstsplit]), int(form["markers"].value[firstsplit+2:secondsplit]), form["markers"].value[secondsplit+2:]])
                    grad = []
                    gradnames = ["black", "green", "yellow", "orange", "red", "blue"]
                    for i in range(0, 6):
                        if form.has_key("g"+gradnames[i]):
                            grad.append(float(form["g"+gradnames[i]].value)/100)
                        else:
                            grad.append(None)
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
                    makeheatmap(windowfilename + '.ppm', homologyvector(query, match, int(window)), features, grad, heatmapsize, 25)
                finally:
                    alignedfile.close()
                retcode = subprocess.call(["../apps/png.py", "-c", "9", windowfilename+'.ppm'])
                os.remove(windowfilename+'.ppm')
                print 'Query ID:' + tempfilename[8:] + '<br>'
                #print 'The aligned sequence can be viewed <a href="../data/servefile.php?file=' + tempfilename[8:] + '"> Save text file</a><br>'
                print 'The unformatted aligned sequence can be viewed <a href=../data/' + tempfilename[8:]+ '.aligned>here.</a><br>'
                print 'The heatmap diagram below can be downloaded by right-clicking on the image and \
                    selecting \"Save Picture As...\" (or \"Save Image As...\" in Firefox)'  
                print '<br><img src=../data/' + windowfilename[8:] + '.png' + '></body></html>'
            else:
                print [alignref]
                print '<br>' 
                print 'Error: Unable to parse reference sequence. Invalid characters detected. \n \
                Click <a href="javascript:history.go(-1)"> here</a> to return to the input page.'
                print '</body></html>'
        else:
            print [alignquery]
            print '<br>'
            print 'Error: Unable to parse query sequence. Invalid characters detected. \n \
                Click <a href="javascript:history.go(-1)"> here</a> to return to the input page.'
            print '</body></html>'
    else:
        print 'Invalid sequences entered. Click <a href="../../"> here</a> to return to the input page.</body></html>'

weboutput()

