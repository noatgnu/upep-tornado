import time
import accessions
import subprocess
import upepsetting
import os
import random
import conserveduPEP

def whmstandalone(body_arg):
    parameters = dict()
    bcodon = body_arg["codon"][0]
    parameters["codon"] = bcodon.decode("utf-8")
    bdbv = body_arg["database_version"][0]
    parameters["dbv"] = bdbv.decode("utf-8")
    parameters['alignquery'] = body_arg["alignquery"][0].decode("utf-8")
    parameters['alignref'] = body_arg["alignref"][0].decode("utf-8")
    parameters["match"] = 0
    #print(body_arg)
    if "storage" in body_arg:
        parameters["markers"] = body_arg['storage'][0]
    if "match" in body_arg:
        parameters["match"] = body_arg['match'][0].decode("utf-8")
    parameters["mismatch"] = 0
    if "mismatch" in body_arg:
        parameters["mismatch"] = body_arg['mismatch'][0].decode("utf-8")
    parameters["existence"] = 0
    if "existence" in body_arg:
        parameters["existence"] = body_arg['existence'][0].decode("utf-8")
    parameters["existension"] = 0
    if "extension" in body_arg:
        parameters["extension"] = body_arg['extension'][0].decode("utf-8")
    parameters["min"] = 0
    if "min" in body_arg:
        bmin = body_arg["min"][0]
        parameters["min"] = bmin.decode("utf-8")
    parameters["max"] = 0
    if "min" in body_arg:
        bmax = body_arg["max"][0]
        parameters["max"] = bmax.decode("utf-8")
    parameters["grace"] = 0
    if "min" in body_arg:
        bgrace = body_arg["grace"][0]
        parameters["grace"] = bgrace.decode("utf-8")
    gradi = []
    gradnames = ["black", "blue", "green", "yellow", "red"]
    for i in range(0, 5):
        if "g" + gradnames[i] in body_arg:
            if not body_arg["g" + gradnames[i]][0] == b"":
                bcolor = body_arg["g" + gradnames[i]][0]
                color = bcolor.decode("utf-8")
                gradi.append(float(color) / 100)
            else:
                gradi.append(None)
    parameters["grad"] = gradi
    parameters["heatmapsize"] = 0
    if "heatmapsize" in body_arg:
        bhms = body_arg["heatmapsize"][0]
        parameters["heatmapsize"] = bhms.decode("utf-8")
    parameters["window"] = 0
    if "window" in body_arg:
        bws = body_arg["window"][0]
        parameters["window"] = bws.decode("utf-8")
    upep_loc = upepsetting.APP_ROOT
    apps_loc = upepsetting.APPS
    lagan_loc = upepsetting.LAGAN
    
    upeplib_loc = upepsetting.UPEPLIB
    data_loc = upepsetting.DATA
    results = dict()
    unid = '%012x%016x' %(int(time.time()*1000), random.randint(0,0xFFFFFFFFFFFFFFFF))
    
    results['unid'] = unid
    tempfilename = os.path.join(data_loc, unid)
    results['tfn'] = tempfilename
    dbstuff = parameters["dbv"] + upepsetting.DBRE["Complete"]
    print(parameters['markers'])
    features = []
    if "alignquery" in parameters and "alignref" in parameters:
        if (parameters["alignquery"][:3] == "NM_") or (parameters["alignquery"][:3] == "XM_") or (parameters["alignquery"][:3] == "GI:"):
            mRNAparams = [0, 0, 0]
            if "min" in parameters:
                mRNAparams[0] = parameters["min"]
            if "max" in parameters:
                mRNAparams[1] = parameters["max"]
            if "grace" in parameters:
                mRNAparams[2] = parameters["grace"]
            details = accessions.getmRNA(parameters["alignquery"], dbstuff, parameters["codon"], mRNAparams[0]*3, mRNAparams[1]*3, mRNAparams[2])
            if details:
                alignquery = details[0]
                if "drawcds" in parameters:
                    features.append(details[1] + ['CDS'])
                if "uorfs" in parameters:
                    for i in details[2]:
                        features.append(i + ['uORF'])
                #print('<b>QUERY: ' + details[3] + '</b><br>')
                results['query'] = details[3]
            else:
                results['error'] = 'Unable to find query mRNA sequence in selected Refseq database.'
                #print('<br>') 
                #print('Error: Unable to find query mRNA sequence in selected Refseq database. ' + parameters["alignquery"] + ' not found \n \
                #Click <a href="javascript:history.go(-1)"> here </a> to return to the input page.')
                #print('</body></html>')
                return results
        else:
            #print('<b>QUERY: User Entered</b><br>')
            results['query'] = 'User Entered Sequence'
            alignquery = parameters["alignquery"].upper().replace('\n','').replace('\r','').replace(' ','')
        if (parameters["alignref"][:3] == "NM_") or (parameters["alignref"][:3] == "XM_") or (parameters["alignref"][:3] == "GI:"):
            details = accessions.getmRNA(parameters["alignref"], dbstuff, parameters["codon"], 0, 0, 0)
            if details:
                alignref = details[0]
                #print('<b>REFERENCE: ' + details[3] + '</b><br>')
                results['reference'] = details[3]
            else:
                results['error'] = 'Unable to find query mRNA sequence in selected Refseq database.'
                return results
                #print('<br>') 
                #print('Error: Unable to find reference mRNA sequence in selected Refseq database. ' + form["alignref"].value + ' not found \n \
                #Click <a href="javascript:history.go(-1)"> here </a> to return to the input page.')
                #print('</body></html>')
                #return 0
        else:
            #print('<b>REFERENCE: User Entered</b><br>')
            results['reference'] = 'User Entered Sequence'
            alignref = parameters["alignref"].upper().replace('\n','').replace('\r','').replace(' ','')
        if conserveduPEP.seqcheck(alignquery, ('A', 'T', 'C', 'G')):
            if conserveduPEP.seqcheck(alignref, ('A', 'T', 'C', 'G')):
                alignfile1 = open(tempfilename+"_1",'wt')
                alignfile2 = open(tempfilename+"_2",'wt')
                try:
                    alignfile1.write(">Query\n" + alignquery + "\n")
                    alignfile2.write(">Reference\n" + alignref + "\n")
                finally:
                    alignfile1.close()
                    alignfile2.close()
                os.environ["LAGAN_DIR"] = lagan_loc
                retcode = subprocess.call(["perl", lagan_loc + "/lagan.pl", tempfilename+"_1", tempfilename+"_2", "-out", tempfilename+".mfa" ,"-mfa"])
                ##print [form["match"].value, form["mismatch"].value, form["existence"].value, form["extension"].value]
                retcode = subprocess.call(["perl", lagan_loc + "/mf_to_align.pl", "-f", tempfilename+".mfa", "-out", tempfilename+"_aligneduf"])
                window = parameters["window"]
                windowfilename = tempfilename+ 'w'+ window
                results['windowpng'] = unid + 'w'+ window + '.png'
                with open(tempfilename + "_aligneduf", "rt") as aligneduffile, open(tempfilename + "_aligned", "wt") as alignedfile: 
                    query = "".join(aligneduffile.readline().rstrip("\r\n"))
                    match = "".join(aligneduffile.readline().rstrip("\r\n"))
                    ref = "".join(aligneduffile.readline().rstrip("\r\n"))
                    query2 = conserveduPEP.split_sequence(query)
                    match2 = conserveduPEP.split_sequence(match)
                    ref2 = conserveduPEP.split_sequence(ref)
                    for i1, i2, i3 in zip(query2, match2, ref2):
                        alignedfile.write("%s\n%s\n%s\n\n" % (i1, i2, i3))
                        
                    if "markers" in parameters:
                        markers = parameters["markers"].decode('utf-8').split('\r\n')
                        print(parameters["markers"].decode('utf-8').split('\r\n'))
                        print(type(markers))
                        if type(markers) == list:
                            for i in markers:
                                firstsplit = i.index(',')
                                secondsplit = i.index(']')
                                temp = [int(i[1:firstsplit]), int(i[firstsplit+2:secondsplit]), i[secondsplit+2:]]
                                print(temp)
                                features.append(temp)
                        else:
                            firstsplit = parameters["markers"].index(',')
                            secondsplit = parameters["markers"].index(']')
                            features.append([int(parameters["markers"][1:firstsplit]), int(parameters["markers"][firstsplit+2:secondsplit]), parameters["markers"][secondsplit+2:]])
                    if "heatmapsize" not in parameters:
                        results['error'] = 'The width of the heatmap was not specified.'
                        #print('<br>') 
                        #print('Error: The width of the heatmap was not specified.  \
                        #Click <a href="javascript:history.go(-1)"> here </a> to return to the input page.')
                        #print('</body></html>')
                        return results
                    heatmapsize = int(parameters["heatmapsize"])
                    if (heatmapsize < 400) or (heatmapsize > 10000):
                        results['error'] = 'The heapmap width must be between 400 and 10000 pixels.'
                        #print('<br>') 
                        #print('Error: The heapmap width must be between 400 and 10000 pixels.  \
                        #Click <a href="javascript:history.go(-1)"> here </a> to return to the input page.')
                        #print('</body></html>')
                        return results
                    print(features)
                    conserveduPEP.makeheatmap(windowfilename + '.ppm', conserveduPEP.homologyvector(query, match, int(window)), features, parameters["grad"], heatmapsize, 25)
                
                retcode = subprocess.call([apps_loc+"/png.py", "-c", "9", windowfilename+'.ppm'])
                results['windowpng'] = unid + 'w'+ window + '.png'
                os.remove(windowfilename+'.ppm')
                return results
                #print('Query ID:' + unid + '<br>')
                #print 'The aligned sequence can be viewed <a href="../data/servefile.php?file=' + tempfilename[8:] + '"> Save text file</a><br>'
                #print('The unformatted aligned sequence can be viewed <a href=' + os.path.join(data_loc, unid) + '_aligned>here.</a><br>')
                #print('The heatmap diagram below can be downloaded by right-clicking on the image and \
                    #selecting \"Save Picture As...\" (or \"Save Image As...\" in Firefox)')  
                #print('<br><img src='+ os.path.join(data_loc, unid) + '.png' + '></body></html>')
            else:
                
                results['error'] = 'Unable to parse reference sequence.'
                #print([alignref])
                #print('<br>')
                return results
                #print('Error: Unable to parse reference sequence. Invalid characters detected. \n \
                #Click <a href="javascript:history.go(-1)"> here</a> to return to the input page.')
                #print('</body></html>')
        else:
            results['error'] = 'Unable to parse reference sequence.'
            #print([alignquery])
            #print('<br>')
            #print('Error: Unable to parse query sequence. Invalid characters detected. \n \
                #Click <a href="javascript:history.go(-1)"> here</a> to return to the input page.')
            #print('</body></html>')
            return results
    else:
        results['error'] = 'Invalid sequences entered'
        #print('Invalid sequences entered. Click <a href="../../"> here</a> to return to the input page.</body></html>')
        return results

