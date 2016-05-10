from django.http import Http404
from django.shortcuts import render
from django.template import loader
from django.http import HttpResponse, HttpResponseRedirect
from .models import Refseqdb, Starting_codon, Refseqdb_ACC_GI_build_log, Refseqdb_blast_db_build_log
from . import helper
from django.core.urlresolvers import reverse
from django.conf import settings
from django.utils import timezone
from .tasks import upephelper_processing
import os
import _mysql



def index(request):
    latest_db_list = Refseqdb.objects.order_by('-input_date')
    codon_list = Starting_codon.objects.all()
    context = {'latest_db_list': latest_db_list, 'codons': codon_list}
    return render(request, 'index.html', context)
# Create your views here.

def processing(request):
    response = HttpResponse()
    key = request.POST['database']    
    codons = Starting_codon.objects.all()
    response.write('Starting codon list:')
    cl = []
    for codon in codons:
        if codon.input_text in request.POST:
            response.write(codon.input_text)
            cl.append(codon.input_text.replace("U", "t").replace("A", "a").replace("G", "g").replace("C", "c"))
    if cl is None:
        return response.write("Invalid request.")
    outpath = settings.UPEPHELPER_STAGING
    data_loc = settings.UPEPHELPER_DATABASE

    override_condition = []
    if 'override' in request.POST:
        override_condition.append('True')
    else:
        override_condition.append('False')
    r = upephelper_processing.delay(key, cl, override_condition[0])
    return HttpResponse("%s" % (override_condition[0]))

def result(request, refseqdb_input_text):
    db = refseqdb_input_text
    return HttpResponse(db)
