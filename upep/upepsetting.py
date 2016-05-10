import os

APP_ROOT = os.path.dirname(os.path.abspath(__file__))
APPS = os.path.join(APP_ROOT, "apps")
LAGAN = os.path.join(APPS, "lagan")
UPEPCGI = os.path.join(APP_ROOT, "cgi-bin")
UPEPLIB = os.path.join(UPEPCGI, "upeplib")
DATA = os.path.join(APP_ROOT, "data")

UPEPHELPER_STAGING = '/media/refseqs/SCMB/Research/Rothnagel_J/uPEPdb/RefSeq/Staging/'
UPEPHELPER_DATABASE = '/media/refseqs/SCMB/Research/Rothnagel_J/uPEPdb/RefSeq/'
DATABASES = {
    'default': {
        'NAME': 'upep',
        'USER': 'root',
        'PASSWORD': 'upep2016',
        'HOST': 'localhost',
        'PORT': '',
        'DB' : 'uPEP',
    }
}
STARTING_CODONS = ['aug','uug','cug','auu','auc','aua','gug','acg','aag','agg']
REFSEQ_DBS = ['RefSeq-complete', 'RefSeq-fungi', 'RefSeq-invertebrate', 'RefSeq-plant', 'RefSeq-vertebrate_mammalian', 'RefSeq-vertebrate_other']
MAMMAL = 'vertebrate_mammalian'
NON_MAMMALIAN_VERTEBRATES = 'vertebrate_other'
INVERTEBRATES = 'invertebrate'
PLANTS = 'plant'
FUNGI = 'fungi'
COMPLETE = 'complete' 
DBWEB = ['Complete','Mammals','Plants','Non-mammalian vertebrates','Invertebrates','Fungi']
DBRE = {'Mammals': MAMMAL,
               'Non-mammalian vertebrates': NON_MAMMALIAN_VERTEBRATES,
               'Invertebrates': INVERTEBRATES,
               'Plants': PLANTS,
               'Fungi': FUNGI,
               'Complete': COMPLETE}
