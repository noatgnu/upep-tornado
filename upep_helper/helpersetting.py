import os

APP_ROOT = os.path.dirname(os.path.abspath(__file__))
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
STARTING_CODONS = ['ttg','ctg','att','atc','ata','gtg','atg','acg','aag','agg']
REFSEQ_DBS = ['RefSeq-complete', 'RefSeq-fungi', 'RefSeq-invertebrate', 'RefSeq-plant', 'RefSeq-vertebrate_mammalian', 'RefSeq-vertebrate_other']
