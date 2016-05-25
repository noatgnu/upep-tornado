uPEP notes
==========

Version 2.0:
Requirements:
    * Python >= 2.7
    * MySQL >= 5.6
    * Redis
    * LAGAN Dependencies:
        + libgetopt-long*
        + libautobox-list-util-perl
        + bioperl
        + blast2
    * Python lib:
        + tornado
        + mysqlclient
        + futures
        + ftputil
        + redis
        + rq
        + tqdm

Re-compiled lagan::

    cd <lagan directory>
    make

Database building info:
    - Change the file location in the variables in /upep/upepsetting.py and /upep_helper/helpersetting.py:
        + UPEPHELPER_STAGING (temporary location to download the database)
        + UPEPHELPER_DATABASE (final location)
    - To build database:
        + cd <upep-tornado directory>
        + python
        + from upep_helper import helpersetting
        + codons = helpersetting.STARTING_CODONS
        + key = <one of the following 'RefSeq-complete', 'RefSeq-fungi', 'RefSeq-invertebrate', 'RefSeq-plant', 'RefSeq-vertebrate_mammalian', 'RefSeq-vertebrate_other'>
        + from upep_helper import helpertask
        + helpertasks.upephelper_processing(key, codons, 'True')
                                    
