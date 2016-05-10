uPEP notes
==========

Overview
--------

See: http://upep-scmb.biosci.uq.edu.au

With this document and the following commands you should be close to getting a functioning uPEP install:
    * git clone https://github.com/mscook/uPEP-helpers.git
    * git clone https://github.com/mscook/uPEP.git


Server configuration (UQ install)
---------------------------------

Performed the following:
    * Provisioned up a Virtual Machine on the Beatson r720 (1 vCPU, 2GB RAM, 250 GB space, running Ubuntu 12.04.2 LTS)
    * At present only accounts are: sciitadmin (Leslie Elliot) and uqmstan1
    * Had ITS/ScienceIT register http://upep-scmb.biosci.uq.edu.au/
    * Had http://upep-scmb.biosci.uq.edu.au/ added to the UQ exclusions list (world accessible)


Webserver (http://upep-scmb.biosci.uq.edu.au/) configuration
------------------------------------------------------------

Performed the following:
    * Renamed default vitrualhost to upep
    
Configuration looks like this::

    <VirtualHost *:80>                                                                                                                                                                                                                          
        ServerAdmin m.stantoncook@uq.edu.au                                     
                                                                                
        DocumentRoot /var/www                                                   
        <Directory />                                                           
                Options FollowSymLinks                                          
                AllowOverride None                                              
        </Directory>                                                            
        <Directory /var/www/>                                                   
                Options Indexes FollowSymLinks MultiViews                       
                AllowOverride None                                              
                Order allow,deny                                                
                allow from all                                                  
        </Directory>                                                            
        <Directory /var/www/data/>                                              
                Options -Indexes -FollowSymLinks -MultiViews                    
                AllowOverride None                                              
                Order allow,deny                                                
                allow from all                                                  
        </Directory>                                                            
                                                                                
        ScriptAlias /cgi-bin/ /var/www/cgi-bin/                                 
        <Directory "/var/www/cgi-bin/">                                         
                AllowOverride None                                              
                Options +ExecCGI -MultiViews +SymLinksIfOwnerMatch              
                Order allow,deny                                                
                Allow from all                                                  
        </Directory>                                                            
                                                                                
        ErrorLog /var/log/apache2/error.log                                     
                                                                                
        # Possible values include: debug, info, notice, warn, error, crit,         
        # alert, emerg.                                                         
        LogLevel warn                                                           
                                                                                
        CustomLog /var/log/apache2/access.log combined                          
                                                                                
    Alias /doc/ "/usr/share/doc/"                                               
    <Directory "/usr/share/doc/">                                               
        Options Indexes MultiViews FollowSymLinks                               
        AllowOverride None                                                      
        Order deny,allow                                                        
        Deny from all                                                           
        Allow from 127.0.0.0/255.0.0.0 ::1/128                                  
    </Directory>                                                                
                                                                                
    </VirtualHost>     

uPEP lives in:
    * /var/www

RefSeq dbs (uPEP specific) live in:
    * /var/RefSeq


uPEP installation/configuration
-------------------------------

Installed (these seemed to be required for lagan app):
    * libgetopt-long*
    * libautobox-list-util-perl
    * bioperl
    * blast2 (should eventually move off legacy BLAST)


Re-compiled lagan::

    cd /var/www/apps/lagan
    make


uPEP modifications
------------------

Please refer to commit history (diffs for details). Summary:
    * Updated paths, report RefSeq version and clean whitespace
    * Enabled all databases
    * Possible missing dependency. Listed could not find in logs.. May need to be looked into more
    * Refactored ACLoc & GILoc methods & added debugging calls


uPEP helpers
------------

Scripts/helpers to automatically update the RefSeq/uPEP databases. See: http://github.com/mscook/uPE-helpers

**Setup as a monthly (1st day) cron job**
