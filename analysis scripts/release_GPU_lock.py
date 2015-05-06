# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 14:24:58 2014

@author: gu98big
"""

#cd /
#cd /tmp
#ls -l gpus.lock


from subprocess import call
LOCK_FILE          = '/tmp/gpus.lock'
call(['rm', '-f', LOCK_FILE])
