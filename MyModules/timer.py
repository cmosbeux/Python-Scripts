#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  3 12:25:43 2017

@author: cmosbeux
"""

import time

class timer:
    @staticmethod
    def update(iteration, ntot, output_file='output.progress'):
        with open(output_file, 'w') as f:
            start_time = tinit
            message='Computing progress: (00000 %) [time running : 00000 hrs] [estimated : 00000 hrs]'
            len_message=len(message)
            f.write(message)
            f.seek(21)
            ratio=float(iteration)/ntot
            f.write(str(ratio*100)[0:5])
            f.seek(44)
            spent_time=(time.time() - start_time)/3600.
            f.write(str(spent_time)[0:5])
            f.seek(68)
            estim=spent_time/(ratio+1e-20)-spent_time
            f.write(str(estim)[0:5])

    @staticmethod
    def init():
        global tinit
        tinit = time.time()
        return tinit
