# -*- coding: utf-8 -*-
"""
Created on Mon May 13 20:40:59 2013

@author: root
"""
# import pdb
class FeaturesContainer(dict):
    def computeKernel(self,f):
        return self['mi']-f['mi']        
    
if __name__=="__main__":
    F=FeaturesContainer()
    F['mi']=2
    F2=FeaturesContainer()
    F2['mi']=2