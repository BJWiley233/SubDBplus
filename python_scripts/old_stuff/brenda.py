#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 28 09:40:16 2020

@author: coyote
"""


from zeep import Client
import hashlib

wsdl = "https://brenda-enzymes.org/soap/brenda_zeep.wsdl"
password = hashlib.sha256("**".encode("utf-8")).hexdigest()
client = Client(wsdl)
parameters = ( "bwiley4@jh.edu",password,"ecNumber*1.1.1.1","organism*Homo sapiens","kmValue*",
              "kmValueMaximum*","substrate*","commentary*","ligandStructureId*","literature*" )
resultString = client.service.getKmValue(*parameters)
print (resultString)
len(resultString)
resultString[10]
