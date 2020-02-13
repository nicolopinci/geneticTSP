#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 12:52:26 2020

@author: Nicolò Pinciroli
"""

from tkinter import Tk
from tkinter.filedialog import askopenfilename
import csv
import math
from random import randrange

def parseDataset(ds):
    csvReader = csv.reader(ds, delimiter=' ')
    parsedDataset = []
    for row in csvReader:
        parsedDataset.append({'num': row[0], 'lat': row[1], 'lon': row[2]})
    
    return parsedDataset

def generateChromosomes(ds, factor):
    numberCities = len(ds)
    numberChromosomes = math.ceil(numberCities/factor)
    listChromosomes = []
    
    for c in range (numberChromosomes):
        randomChromosome = newChromosome(ds)
        listChromosomes.append(randomChromosome)
    
    return listChromosomes

def newChromosome(ds):
    chromosomeLength = len(ds)
    baseChromosome = list(range(1, chromosomeLength))
    for l in range(chromosomeLength):
        index1 = randrange(chromosomeLength-1)
        index2 = randrange(chromosomeLength-1)
        temp = baseChromosome[index2]
        baseChromosome[index2] = baseChromosome[index1]
        baseChromosome[index1] = temp
    return baseChromosome

# Upload a file
Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
filename = askopenfilename() # show an "Open" dialog box and return the path to the selected file

dataset = open(filename, "r")
parsedDataset = parseDataset(dataset)

print(generateChromosomes(parsedDataset, 4))