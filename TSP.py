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
import itertools
import numpy

def parseDataset(ds):
    csvReader = csv.reader(ds, delimiter=' ')
    parsedDataset = []
    for row in csvReader:
        parsedDataset.append({'num': row[0], 'lat': row[1], 'lon': row[2]})
    
    return parsedDataset

def generateChromosomes(ds, number):
    listChromosomes = []
    
    for c in range (number):
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

def fitness(chromosome, ds):
    dist = 0;
    for c in range(len(chromosome)):
        currentIndex = int(chromosome[c])
        latC = float(ds[currentIndex].get('lat'))
        lonC = float(ds[currentIndex].get('lon'))
        if(c == 0):
            dist += math.sqrt(math.pow(latC, 2)+math.pow(lonC, 2))
        else:
            previousIndex = int(chromosome[c-1])
            dist += math.sqrt(math.pow(latC-float(ds[previousIndex].get('lat')), 2)+math.pow(lonC-float(ds[previousIndex].get('lon')), 2))
    return 1/dist
    
def mutation(chromosome):
    chromosomeLength = len(chromosome)
    index1 = randrange(chromosomeLength-1)
    index2 = randrange(chromosomeLength-1)
    temp = chromosome[index2]
    chromosome[index2] = chromosome[index1]
    chromosome[index1] = temp
#    print(len(chromosome))
    return chromosome

def crossover(chromosome1, chromosome2):
    
    pivot = randrange(len(chromosome1)-2)
    newChromosome = []
    
   
    for c1 in range(pivot):
        newChromosome.append(chromosome1[c1])
        
    for c2 in range(len(chromosome1)-pivot):
        newChromosome.append(chromosome2[pivot+c2])
    
    return newChromosome
    
def generateFitnessList(chromosomeList, parsedDataset):
    fitnessList = []
    for f in range(len(chromosomeList)):
        fitnessList.append(fitness(chromosomeList[f], parsedDataset))
    return fitnessList
        
def extractBestN(chromosomeList, fitnessList, n):
    bestChromosomesList = []
    indexList = []
    
    highestList = sorted(fitnessList, reverse=True)[:n]
  
    for i in range(len(highestList)):
        indexList.append(fitnessList.index(highestList[i]))
        
    for j in indexList:
        bestChromosomesList.append(chromosomeList[j])

    return bestChromosomesList
        
def mixList(myList):
    
    listLength = len(myList)
    index1 = randrange(listLength-1)
    index2 = randrange(listLength-1)
    temp = myList[index2]
    myList[index2] = myList[index1]
    myList[index1] = temp
    
    return myList
    
def listCrossover(chromosomeList):
#    newGeneration = []
#    if(len(myList)>=2):
#        for a in range(len(myList)-2):
#            newGeneration.append(crossover(myList[a], myList[a+1]))
#    return newGeneration
#    print(len(chromosomeList))
    totDist = 0;
    crossList = []
    
    for fitness in fitnessList:
        totDist += 1/fitness
        
    for c in range(len(chromosomeList)-1):
        if((1/fitnessList[c])/totDist <= 2*randrange(math.floor(1000))/1000):
            crossList.append(crossover(chromosomeList[c], chromosomeList[c+1]))
#            crossList.append(crossover(chromosomeList[c+1], chromosomeList[c]))
        crossList.append(chromosomeList[c])
#            crossList.append(chromosomeList[c+1])
    crossList.append(chromosomeList[len(chromosomeList)-1])
    return crossList

def mutateGroup(chromosomeList, fitnessList, amount):
    
    totDist = 0;
    mutateList = []
    
    for fitness in fitnessList:
        totDist += 1/fitness
        
    for c in range(len(chromosomeList)):
        if((amount/fitnessList[c]) >= randrange(math.floor(1000))/1000):
            mutateList.append(mutation(chromosomeList[c]))
        
        mutateList.append(chromosomeList[c])

    return mutateList

#    mutatedChromosomes = []
#    for c in range(numberMutations):
#        mutatedChromosomes.append(mutation(chromosomeList[c]))
#    for d in range(len(chromosomeList)-numberMutations):
#        mutatedChromosomes.append(chromosomeList[d+numberMutations])
##    print(len(chromosomeList))
##    print(len(mutatedChromosomes))
#    return mutatedChromosomes
        
def kill(chromosomeList, parsedDataset, threshold):
    killList = []
    for c in range(len(chromosomeList)-1):
        if(fitness(chromosomeList[c], parsedDataset)<threshold):
            killList.append(chromosomeList[c])
            
    return killList
    
def calculateThreshold(fitnessList):
    fitnessList.sort()
    return fitnessList[math.floor(0.5*len(fitnessList))]

def removeEmpty(chromosomeList):
    return [x for x in chromosomeList if len(x)!=0]

def killWeak(chromosomeList, fitnessList):
    totDist = 0;
    killList = []
    
    for fitness in fitnessList:
        totDist += 1/fitness
        
    for c in range(len(chromosomeList)-1):
        if((1/fitnessList[c])/totDist >= randrange(math.floor(1000))/1000):
            killList.append(chromosomeList[c])
    
    return killList
    
def killNWeakest(chromosomeList, ds, n):
    fitnessList = []
    cleanedList = []
    for c in range(len(chromosomeList)):
        fitnessList.append(fitness(chromosomeList[c], ds))
    sortedInd = numpy.argsort(fitnessList)    
    
    for c in range(len(sortedInd)):
        if(sortedInd[c] > n):
            cleanedList.append(chromosomeList[sortedInd[c]])
        
    return cleanedList
    
    
# Upload a file
Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
filename = askopenfilename() # show an "Open" dialog box and return the path to the selected file

dataset = open(filename, "r")
parsedDataset = parseDataset(dataset)

chromosomeList = generateChromosomes(parsedDataset, 1000)

evolve = 1

bestChromosomeBefore = 0
bestChromosomeAfter = 0
amount = 0
cumulateSaved = 0

while(evolve):
    delta = -bestChromosomeBefore+bestChromosomeAfter

    numberBefore = len(chromosomeList)
    fitnessList = generateFitnessList(chromosomeList, parsedDataset)
    chromosomeList = mutateGroup(chromosomeList, fitnessList, amount)
    fitnessList = generateFitnessList(chromosomeList, parsedDataset)

#    print(len(chromosomeList))

#    killList = killWeak(chromosomeList, fitnessList)
##    print(killList)
#    set1 = {tuple(item) for item in chromosomeList}
#    chromosomeList = [list(item) for item in set1 if tuple(item) not in killList]
#    
##    chromosomeList = set(chromosomeList) ^ set(killList)
    
    bestChromosomes = extractBestN(chromosomeList, fitnessList, math.ceil(len(chromosomeList)/1))
    
    bestChromosomeBefore = bestChromosomeAfter
    bestChromosomeAfter = 1/fitness(bestChromosomes[0], parsedDataset)
    amount = 100*math.exp(10*(-delta)/(bestChromosomeBefore+1))
    

#    print("AMOUNT: " + str(amount) + " ... DELTA " + str(-bestChromosomeBefore+bestChromosomeAfter))
#    print(amount)
#    amount = 100000/(bestChromosomeAfter - bestChromosomeBefore + 1)
    
    print(1/fitness(bestChromosomes[0], parsedDataset))
    
#    print(bestChromosomes[0])
    
#    threshold = calculateThreshold(fitnessList)
#    killList = []
#    killList = kill(chromosomeList, parsedDataset, threshold)
#    
#    set1 = {tuple(item) for item in chromosomeList}
#    chromosomeList = [item for item in set1 if tuple(item) not in killList]
#    
#    print(len(chromosomeList))
    
#    mixedBestChromosomes = mixList(bestChromosomes)
    chromosomeList = listCrossover(chromosomeList)  
    numberAfter = len(chromosomeList)
    deltaAlive = numberAfter - numberBefore - 1
    numberKill = max(0.9*deltaAlive, min(deltaAlive * math.pow(abs(delta/(bestChromosomeAfter+1)),0.5), deltaAlive))
    cumulateSaved += deltaAlive - numberKill
    if(cumulateSaved*delta/(bestChromosomeBefore+1)>1):
        numberKill += cumulateSaved
        cumulateSaved = 0

    chromosomeList = killNWeakest(chromosomeList, parsedDataset, numberKill)
    
#    print(len(chromosomeList))
#    chromosomeList = list(dict.fromkeys(chromosomeList))
    
#    chromosomeList = chromosomeList + newGeneration
    
#    chromosomeList = mixList(chromosomeList)
        
#    numberMutations = math.ceil(randrange(len(chromosomeList))/2)
    
#    print(len(bestChromosomes))
#    print(len(chromosomeList))
#    print(len(newGeneration))
#    print(numberMutations)
#    print("\n")
#    print(numberMutations)
#    print(len(chromosomeList))
    
    

    
    #print(chromosomeList)
    
    
#    print(bestChromosome[1])
#    print(1/fitness(bestChromosome[1], parsedDataset))
#    
#    


#print(chromosomeList[0])
#print("\n")
#print(chromosomeList[1])
#print("\n\n")
#print(crossover(chromosomeList[0], chromosomeList[1]))