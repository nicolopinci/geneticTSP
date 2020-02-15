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
import numpy
import random
import pygame

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

def generateGreedyChromosomes(ds, number):
    listChromosomes = []
    
    number = min(number, len(ds))
    for k in range(1, number):
        city = ds[k-1]
        listChromosomes.append(generateGreedyChromosome(city.get('num'), ds))
    
    return listChromosomes

def generateGreedyChromosome(startingPoint, ds):
    numberCities = len(ds)+1
    baseChromosome = list(range(1, numberCities))
    greedyChromosome = []
    
    greedyChromosome.append(int(startingPoint))
    baseChromosome.remove(int(startingPoint))
    
    for g in range(numberCities-2):
        closestPoint = findClosest(greedyChromosome[-1], baseChromosome, ds)[0]
        greedyChromosome.append(closestPoint)
      
        baseChromosome.remove(int(closestPoint))
        
#    print(greedyChromosome)
    return greedyChromosome


def generateAlmostGreedyChromosomes(ds, number):
    listChromosomes = []
    
    for l in range(0, min(1, math.ceil(number/len(ds)))):
        for k in range(1, len(ds)-1):
            city = ds[k-1]
            listChromosomes.append(generateAlmostGreedyChromosome(city.get('num'), ds))
    
    return listChromosomes

def generateAlmostGreedyChromosome(startingPoint, ds):
    numberCities = len(ds)+1
    baseChromosome = list(range(1, numberCities))
    greedyChromosome = []
    
    greedyChromosome.append(int(startingPoint))
    baseChromosome.remove(int(startingPoint))
    
    for g in range(numberCities-2):
        randomNumber = min(math.floor(1+min((2.99*random.random()), 2)), len(greedyChromosome))
        closestPoint = findClosest(greedyChromosome[-randomNumber], baseChromosome, ds)[0]
        greedyChromosome.append(closestPoint)
      
        baseChromosome.remove(int(closestPoint))
        
#    print(greedyChromosome)
    return greedyChromosome


def findClosest(element, baseChromosome, ds):
    minDist = float("inf")
    minElem = 0
    currentLat = -1
    currentLon = -1
    
    for el in ds:
        if(int(el.get('num'))==int(element)):
            currentLat = float(el.get('lat'))
            currentLon = float(el.get('lon'))
            
    for elB in baseChromosome:
        newLat = -1
        newLon = -1
        for elem in ds:
           if(int(elem.get('num'))==elB):
               newLat = float(elem.get('lat'))
               newLon = float(elem.get('lon'))
            
        if(math.sqrt(math.pow(newLat-currentLat,2)+math.pow(newLon-currentLon,2))<minDist):
            minElem = int(elB)
            minDist = math.sqrt(math.pow(newLat-currentLat,2)+math.pow(newLon-currentLon,2))
    
    return [minElem, minDist]
    

def newChromosome(ds):
    chromosomeLength = len(ds)+1
    baseChromosome = list(range(1, chromosomeLength))
    for l in range(chromosomeLength):
        index1 = randrange(chromosomeLength-1)
        index2 = randrange(chromosomeLength-1)
        temp = baseChromosome[index2]
        baseChromosome[index2] = baseChromosome[index1]
        baseChromosome[index1] = temp
    return baseChromosome

def fitness(chromosome, ds):
    return 1/distance(chromosome, ds, False)
    
def mutation(chromosome):
    chromosomeLength = len(chromosome)
    index1 = randrange(chromosomeLength-1)
    index2 = randrange(chromosomeLength-1)
    temp = chromosome[index2]
    chromosome[index2] = chromosome[index1]
    chromosome[index1] = temp
    return chromosome

def multipleMutation(chromosome, percentage):
    multipleMutatedChromosome = chromosome
    for m in range(0, max(math.floor(percentage*len(chromosome)-1), 1)):
        multipleMutatedChromosome = mutation(multipleMutatedChromosome)
    return multipleMutatedChromosome
        
def crossover(parent1, parent2): # source: https://towardsdatascience.com/evolution-of-a-salesman-a-complete-genetic-algorithm-tutorial-for-python-6fe5d2b3ca35
    

    child = []
    childP1 = []
    childP2 = []
    
    geneA = int(random.random() * len(parent1))
    geneB = int(random.random() * len(parent1))
    
    startGene = min(geneA, geneB)
    endGene = max(geneA, geneB)

    for i in range(startGene, endGene):
        childP1.append(parent1[i])
        
    childP2 = [item for item in parent2 if item not in childP1]

    child = childP1 + childP2
        
    return [child, child]


def oldCrossover(chromosome1, chromosome2):
    
    pivot = randrange(len(chromosome1)-2)
    newChromosome = []
    
   
    for c1 in range(pivot):
        newChromosome.append(chromosome1[c1])
        
    for c2 in range(len(chromosome1)-pivot-1):
        newChromosome.append(chromosome2[pivot+c2])
    

            
    if(isValid(newChromosome) and len(newChromosome) == len(chromosome1)):
        return [newChromosome, newChromosome]
    else:
        return [chromosome1, chromosome2]
    
def isValid(chromosome):
    baseChromosome = list(range(1, len(chromosome)))
    return all(elem in chromosome  for elem in baseChromosome)

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
    
def listCrossover(chromosomeList, fitnessList, amountCrossover):
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
        
#    print(len(chromosomeList))
#    print(len(fitnessList))
    for c in range(len(chromosomeList)-1):
#        print("xxxxxxxxxxxxxxx")
#        print((len(fitnessList)/fitnessList[c])/totDist)
#        print("xxxxxxxxxxxx")
        hitProb = 1-(1/fitnessList[c])/totDist
        randomNumber = randrange(1000)/1000
        
        if(hitProb > randomNumber):
            sMax = len(chromosomeList)
            s = 0
            found = 0
            while(found == 0 and s < sMax):
                s = s+1
                newChromosomes = crossover(shiftChromosome(chromosomeList[c], s), chromosomeList[c+1])
                if(newChromosomes[1] != chromosomeList[c+1]):
                    crossList = crossList + newChromosomes
                    found = 1
#            crossList.append(crossover(chromosomeList[c+1], chromosomeList[c]))
        crossList.append(chromosomeList[c])
#            crossList.append(chromosomeList[c+1])
    crossList.append(chromosomeList[len(chromosomeList)-1])
    return crossList

def shiftChromosome(chromosome, s):
    return numpy.roll(chromosome, s)

def mutateGroup(chromosomeList, fitnessList, amount, multiple):
    
    totDist = 0;
    mutateList = []
    
    for fitness in fitnessList:
        totDist += 1/fitness
        
    for c in range(len(chromosomeList)):
#        print("======")
#        print(amount/fitnessList[c])
#        print(randrange(1, 500)/1000*totDist/len(fitnessList))
#        print("=====")
        hitProb = amount*(1/fitnessList[c])/totDist
        randomNumber = randrange(1000)/1000
        
        if(hitProb > randomNumber):
#            print("!!!")
            if(chromosomeList[c]!=extractBestN(chromosomeList, fitnessList, 1)[0]):
                if(multiple == False):
                    mutateList.append(mutation(chromosomeList[c]))
                else:
                    mutateList.append(multipleMutation(chromosomeList[c], amount))
        
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
    
  
def epidemy(chromosomeList, parsedDataset):
    
#    for l in range(0, len(chromosomeList)-1):
#        s = 0
#        for gene in chromosomeList[l]:
#            if(gene == 1):
#                chromosomeList[l] = numpy.roll(chromosomeList[s], s)
#            s = s+1
        
    if(len(chromosomeList)>500):
        chromosomeList = killNWeakest(chromosomeList, parsedDataset, len(chromosomeList) - 250)
    return chromosomeList

def distance(chromosome, ds, drawLine):
    dist = 0;
    currentIndex = -1
    previousIndex = -1
    
    minLat = float("inf")
    maxLat = -float("inf")
    minLon = float("inf")
    maxLon = -float("inf")
    
    centerLat = 0
    centerLon = 0
    
    prevLat = []
    prevLon = []
    latC = []
    lonC = []
    
    if(drawLine == True):
        screen.fill((0,0,0))
        pygame.display.update()
            
    for c in range(len(chromosome)):
        if(c!=0):
            currentIndex = int(chromosome[c])
            previousIndex = int(chromosome[c-1])
        if(c == 0):
            currentIndex = int(chromosome[0])
            previousIndex = chromosome[len(chromosome) - 1]
        
        latC.append(float(ds[currentIndex-1].get('lat')))
        lonC.append(float(ds[currentIndex-1].get('lon')))
        
        prevLat.append(float(ds[previousIndex-1].get('lat')))
        prevLon.append(float(ds[previousIndex-1].get('lon')))
        
        dist += math.sqrt(math.pow((latC[-1]-prevLat[-1]), 2)+math.pow((lonC[-1]-prevLon[-1]), 2))
    
    if(drawLine == True):
        for x in range(0, len(latC)-1):
            if(max(latC[x], prevLat[x]) > maxLat):
                maxLat = max(latC[x], prevLat[x])
                
            if(min(latC[x], prevLat[x]) < minLat):
                minLat = min(latC[x], prevLat[x])
                
            if(max(lonC[x], prevLon[x]) > maxLon):
                maxLon = max(lonC[x], prevLon[x])
                
            if(min(lonC[x], prevLon[x]) < minLon):
                minLon = min(lonC[x], prevLon[x])
                
        centerLat = (minLat+maxLat)/2
        centerLon = (minLon+maxLon)/2
        
        scaleLat = 450/(maxLat - minLat)
        scaleLon = 450/(maxLon - minLon)
        
        trasLat = 250
        trasLon = 250
        
        for k in range(0, len(latC)):
            latC[k] = -latC[k] + centerLat
            prevLat[k] = -prevLat[k] + centerLat
            lonC[k] = -lonC[k] + centerLon
            prevLon[k] = -prevLon[k] + centerLon
            
            latC[k] = latC[k]*scaleLat + trasLat
            prevLat[k] = prevLat[k]*scaleLat + trasLat
            lonC[k] = lonC[k]*scaleLon + trasLon
            prevLon[k] = prevLon[k]*scaleLon + trasLon
        
            pygame.draw.line(screen, (255,0,0), (prevLon[k], prevLat[k]), (lonC[k], latC[k]), 1)

    pygame.display.flip()

#        pygame.display.update()

    return dist
    
amountCrossover = 0
# Upload a file
Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
filename = askopenfilename() # show an "Open" dialog box and return the path to the selected file

dataset = open(filename, "r")
parsedDataset = parseDataset(dataset)

chromosomeList = generateChromosomes(parsedDataset, 120)
chromosomeList += generateGreedyChromosomes(parsedDataset, 30)
chromosomeList += generateAlmostGreedyChromosomes(parsedDataset, 30)


evolve = 1

bestChromosomeBefore = 0
bestChromosomeAfter = 0
probabilityMutation = 0.1
probabilityMultipleMutation = 0.01

numElite = len(chromosomeList)*0.2
cumulateSaved = 0
numDelta0 = 0
eccedent = 0
first = 1

count = 0
pygame.init()
screen = pygame.display.set_mode([500, 500])


while(evolve):
    screen.fill((0, 0, 0))

    count = count+1
        
    fitnessList = generateFitnessList(chromosomeList, parsedDataset)
    chromosomeList = mutateGroup(chromosomeList, fitnessList, probabilityMutation, False)
    fitnessList = generateFitnessList(chromosomeList, parsedDataset)
    
    bestChromosomes = extractBestN(chromosomeList, fitnessList, len(chromosomeList))
        
    bestChromosomeBefore = bestChromosomeAfter
    bestChromosomeAfter = 1/fitness(bestChromosomes[0], parsedDataset)
    delta = bestChromosomeBefore - bestChromosomeAfter

    if(delta <= 100):
        numDelta0 = numDelta0 + 1
    else:
        numDelta0 = 0

    if(delta > 50):
        numElite = min(math.floor(len(chromosomeList)*0.5), numElite+math.floor(delta/100))
        chromosomeList = killNWeakest(chromosomeList, parsedDataset, len(chromosomeList) - max(math.ceil(20000/delta), 100))

        if(numDelta0 > 3):
            probabilityMutation = probabilityMutation*(numDelta0+1)
            
            
        
    elif(delta < -500 and first == 0):
        print("*****")
        numElite = max(2, numElite-1)
        probabilityMutation = probabilityMutation/1.3
        if(delta < -1500):
            chromosomeList += generateChromosomes(parsedDataset, math.floor(-delta/100))
        chromosomeList = epidemy(chromosomeList, parsedDataset)
            
    else:
        probabilityMutation = min(pow(probabilityMutation, 0.5), 1.2)
#        probabilityMultipleMutation = min(pow(probabilityMultipleMutation, 0.5), 1.2)
        

        if(numDelta0 > 2):
            probabilityMutation = probabilityMutation*(numDelta0+1)
            chromosomeList += generateChromosomes(parsedDataset, numDelta0)
            chromosomeList += generateAlmostGreedyChromosomes(parsedDataset, numDelta0)

            if(numDelta0 > 5):
                probabilityMultipleMutation = probabilityMultipleMutation+0.01*max(3*numDelta0, 50)
                if(numDelta0 > 15):
                    probabilityMultipleMutation = probabilityMultipleMutation+0.01*max(2*numDelta0, 70)
                    if(numDelta0 > 20):
                        probabilityMultipleMutation = probabilityMultipleMutation+0.01*max(2*numDelta0, 80)
                        if(numDelta0 > 30):
                            probabilityMultipleMutation = probabilityMultipleMutation+0.01*max(2*numDelta0, 90)

                fitnessList = generateFitnessList(chromosomeList, parsedDataset)
                chromosomeList = mutateGroup(chromosomeList, fitnessList, probabilityMultipleMutation, True)


    print(str(count) + " generations (distance: " + str(distance(bestChromosomes[0], parsedDataset, True)) +" and best path: " + str(bestChromosomes[0]) + ")")
    
    fitnessList = generateFitnessList(chromosomeList, parsedDataset)
    chromosomeList = listCrossover(chromosomeList, fitnessList, numElite)  
    chromosomeList = epidemy(chromosomeList, parsedDataset)
    
    first = 0
    
pygame.quit()