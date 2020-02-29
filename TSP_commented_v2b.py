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
import copy

showGraphicalOutput = True

def parseDataset(ds): # Parsing of the given dataset
    csvReader = csv.reader(ds, delimiter=' ')
    parsedDataset = []
    for row in csvReader:
        parsedDataset.append({'num': row[0], 'lat': row[1], 'lon': row[2]})
    
    return parsedDataset

def generateChromosomes(ds, number):
    # Generates a list of random chromosomes and returns it
    listChromosomes = []
    
    for c in range (number):
        randomChromosome = newChromosome(ds)
        listChromosomes.append(randomChromosome)
    
    return listChromosomes

def generateGreedyChromosomes(ds, number):
    # Gnerates a list of greedy chromosomes and returns it
    listChromosomes = []
    
    number = min(number, len(ds))
    for k in range(1, number):
        index = randrange(1, number)
        city = ds[index-1]
        listChromosomes.append(generateGreedyChromosome(city.get('num'), ds))
    
    return listChromosomes

def generateGreedyChromosome(startingPoint, ds):
    # Gnerates a single greedy chromosome
    # The idea is to have a base chromosome with all the cities' ids
    # Every time a city is chosen, it is removed from the base chromosome and put in the greedy chromosome
    
    numberCities = len(ds)+1
    baseChromosome = list(range(1, numberCities)) # this generates a list [1, 2, ... n]
    greedyChromosome = [] # the greedy chromosome is initially empty
    
    greedyChromosome.append(int(startingPoint)) # add the starting point to the greedy chromosome
    baseChromosome.remove(int(startingPoint)) # and remove it from the base chromosome
    
    for g in range(numberCities-2):
        closestPoint = findClosest(greedyChromosome[-1], baseChromosome, ds)[0] # finds the closest unexplored city
        greedyChromosome.append(closestPoint) # adds it to the greedy chromosome
      
        baseChromosome.remove(int(closestPoint)) # removes it from the base chromosome (-> don't consider it anymore)
        
    return greedyChromosome


def generateAlmostGreedyChromosomes(ds, number):
    # SImilar to the greedy chromosomes generation
    listChromosomes = []
    
    for l in range(0, max(1, math.ceil(number/len(ds)))):
        city = ds[randrange(0, len(ds)-1)]
        listChromosomes.append(generateAlmostGreedyChromosome(city.get('num'), ds))
    
    return listChromosomes

def generateAlmostGreedyChromosome(startingPoint, ds):
    # The only difference with respect to the greedy chromosome generation consists in choosing the next city
    # among the closest, so the ciy is not necessarily the closest
    numberCities = len(ds)+1
    baseChromosome = list(range(1, numberCities))
    greedyChromosome = []
    closestPoint = -1
    
    greedyChromosome.append(int(startingPoint))
    baseChromosome.remove(int(startingPoint))
    
    for g in range(numberCities-2):
        randomValue = random.random()
        scaledRandom = 1
        if(randomValue > 0.9):
            scaledRandom = 1
        else:
            scaledRandom = 2
                
        firstClosest = findClosest(greedyChromosome[-1], baseChromosome, ds)[0]
        secondClosest = findClosest(firstClosest, baseChromosome, ds)[0]
        
        if(scaledRandom == 1 or len(greedyChromosome) == 1):
            closestPoint = firstClosest
        else:
            closestPoint = secondClosest
            
        greedyChromosome.append(closestPoint)
      
        baseChromosome.remove(int(closestPoint))
        
    return greedyChromosome


def findClosest(element, baseChromosome, ds):
    # Given a city and the list of remaining cities, find the closest unexplored city
    minDist = float("inf")
    minElem = 0
    currentLat = -1
    currentLon = -1
        
    for el in ds:
        if(int(el.get('num'))==int(element)):
            currentLat = float(el.get('lat'))
            currentLon = float(el.get('lon'))
            
    for elB in baseChromosome: # finds the minimum in the base chromosome
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
    

def newChromosome(ds): # generates a new chromosome randomly
    chromosomeLength = len(ds)+1
    baseChromosome = list(range(1, chromosomeLength))
    for l in range(chromosomeLength):
        index1 = randrange(chromosomeLength-1)
        index2 = randrange(chromosomeLength-1)
        temp = baseChromosome[index2]
        baseChromosome[index2] = baseChromosome[index1]
        baseChromosome[index1] = temp
    return baseChromosome

def fitness(chromosome, ds): # the fitness function takes the distance and calculates the inverse
    return 1/distance(chromosome, ds, False)
    
def mutation(refChromosome): # swap of two cities inside a chromosome (chromosome = path)
    chromosome = refChromosome.copy()
    chromosomeLength = len(chromosome)
    index1 = randrange(chromosomeLength-1)
    index2 = randrange(chromosomeLength-1)
    temp = chromosome[index2]
    chromosome[index2] = chromosome[index1]
    chromosome[index1] = temp
    return chromosome

def multipleMutation(chromosome, percentage): # performs more than one mutation at the same time
    multipleMutatedChromosome = chromosome.copy()
    for m in range(0, max(math.floor(min(percentage, 1)*len(chromosome)-1), 1)):
        multipleMutatedChromosome = mutation(multipleMutatedChromosome)
    
    return multipleMutatedChromosome
        
def crossover(parent1, parent2): # crossover of two chromosomes (for each position, takes one element from one chromosome)
    
    firstList = copy.deepcopy(parent1)
    secondList = copy.deepcopy(parent2)
    
    crossedList = []
    addedCities = []
    
    for i in range(0, len(firstList)):
        if((firstList[i] not in addedCities) and (secondList[i] not in addedCities)):
            casualElement = randrange(0,1)
            if(casualElement == 0):
                crossedList.append(firstList[i])
            else:
                crossedList.append(secondList[i])
        elif((firstList[i] not in addedCities) or (secondList[i] not in addedCities)):
            if(firstList[i] not in addedCities):
                crossedList.append(firstList[i])
            else:
                crossedList.append(secondList[i])
        else:
            j = 0
            while j<len(firstList):
                if(firstList[j] not in addedCities):
                    crossedList.append(firstList[j])
                    addedCities.add(firstList[j])
                    j = len(firstList)
                elif(secondList[j] not in addedCities):
                    crossedList.append(secondList[j])
                    addedCities.add(secondList[j])
                    j = len(firstList)
                else:
                    j = j+1

    return [crossedList]


    
def isValid(chromosome): # the chromosome contains a list of different cities
    baseChromosome = list(range(1, len(chromosome)))
    return all(elem in chromosome  for elem in baseChromosome)

def generateFitnessList(chromosomeList, parsedDataset): # this function would not be necessary in a real OOP implementatio. It calculates the fitness for all the points.
    fitnessList = []
    for f in range(len(chromosomeList)):
        fitnessList.append(fitness(chromosomeList[f], parsedDataset))
    return fitnessList
        
def extractBestN(chromosomeList, fitnessList, n): # Sorts the chromosomes from the best to the worst
    bestChromosomesList = []
    indexList = []
        
    highestList = sorted(fitnessList, reverse=True)[:n]
    
    for i in range(0, len(highestList)):
        indexList.append(fitnessList.index(highestList[i]))
            
    for j in indexList:
        bestChromosomesList.append(chromosomeList[j])
    
    return bestChromosomesList
        
    
def listCrossover(chromosomeList, fitnessList, amountCrossover):
    # Performs multiple corssovers
    # The probability of a crossover depends on the crossover amount
    
    totDist = 0;
    crossList = []
    
    for fitness in fitnessList:
        totDist += 1/fitness
        
    for c in range(len(chromosomeList)-1):

        hitProb = 1-(1/fitnessList[c])/totDist
        hitProb = hitProb*amountCrossover
        
        randomNumber = randrange(1000)/1000
        
        if(hitProb > randomNumber):
            sMax = len(chromosomeList)
            s = 0
            found = 0
            while(found == 0 and s < sMax):
                s = s+1
                newChromosomes = crossover(shiftChromosome(chromosomeList[c], s), chromosomeList[c+1])
                
                if(newChromosomes[0] != chromosomeList[c+1]):
                    crossList = crossList + newChromosomes
                    found = 1
                    
    return crossList

def shiftChromosome(refChromosome, s): # to define the fact that a path c1, c2, ..., cn is the same as cn, c1, c2 ... c(n-1)
    chromosome = refChromosome.copy()
    return numpy.roll(chromosome, s)

def mutateGroup(chromosomeList, fitnessList, amount, multiple):
    # Performs multiple mutations
    
    totDist = 0;
    mutateList = []
    
    for fitness in fitnessList:
        totDist += 1/fitness
          
    for c in range(0, len(chromosomeList)):

        mutateList.append(chromosomeList[c])

        hitProb = amount*fitnessList[c]*totDist
        
        randomNumber = randrange(1000)/1000
        
        if(hitProb > randomNumber or 1/fitnessList[c] > 5*totDist/len(chromosomeList)):
            if(multiple == False):
                mutateList.append(mutation(chromosomeList[c]))
            else:
                mutateList.append(multipleMutation(chromosomeList[c], amount))
    
    return mutateList
    
def killNWeakest(chromosomeList, ds, n):
    # This allows to perform elitism
    n = max(len(chromosomeList) - n, 4)
    
    fitnessList = generateFitnessList(chromosomeList, ds)
    toBeSaved = extractBestN(chromosomeList, fitnessList, n)
    
    return toBeSaved
    

def distance(chromosome, ds, drawLine):
    # Given two points, find the distance. If drawing is requested, draws the path in pygame
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
    
    if(showGraphicalOutput == True):
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
    
    if(showGraphicalOutput == True):
        if(drawLine == True):
            for x in range(0, len(latC)):
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

    return dist
    
amountCrossover = 0
# Upload a file
Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
filename = askopenfilename() # show an "Open" dialog box and return the path to the selected file

dataset = open(filename, "r")
parsedDataset = parseDataset(dataset)

maxChromosomes = 10 #len(parsedDataset)
maxGreedy = len(parsedDataset)

print("Generating the random chromosomes...")
chromosomeList = generateChromosomes(parsedDataset, maxChromosomes)

print("Generating the greedy chromosomes...")
chromosomeList += generateGreedyChromosomes(parsedDataset, maxGreedy)

print("Generating the almost greedy chromosomes...")
#chromosomeList += generateAlmostGreedyChromosomes(parsedDataset, maxGreedy)

print(maxChromosomes)

evolve = 1

bestChromosomeBefore = 0
bestChromosomeAfter = 0
probabilityMutation = 0.3
probabilityMultipleMutation = 0.01

numElite = int(0.6*maxChromosomes)
cumulateSaved = 0
numDelta0 = 0
eccedent = 0
first = 1

count = 0

if(showGraphicalOutput == True):
    pygame.init()
    pygame.display.set_caption(filename.split("/")[-1])
    screen = pygame.display.set_mode([500, 500])

deltaThreshold = 1

print("Starting evolution...")
while(evolve):
    count = count + 1
    
    fitnessList = generateFitnessList(chromosomeList, parsedDataset)
    chromosomeList = chromosomeList + mutateGroup(chromosomeList, fitnessList, probabilityMutation, False)
    
    fitnessList = generateFitnessList(chromosomeList, parsedDataset)
    chromosomeList = chromosomeList + mutateGroup(chromosomeList, fitnessList, probabilityMultipleMutation, True)

    fitnessList = generateFitnessList(chromosomeList, parsedDataset)

    bestChromosomes = extractBestN(chromosomeList, fitnessList, numElite)
        
    bestChromosomeBefore = bestChromosomeAfter
    bestChromosomeAfter = distance(bestChromosomes[0], parsedDataset, False)
    
    delta = bestChromosomeBefore - bestChromosomeAfter

    if(delta <= deltaThreshold):
        numDelta0 = numDelta0 + 1
    else:
        numDelta0 = 0
    
    probabilityMutation = min(numDelta0/10, 1)
    probabilityMultipleMutation = min(numDelta0/20, 1)
    
    chromosomeList = killNWeakest(chromosomeList, parsedDataset, len(chromosomeList) - maxChromosomes)
    
    if(numDelta0 == 0 or count % 100 == 0):
        bestDistance = 0
        if(showGraphicalOutput == True):
            bestDistance = distance(bestChromosomes[0], parsedDataset, True)
            pygame.display.set_caption(filename.split("/")[-1] + " - distance: " + str(round(bestDistance, 2)) + " - generation " + str(count))
        else:
            bestDistance = distance(bestChromosomes[0], parsedDataset, False)
        print(str(count) + " generations (distance: " + str(bestDistance) +" and best path: " + str(bestChromosomes[0]) + ")")
    
    fitnessList = generateFitnessList(chromosomeList, parsedDataset)
        
    chromosomeList = chromosomeList + listCrossover(chromosomeList, fitnessList, numElite)  
    
    first = 0
            
pygame.quit()