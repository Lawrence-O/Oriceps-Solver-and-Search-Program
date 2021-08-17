
#*************
# Created by: Lawrence Onyango
#*************


#Imports
from sympy import*
import concurrent.futures
import time
import csv
import sys

#Creates a CSV File within the designated data folder in exports
def createCSV(a_initial,a_final):
    with open(f"data_{a_initial}_{a_final}.csv", "a", newline="") as file:
        csv.writer(file, delimiter = ";", quotechar = '"', quoting=csv.QUOTE_MINIMAL).writerow(["h0","h1","delH","a","b","c","c_constraint"])

#Inserts the corresponding values(h0,h1,delH,a,b,c) into a created CSV File in the data folder
def insertValues_CSV(h0,h1,delH,a,b,c,c_constraint,a_initial,a_final):
    with open(f"data_{a_initial}_{a_final}.csv", "a", newline="") as file:
        csv.writer(file, delimiter = ";", quotechar = '"', quoting=csv.QUOTE_MINIMAL).writerow([h0,h1,delH,a,b,c,c_constraint])

#Returns a sympy equation with h as an independent variable
def equation(a,b,c,h,J):
    return Eq((((8*b**3) + (8*b*c**2) - (2*b*h**2) + (2*a*(4*b**2 - 4*c**2 + h**2)))/(4*b**2 + 4*c**2 - h**2)),J)

#Solves the sympy equation
def solveEquation(a,b,c,a_initial,a_final):
    h0,h1 = symbols('h0 h1', real = True)
    eq0 = equation(a,b,c,h0,3.5)
    eq1 = equation(a,b,c,h1,7.5)
    #Used to capture errors (mostly here due to sympy's finicky ness)
    try:
        sol = list(nonlinsolve([eq0,eq1],[h0,h1]))
        parsedSol(sol,a,b,c,a_initial,a_final)
    except:
        pass

#Parses the solutions to get rid of negative h0,h1 and delHs
def parsedSol(sol,a,b,c,a_initial,a_final):
    for h0,h1 in sol: 
        if(h0 > 0  and h1 > 0 and h1-h0 > 0): 
            #Determines whether the values abide by the c-constraint filter
            if(0 < h0 < 2*c and 0 < h1 < 2*c):
                c_constraint = 1
            else:
                c_constraint = 0
            insertValues_CSV(h0,h1,h1-h0,a,b,c,c_constraint,a_initial,a_final) 

#This does the calculation via concurrent futures(multi-processing via CPU cores)
def doCalculations(a_initial,a_final):
    initialTime = time.time()
    print("Starting...")
    with concurrent.futures.ProcessPoolExecutor() as executor: #ProcessPool is used since solving is CPU bound not I/O
        for a in range(a_initial,a_final):
            for nL in range(10,294):#Max of 1472 
                for b in range(10,nL):#Max of 1472
                    c = ((nL*0.34)**2 - (b*0.34)**2)**(0.5) 
                    executor.submit(solveEquation,round(a*0.34,5),round(b*0.34,5),c,a_initial,a_final)
    print("Done in: ",time.time()-initialTime)

#Main Function where the user interacts with the program
def main():
    a_initial = input("Insert Starting A Value: ")
    a_final = input("Insert Final A Value: ")
    createCSV(int(a_initial),int(a_final))
    doCalculations(int(a_initial), int(a_final))
    isContinue = input("Do You Wish to Continue ?")
    if(isContinue.lower() in ["y", "yes", "t", "true"]):
        main()
    else:
        sys.exit()


if __name__ == '__main__':
    main()




