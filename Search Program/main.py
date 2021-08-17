#*************
# Created by: Lawrence Onyango
#*************

#Imports
from sympy import*
import matplotlib.pyplot as plt
import pandas as pd
import sys
import os
import drawSvg as draw

#Solves for the scaffold length       
def scaffoldLengthCalc(a,b,c):
    return (4*((b**2 + c**2)**0.5) + 4*(((a+b)**2 + c**2)**0.5)+ 2*a + 2*b)/0.34

#Searches for corresponding values via Pandas
def searchValues(h0,h1,rang,c_constraintFilter):
    print("Using CSV....")
    df = pd.read_csv('data.csv',delimiter=';')
    #Adds the Scaffold Length column to the data
    df["scaffoldLength(bases)"] = scaffoldLengthCalc(df["a"],df["b"],df["c"])
    if(c_constraintFilter):
        return df[df['h0'].between(h0-rang,h0+rang) & df['h1'].between(h1-rang,h1+rang) & df['c_constraint'] == 1]
    return df[df['h0'].between(h0-rang,h0+rang) & df['h1'].between(h1-rang,h1+rang)]
    
#Plots a AvBvC plot with a Delta H color bar
def threeDPlot(h0,h1,rang,df):
    fig = plt.figure(figsize=(10,10), dpi = 100)
    ax = plt.axes(projection="3d")
    ax.set_xlabel('A')
    ax.set_ylabel('B')
    ax.set_zlabel('C')
    ax.set_title("AvBvC(ColorBar = Delta H)")
    p = ax.scatter3D(df['a'],df['b'],df['c'],c=df['deltaH'],cmap="magma")
    fig.colorbar(p)
    plt.savefig(f"Exports/Data_{h0}_{h1}_{rang}/AvBvC(ColorBar = Delta H)")
    plt.show()

#Plots a h0vh1 graph with a colorbar of deltaH
def h0vh1(h0,h1,rang,df):
    p1 = plt.scatter(df['h0'],df['h1'],c = df['deltaH'], cmap = "magma")
    plt.colorbar(p1)
    plt.title('H0vH1(ColorBar = Delta H)')
    plt.ylabel('H1')
    plt.xlabel('H0')
    plt.savefig(f"Exports/Data_{h0}_{h1}_{rang}/H0vH1(ColorBar = Delta H)")
    plt.show()

#Plots h0 on the x and the difference between the target and actual value on the Y
def h0vTargetDifference(h0,h1,rang,df):
    df["rangeDifferenceH0"] = abs(df['h0'] - float(h0))
    plt.scatter(df['h0'],df['rangeDifferenceH0'])
    plt.title('H0vTargetDifference(h0 - desiredH0)')
    plt.ylabel('Difference')
    plt.xlabel('H0')
    plt.savefig(f"Exports/Data_{h0}_{h1}_{rang}/H0vTargetDifference(h0 - desiredH0)")
    plt.show()

#Plots h1 on the x and the difference between the target and actual value on the Y
def h1vTargetDifference(h0,h1,rang,df):
    df["rangeDifferenceH1"] = abs(df['h1'] - float(h1))
    plt.scatter(df['h1'],df['rangeDifferenceH1'])
    plt.title('H1vTargetDifference(h1 - desiredH1)')
    plt.ylabel('Difference')
    plt.xlabel('H1')
    plt.savefig(f"Exports/Data_{h0}_{h1}_{rang}/H1vTargetDifference(h1 - desiredH1)")
    plt.show()

#Plots an overall graph of h0vh1 if not present within the exports folder
def overallGraph(h0,h1,rang):
    if(not os.path.isfile("Exports/H0vDeltaH(Overall Map).png")):
        df = pd.read_csv('data.csv',delimiter=';')
        plt.scatter(df['h0'],df['deltaH'])
        plt.title('H0vDeltaH(Overall Map)')
        plt.ylabel('Delta H')
        plt.xlabel('H0')
        plt.savefig(f"Exports/H0vDeltaH(Overall Map)")
        plt.show()
    if(not os.path.isfile("Exports/H0vDeltaH(Overall Map Limited to 100).png")):
        df = pd.read_csv('data.csv',delimiter=';')
        df = df[df["h0"].between(0,100)]
        plt.scatter(df["h0"],df['deltaH'])
        plt.title('H0vDeltaH(Overall Map Limited to 100)')
        plt.ylabel('Delta H')
        plt.xlabel('H0')
        plt.savefig(f"Exports/H0vDeltaH(Overall Map Limited to 100)")
        plt.show()

#Exports the chosen values to a csv file within the corresponding data folder in exports
def exportCSV(h0,h1,rang,df):
    if(not os.path.exists(f"Exports/Data_{h0}_{h1}_{rang}/Data_{h0}_{h1}_{rang}.csv")):
        os.makedirs(f"Exports/Data_{h0}_{h1}_{rang}",exist_ok= True)
        df.to_csv(f"Exports/Data_{h0}_{h1}_{rang}/Data_{h0}_{h1}_{rang}.csv",index=False)

#Calls all the graphing functions
#Note all graphs are exported to their corresponding data folder within exports
def graph(h0,h1,rang,df):
    h0vh1(h0,h1,rang,df)
    h0vTargetDifference(h0,h1,rang,df)
    h1vTargetDifference(h0,h1,rang,df)
    threeDPlot(h0,h1,rang,df)
    overallGraph(h0,h1,rang)

#Returns a sympy function for J
def equation(a,b,c,h,J):
    return Eq((((8*b**3) + (8*b*c**2) - (2*b*h**2) + (2*a*(4*b**2 - 4*c**2 + h**2)))/(4*b**2 + 4*c**2 - h**2)),J)

#Solves the sympy equation
#Used to find a h0/h1 for the oriceps creation
def findH0H1(a,b,c):
    h0,h1 = symbols('h0 h1', real = True)
    eq0 = equation(a,b,c,h0,3.5)
    eq1 = equation(a,b,c,h1,7.5)
    try:
        sol = list(nonlinsolve([eq0,eq1],[h0,h1]))
        return parsedSol(sol,a,b,c)
    except:
        return None,None

def parsedSol(sol):
    for h0,h1 in sol: 
        if(h0 > 0  and h1 > 0 and h1-h0 > 0): 
            return h0,h1
    return None,None

#Uses drawSVG to draw the oriceps and then exports a png and svg file
def drawForeceps(a,b,c,scale,gap,width,height,h0,h1):
    if(not os.path.exists(f"Exports/Oriceps")):
        os.makedirs(f"Exports/Oriceps",exist_ok= True)
    d = draw.Drawing(width,height, origin = "center", displayInline= False)
    background = draw.Rectangle(-width//2,height//2,width,-height, fill='#FFFFFF')
    d.append(background)
    legendBox = draw.Rectangle(c*scale//2 + gap,-b*scale - a*scale//2,width//2,-height//2, fill='#FFFFFF', stroke = "Black", stroke_width = 12)
    d.append(legendBox)
    d.append(draw.Text('Legend:',width//30, c*scale//2 + width//30 + gap, -b*scale - a*scale//2 - width//30, fill='Black'))
    d.append(draw.Text(f'A:{a}', width//30, c*scale//2 + width//30 + gap, -b*scale - a*scale//2 - 2*width//30, fill='Black'))
    d.append(draw.Text(f'B:{b}', width//30, c*scale//2 + width//30 + gap, -b*scale - a*scale//2 - 3*width//30, fill='Black'))
    d.append(draw.Text(f'C:{c}', width//30, c*scale//2 + width//30 + gap, -b*scale - a*scale//2 - 4*width//30, fill='Black'))
    if(h0 == None or h1 == None):
        d.append(draw.Text(f'H0:N/A', width//30, c*scale//2 + width//30 + gap, -b*scale - a*scale//2 - 5*width//30, fill='Black'))
        d.append(draw.Text(f'H1:N/A', width//30, c*scale//2 + width//30 + gap,-b*scale - a*scale//2 - 6*width//30, fill='Black'))
    else:
        d.append(draw.Text(f'H0:{round(h0,3)}', width//30,c*scale//2 + width//30 + gap, -b*scale - a*scale//2 - 10*width//30, fill='Black'))
        d.append(draw.Text(f'H1:{round(h1,3)}', width//30, c*scale//2 + width//30 + gap,-b*scale - a*scale//2 - 12*width//30, fill='Black'))
    d.append(draw.Line(0,0,0,b*scale,fill = '#000000',stroke = "Black", stroke_width = 12)) # B ***
    d.append(draw.Line(gap,b*scale+gap//4,c*scale+gap,gap//4,fill = '#000000',stroke = "Black",stroke_width = 12)) #nL
    d.append(draw.Line(-gap,b*scale+gap//4,-c*scale-gap,gap//4,fill = '#000000',stroke = "Black", stroke_width = 12)) #nL
    d.append(draw.Line(0,b*scale+gap,0,b*scale+a*scale+gap,fill = '#000000',stroke = "Black", stroke_width = 12)) #a
    d.append(draw.Line(gap,b*scale+a*scale+gap,c*scale+gap,gap,fill = '#000000',stroke = "Black", stroke_width = 12)) #K
    d.append(draw.Line(-gap,b*scale+a*scale+gap,-c*scale-gap,gap,fill = '#000000',stroke = "Black", stroke_width = 12)) #K
    d.append(draw.Line(0,0,0,-b*scale,fill = '#000000',stroke = "Black", stroke_width = 12)) # B ***
    d.append(draw.Line(gap,-b*scale-gap//4,c*scale+gap,-gap//4,fill = '#000000',stroke = "Black", stroke_width = 12)) #nL
    d.append(draw.Line(-gap,-b*scale-gap//4,-c*scale-gap,-gap//4,fill = '#000000',stroke = "Black", stroke_width = 12)) #nL
    d.append(draw.Line(0,-b*scale-gap,0,-b*scale-a*scale-gap,fill = '#000000',stroke = "Black", stroke_width = 12)) #a
    d.append(draw.Line(gap,-b*scale-a*scale-gap,c*scale+gap,-gap,fill = '#000000',stroke = "Black", stroke_width = 12)) #K
    d.append(draw.Line(-gap,-b*scale-a*scale-gap,-c*scale-gap,-gap,fill = '#000000',stroke = "Black", stroke_width = 12))
    d.saveSvg(f'Exports/Oriceps/Oriceps_{a}_{b}_{c}.svg')
    d.savePng(f'Exports/Oriceps/Oriceps_{a}_{b}_{c}.png')

#The mainMenu where the user interacts with the program
def mainMenu():
    print("Welcome to The Oriceps Program")
    print("Press 1 to run the search program.")
    print("Press 2 to graph the oriceps.")
    print("Press 3 exit.")
    selection = input("Selection: ")
    if(selection in ["1", " 1", "1 "]):
        runSearchProgram()
    elif(selection in ["2", " 2", "2 "]):
        oricepsAsker()
    elif(selection in ["3", " 3", "3 "]):
        sys.exit()
    else:
        mainMenu()

#Returns the user input for the search program
def runAsker():
    print("Welcome To The Search Program.")
    h0 = input("Insert H0: ")
    h1 = input("Insert H1: ")
    rang = input("Insert Value Range Above Number: ")
    c_constraintFilter = input("C-Constraint Filter ? [y/n]")
    if(c_constraintFilter.lower() in ["y","yes","t","true"]):
        c_constraintFilter = True
    elif c_constraintFilter.lower() in ["n","no","f","false"]:
        c_constraintFilter = False
    else:
        runAsker()
    return float(h0),float(h1),float(rang),c_constraintFilter

#Handles the search program
def runSearchProgram():
    h0,h1,rang,c_constrantFilter = runAsker()
    df = searchValues(float(h0),float(h1),float(rang),c_constrantFilter)
    exportCSV(h0,h1,rang,df)
    graph(h0,h1,rang,df)
    mainMenu()

#Handles the oricep graphing functionality
def oricepsAsker():
    a = float(input("Insert A: "))
    b = float(input("Insert B: "))
    c = float(input("Insert C: "))
    scale = input("Insert A Scale Multiplier[x for automatic scale]: ")
    gap = float(input("Insert The Gap Size Between Bundles: "))
    width = float(input("Width of Paper(px): "))
    height = float(input("Height of Paper(px): "))
    isCalculate = input("Do you wish to calculate a corresponding(if any) h0 and h1 value? [y/n] ")
    if(scale.lower() == "x"):
        scaleH = height//(2*b + 2*a + gap//2)
        scaleW = width//(2*c + gap//2)
        scaleNew = min(scaleH,scaleW)
    else:
        scaleNew = float(scale)
    if(isCalculate.lower() in ["y","yes","t","true"]):
        h0,h1 = findH0H1(a,b,c)
        drawForeceps(a,b,c,scaleNew,gap,width,height,h0,h1)
    else:
        drawForeceps(a,b,c,scaleNew,gap,width,height,None,None) 
    mainMenu()

#Calls the Main Menu
def main():
    mainMenu()

'''
helixies
'''
if __name__ == '__main__':
    main()