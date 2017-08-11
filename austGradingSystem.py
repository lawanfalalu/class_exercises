# -*- coding: utf-8 -*-
"""
Created on Sat May 27 17:12:48 2017

@author: FALALU
"""
###################################### AUST GRADING SYSTEM #############################################
def GradePoints(mark):
    gSys = {"A":4.0,"A-":3.75,"B+":3.25,"B":3.0,"B-":2.75,"C+":2.25,"C":2.0,"C-":1.75,"D":1.0,"F":0.0,}
#    print(gSys[mark])
    return gSys[mark]
def GradeMe(mks):
    if mks>=95 and mks<=100:
        return "A"
    elif mks>=89:
        return "A-"
    elif mks>=83:
        return "B+"
    elif mks>=77:
        return "B"
    elif mks>=71:
        return "B-"
    elif mks>=65:
        return "C+"
    elif mks>=59:
        return "C"
    elif mks>=53:
        return "C-"
    elif mks>=48:
        return "D"
    else:
        return "F"
   
def GetGrade():
    trans={}
    points=0
#    credit= float(input("WHAT'S THE CREDIT UNIT ? :"))
#    print("\n")
    total= int(input("HOW MANY COURSES ? :"))
    if total > 0:
        for c in range(total):
            code = input("COURSES CODE %d ? :"%(c+1))
            mark = float(input("MARKS OPTAINED IN \"%s\" ? :"%(code)))
            trans[code]=mark

        for key,value in trans.items():
            print(key,value,GradeMe(value))
            points = points + GradePoints(GradeMe(value))
            print("\tCGPA = #",round((points/total),3))
    else:
        print("CANNOT CALCULATE %d COURSES !"%(total))
######################################END OF AUST GRADING SYSTEM #############################################

############################  MAIN ENTRY POINT ###############################
#complete menu
def startMenu():
	while 1==1:
		print("<<<< Module for solving Equations >>>>")
		print("Press 1 for Perfect EGG solutions :")
		print("Press 2 for Simultaneous Equations :")
		print("Press 3 for Quaratic Equations :")
		print("Press 4 for Bisection Method :")
		print("Press 5 for NewtonRaphson Method :")		
		print("Press 6 for Gause-Seidel :")
		print("Press 7 for Gause-Jacobi :")
		print("Press 8 for LU Decomposition :")
		print("Press 9 for AUST CGPA Computations :")
		print("Press 0 TO QUIT :")
		CHOICE=int(input("PLESE ENTER YOUR CHOICE :"))
		if CHOICE == 0:
			break
		elif CHOICE == 1:
			print("calling",CHOICE)
		elif CHOICE == 2:
			print("calling",CHOICE)
		elif CHOICE == 3:
			print("calling",CHOICE)
		elif CHOICE == 4:
			print("calling",CHOICE)
		elif CHOICE == 5:
			print("calling",CHOICE)
		elif CHOICE == 6:
			print("calling",CHOICE)
		elif CHOICE == 7:
			print("calling",CHOICE)
		elif CHOICE == 8:
			print("calling",CHOICE)
		elif CHOICE == 9:
			print("calling",CHOICE)
			GetGrade()  #ENTRY POINTS 
		else:
			print("\n",CHOICE," is INVALID CHOICE!!!\nPLEASE ENTER 0...9 Only\n")
	print("\n<<<<<<<<<<<<----------GOOD BYE---------->>>>>>>>>>>>\n")
startMenu()
#####################################END OF MAIN ENTRY POINT#########################################
            
