# -*- coding: utf-8 -*-

import math
import urllib.request
import time
import os
#######    PDB MENU    #######
def open_files(f):
    return urllib.request.urlopen("http://www.pdb.org/pdb/files/%s.pdb" %f)

# dictionnaire pour le menu de visualisation ou d'enregistrement des données
menu =  {
            1:"Visualiser le document",
            2:"Enregistrer le document",
            3:"Récupérer la séquence protéique",
            4:"Visualiser l'analyse du fichier PDB",
            5:"Enregistrer l'analyse du fichier PDB",
            6:"Profil d'hydrophobicité",
            7:"Calcul de distance entre deux carbones alpha",
            8:"Analyser un autre fichier",
            9:"Quitter"
        }
 
# dictionnaire des actions à effectuer pour récupérer la séquence protéique
recup_sequence =   {
                        1:"Récupérer la séquence protéique au format ' FASTA '",
                        2:"Récupérer la séquence protéique au format 3 lettres",
                        3:"Retourner au menu prinicipal"
                    }
 
# dictionaire affichage/enregistrement de la sequence proteique
action_sequence =   {
                         1:"Visualiser la séquence",
                         2:"Enregistrer la séquence",
                         3:"Retourner au choix du format",
                         4:"Retourner au menu principal"
                     }
# 
#dictionnaire affichage/enregistrement des analyses
analysis_display    =   {
                            1:"Visualiser le résultat",
                            2:"Enregistrer le résultat"
                         }
# variable pour l'enregistrement des fichiers
path = "./documents"
# 
# Importer le fichier de PDB - Faire les vérifications (bon code PDB)    
def verif_num(num):
    # On vérifie que le numero pdb contienne bien 4 caractères
    while len(num) !=4:
        print("/!\/!\/!\ Erreur, le code pdb est constitue de 4 caracteres /!\/!\/!\\")
        num = input("Veuillez saisir un code pdb valide svp: ")
    # on tente la récupération du fichier pdb
    try:
        opened = urllib.request.urlopen("http://www.pdb.org/pdb/files/%s.pdb" % num.upper())
        lines = opened.readlines()
        for i in range(0, len(lines)):
            lines[i] = lines[i][:-1].decode("utf8").strip()
        return [lines,str(num.upper())]
    # si la récupération échoue, on demande a l'utilisateur de vérifier son code et de saisir de nouveau
    except:
        print("Veuillez verifier votre code pdb ou votre connexion internet")
        num = input("Veuillez saisir, de nouveau, votre code pdb: ")
        return verif_num(num)
# 
def display_menu(num):
        # on demande à l'utilisateur s'il veut visualiser ou enregistrer le fichier
    print ("------------------------------  MENU  PRINCIPAL  ------------------------------")
    print("| Fichier : %s |" % num.upper())
    print ("Que voulez-vous faire? \n")
    for n,m in menu.items():
        print("--   "+str(n)+"-",m)
    print ("--------------------------------------------------------------------------------")

def display_submenu(name,myMenu):
    print("\n-----------------------  %s  ----------------------\nQuel format souhaitez vous?\n" % name)
    for n,m in myMenu.items():
        print("--   "+str(n)+"-",m)
    print("--------------------------------------------------------------------------------")
    
def display_sub_submenu(choice):
    print("------------  %s  -------------\nQue souhaitez-vous faire?" % recup_sequence[choice])
    for n,m in action_sequence.items():
        print("--   "+str(n)+"-",m)
    print("--------------------------------------------------------------------------------")
    
# Demander à l'utilisateur de rentrer le code PBD et si il est bon, affichage du menu
def file_action():
    #demande du numéro pdb
    num = input("\nSaisissez le numero du document que vous voulez consulter : ")
    # on appelle la fonction de vérification du numéro
    datas,num = verif_num(num)
    while True:
        display_menu(num)
        choice = control_number(menu)
        # on analyse un nouveau fichier
        if choice == 8:
            return 1
        # on quitte le fichier
        elif choice == 9:
            return 0
        # on visualise le fichier
        elif choice == 1:
            for line in datas:
                print(line)
        # on enregistre le fichier
        elif choice == 2:
            text=""
            for line in datas:
                text+=line+"\n"
            save_file(text,'txt')
        # on récupère la sequence protéique
        elif choice == 3:
            data = retrieve_data(datas)
            protein_sequence(menu[choice],data,num)
        # on visualise les details des analyses (PM ...)
        elif choice == 4:
            content = record_analysis(menu[choice],datas,num)
            for c in content:
                print(c)
        # on enregistre les details des analyses
        elif choice == 5:
            content = record_analysis(menu[choice],datas,num)
            save_file(content, 'txt')
        # On enregistre le profil d hydrophobicité
        elif choice == 6:
            content = getHydrophobicite(datas)
            save_file(content, "csv")
        # On calcule la distance entre deux résiduss
        elif choice == 7:
            print(getDistance(datas,num))
        time.sleep(1)

# Fonction qui retourne la distance entre deux AA
def getDistance(datas,num):
    sequence3L,exp = retrieve_data(datas)
    sequence1L = convert_3L_to_1L(sequence3L)
    res1,res2 = verif_residus(sequence1L)
    if res1 == res2:
        distance = 0
    else:
        distance = calcul_distance(str(res1), str(res2), datas)
    return "\nProtéine "+str(num)+"\nExperience : "+str(exp)+"\nDistance entre les résidus "+str(res1)+" et "+str(res2)+" :  %.3f A\n" % distance

# Fonction qui vérifie qu on ne mette pas n'importe quoi comme numéro de résidu
def verif_residus(sequence):
    while True:
        try:
            res1 = int(input("Entrez le numéro du résidu 1 : "))
            # saisie d'un caractère autre qu'un entier
            if res1 <= 0 or res1 > len(sequence):
                print("\n   /!\/!\/!\/!\ Erreur, l'entier doit etre compris entre 1 et %s /!\/!\/!\/!\\\n" % len(sequence))
                continue
            # on retourne l'entier
            else:
                res2 = int(input("Entrez le numéro du résidu 2 : "))
                if res2 <= 0 or res2 > len(sequence):
                    print("\n   /!\/!\/!\/!\ Erreur, l'entier doit etre compris entre 1 et %s /!\/!\/!\/!\\\n" % len(sequence))
                    continue
                return res1,res2
        # si l'utilisateur ne rentre pas un entier, on affiche une erreur et on lui demande de saisir de nouveau
        except ValueError:
            print ("\n   /!\/!\/!\/!\ Erreur. Veuillez saisir un entier /!\/!\/!\/!\\\n")
            continue

# Fonction qui va nous permettre d'enregistrer les analyses (PM...)
def record_analysis(action,datas,num):
    sequence3L,exp = retrieve_data(datas)
    sequence1L = convert_3L_to_1L(sequence3L)
    frequence = getFrequenceAA(sequence1L)
    return( ["--------------  Analyse du fichier "+num+"  --------------\n",
                "Expérience : "+str(exp)+"\n",
                "Poids moléculaire : "+str(poids_moleculaire(sequence1L))+" Da\n",
                "Pourcentage de résidus hydrophobes : %.2f " % (pourcent_residus_hydrophobes(sequence1L))+"%\n",
                "Pourcentage de résidus chargés : %.2f" % (pourcent_residus_charges(sequence1L))+"%\n",
                "Charge net de la protéines : "+str(charge_nette(sequence1L))+"\n",
                "\nFréquence des Acides Aminés :\n",
                " || ALA: %.2f" %(frequence["A"]), " || ARG: %.2f" %(frequence["R"]), " || ASN: %.2f\n" %(frequence["N"]),
                " || ASP: %.2f" %(frequence["D"]), " || CYS: %.2f" %(frequence["C"]), " || GLN: %.2f\n" %(frequence["Q"]),
                " || GLU: %.2f" %(frequence["E"]), " || GLY: %.2f" %(frequence["G"]), " || HIS: %.2f\n" %(frequence["H"]),
                " || ILE: %.2f" %(frequence["I"]), " || LEU: %.2f" %(frequence["L"]), " || LYS: %.2f\n" %(frequence["K"]),
                " || MET: %.2f" %(frequence["M"]), " || PHE: %.2f" %(frequence["F"]), " || PRO: %.2f\n" %(frequence["P"]),
                " || SER: %.2f" %(frequence["S"]), " || THR: %.2f" %(frequence["T"]), " || TRP: %.2f\n" %(frequence["W"]),
                " || TYR: %.2f" %(frequence["Y"]), " || VAL: %.2f" %(frequence["V"])
            ])

# Fonction qui permet d'afficher le menu de choix du format
# pour la récupération de la séquence protéique
def protein_sequence(name,datas,num):
    exit_ = 0
    while not exit_:
        display_submenu(name,recup_sequence)
        choice = control_number(recup_sequence)
        if choice != 3:
            exit_ = printer_sub_submenu(choice,datas,num)
        else:
            exit_ = 1

# Fonction qui permet de réaliser les différentes action sur les séquences protéiques
def printer_sub_submenu(choice,datas,num):
    sequence3L = datas[0]
    exp = datas[1]
    sequence1L = convert_3L_to_1L(sequence3L)
    fasta = convert_1L_to_fasta(sequence1L, num, exp)
    while True:
        display_sub_submenu(choice)
        action = control_number(action_sequence)
        #retour au menu des formats
        if action == 3:
            break
        # retour au menu principal
        elif action == 4:
            return (1)
        # fasta et on visualise
        elif action == 1 and choice == 1:
            print ("\n"+fasta+"\n")
        #fasta et on enregistre
        elif action == 2 and choice == 1:
            save_file(fasta,'txt')
        # AAA et visualisation
        elif action == 1 and choice == 2:
            print ("\nProtéine : "+num+"\nExperience : "+exp+"\n\n"+sequence3L+"\n")
        # AAA et enregistrement
        elif action == 2 and choice == 2:
            save_file(sequence3L,'txt')

# fonction de vérification des actions choisies par l'utilisateur
def control_number(myMenu):
    while True:
        try:
            c = int(input("Entrez le numéro que vous souhaitez (entre 1 et %s): " % len(myMenu)))
            # saisie d'un caractère autre qu'un entier
            if c <= 0 or c > len(myMenu):
                print("\n   /!\/!\/!\/!\ Erreur, l'entier doit etre compris entre 1 et %s /!\/!\/!\/!\\\n" % len(myMenu))
                continue
            # on retourne l'entier
            else:
                return c
        # si l'utilisateur ne rentre pas un entier, on affiche une erreur et on lui demande de saisir de nouveau
        except ValueError:
            print ("\n   /!\/!\/!\/!\ Erreur. Veuillez saisir un entier /!\/!\/!\/!\\\n")
            continue
        
# Fonction de vérification de l'existence d'un fichier
def file_exist(s):
    while True:
        # le fichier existe déjà? On demande à l'utilisateur s'il veut l'écraser ou non
        rep = input("ATTENTION! Ce fichier existe deja, voulez-vous l'ecraser?[non]/oui: ")
        if rep.lower()=="o" or rep.lower()=="oui" or rep.lower()=="y" or rep.lower()=="yes":
            name = s
            break;
        elif rep == "" or rep.lower()=="n" or rep.lower()=="non" or rep.lower()=="no":
            print("------------------------------------------")
            print("Vous avez décidé de renommer le fichier")
            name = input("Veuillez saisir un nouveau nom de fichier: ")
            # on coupe la chaîne de caractère en fonction des points
            # ainsi, si l'utiisateur saisit "test.pic" ou "test.txt", on ne récupère que "test"
            # puis on ajoute au moment de la création du fichier l'extension .txt
            name = name.split(".")
            # on vérifie que l'utilisateur a bien saisi un nom différent du nom initial
            if name[0] == s:
                continue
            else:
                break
        # si l'utilisateur ne repond pas par un des choix ci-dessus, on lui demande de recommencer la saisie
        else:
            print("Veuillez repondre par oui ou non")
    # on retourne le nom de fichier choisi par l'utilisateur
    return name
 
# Fonction qui permet de réaliser l'enregistrement des fichiers     
def save_file(content,extension):
    try:
        name = input("Sous quel nom voulez-vous enregistrer le fichier? ")
        name = name.split(".")
        # si le répertoire n'est pas crée on le crée
        if(not os.path.isdir(path)):
            print("\nStockage des fichiers à l'endroit de l'éxecution du programme")
            print("\nCréation du repertoire 'Documents' ...")
            time.sleep(5)
            os.mkdir(path)
        # Si le nom de fichier existe déjà, on appelle la fonction file_exist()
        if(os.path.exists(path+"/"+name[0]+"."+extension)):
            name[0] = file_exist(name[0])
        if extension == "txt":
            save = open(path+"/"+name[0]+".txt",'w')
#             for ligne in content:
        # on ecrit dans le fichier
            save.write(content)
        # fermeture du fichier
            save.close()
        elif extension == "csv":
            i = 4
            save = open(path+"/"+name[0]+".csv",'w')
            save.write("Position;Mesure\n")
            for line in content:
                save.write("%i;%.3f\n" %(i, line))
                i += 1
            save.close
        # message d'avertissement à l'utilisateur
        print("Fichier enregistré!")
        time.sleep(1)
    except:
        print("Nom de fichier incorrect.\nVeuillez saisir un nom de fichier valide")
        save_file(content,extension)

##########     ANALYSES    ###########
# -*- coding: utf-8

codeAA={"ALA":"A", "CYS":"C", "ASP":"D", "GLU":"E", "PHE":"F", "GLY":"G", "HIS":"H", 
        "ILE":"I", "LYS":"K", "LEU":"L", "MET":"M", "ASN":"N", "PRO":"P", "GLN":"Q", 
        "ARG":"R", "SER":"S", "THR":"T", "VAL":"V", "TRP":"W", "TYR":"Y"}

codeA = ["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"]
def retrieve_data(datas):
    sequence3L = ""
    exp = "No Experimental Data retrieved"
    # On parcourt toutes les lignes du fichier (un peu trop violent)
    for line in datas:
        # On ne s'interesse qu'aux lignes commençant par "SEQRES"
        if "SEQRES" in line[0:6]:
            # On ne récupère que les résidus (pas les 4 premiers "mots")
            for residu in line.split()[4:]:
                if residu != "ACE":
                    sequence3L = sequence3L + "-" + residu
            #On enleve le premier tiret
        if line[0:6] == "EXPDTA":
            line = line.split()
            exp = ""
            # on récupère les chaines de caractère de la liste créée par le "split"
            for mot in line[1:]:
                exp = exp + mot + " "
            # On enlève le dernier espace
            exp = exp[:-1]
    sequence3L = sequence3L[1:]
    return [sequence3L,exp]


def convert_3L_to_1L(sequence3L):
    sequence1L = ""
    # On récupère les résudus entre les "tirets"
    for residu in sequence3L.split("-"):
        sequence1L = sequence1L + codeAA[residu]
    return sequence1L

def convert_1L_to_fasta(sequence1L,num,exp):
    # La premiere ligne contient le code PDB et la methode expérimentale
    sequenceFasta = ">"+str(num)+" "+str(exp)+"\n"
    compteur = 0
    for l in sequence1L:
        if compteur != 80:
            sequenceFasta = sequenceFasta + l
            compteur = compteur + 1
        else:
            sequenceFasta = sequenceFasta + "\n" + l
            compteur = 1
    return sequenceFasta

def getFrequenceAA(sequence1L):
    dicoFreq = {}
    for aa in sequence1L :
        if aa in codeA:
            if aa in dicoFreq :
                dicoFreq[aa] = dicoFreq[aa] + 1
            else :
                dicoFreq[aa] = 1
    for aa in dicoFreq :
        dicoFreq[aa] = dicoFreq[aa] / len(sequence1L)
    for aa in codeA:
        if aa not in sequence1L:
            dicoFreq[aa] = 0
    return(dicoFreq)

def getHydrophobicite(datas):
    sequence,exp = retrieve_data(datas)
    return (hydrophobicite(sequence.split("-")))
# Prend en entrée une liste de résidus 3 lettres
# sequence3L.split("-")
def hydrophobicite(sequence):
    somme = 0
    moyenne = 0
    liste_hydro = []
    dicoAA_hydro={"ALA":0.310, "CYS":1.540, "ASP":-0.770, "GLU":-0.640, "PHE":1.790, "GLY":0.000, "HIS":0.130, "ILE":1.800, "LYS":-0.990, "LEU":1.700, "MET":1.230, "ASN":-0.600, "PRO":0.720, "GLN":-0.220, "ARG":-1.010, "SER":-0.040, "THR":0.260, "VAL":1.220, "TRP":2.250, "TYR":0.960}
    i=0
    # On parcourt toute la séquence
    while i<len(sequence)-9:
        somme = 0
        j=0
        # On prend
        while j<9:
            somme = somme + dicoAA_hydro[sequence[i+j]]
            j = j + 1
        moyenne = round(somme/9,4)
        liste_hydro.append(moyenne)
        i = i + 1
    return liste_hydro

def calcul_distance(res1,res2,datas):
    distance = 0
    [x1,x2,y1,y2,z1,z2] = [0,0,0,0,0,0]
    for ligne in datas:
        ligne = ligne.split()
        if ligne[0].upper() == "ATOM":
            if ligne[5] == res1 and ligne[2] == "CA":
                x1 = float(ligne[6])
                y1 = float(ligne[7])
                z1 = float(ligne[8])
            elif ligne[5] == res2 and ligne[2] == "CA":
                x2 = float(ligne[6])
                y2 = float(ligne[7])
                z2 = float(ligne[8])

    distance = math.sqrt( (x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2 )
    distance = round(distance,3)

    return(distance)

# Calcul du poids moleculaire de la structure
poids_moleculaireAA = {"A":89 , "R":174 , "N":132 , "D":133 , "C":121 , "Q":146 , "E":147 , "G":75 , "H":155 , "I":131 , "L":131 , "K":146 , "M":149 , "F":165 , "P":115 , "S":105 , "T":119 , "W":204 , "Y":181 , "V":117 }
    
def poids_moleculaire(sequence1L):
    l=0
    for AA in sequence1L:
        l=l+poids_moleculaireAA[AA]
    return (l)

# Calcul du pourcentage de résidus hydrophobes
acides_amines_hydrophobes = ["M", "A", "V", "L", "I", "P", "F", "W"]

def pourcent_residus_hydrophobes (sequence1L):
    l=0
    for AA in sequence1L:
        if AA in acides_amines_hydrophobes:
            l=l+1
    return(l / len(sequence1L) * 100)


# Calcul du pourcentage de residus chargés
acides_amines_charges = ["D", "E", "R", "K", "H"]

def pourcent_residus_charges (sequence1L):
    l=0
    for AA in sequence1L:
        if AA in acides_amines_charges:
            l=l+1
    return(l / len(sequence1L) * 100)


# Calcul de la charge nette
acides_amines_positifs = ["R", "K", "H"]
acides_amines_negatifs = ["D", "E"]

def charge_nette(sequence1L):
    l=0
    for AA in sequence1L:
        if AA in acides_amines_positifs:
            l=l+1
        elif AA in acides_amines_negatifs:
            l=l-1
    return(l)
    
##############    PROGRAMMES   ##########
print ("-------------------------------")
print ("|         PDB analysis        |")
print ("-------------------------------")
analysis = 1
while analysis:
    #     Affichage du menu
    analysis = file_action()
input("Appuyer sur entrer pour quitter")