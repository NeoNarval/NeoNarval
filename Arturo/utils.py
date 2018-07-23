import collections

class CfgFile():

    def __init__(self, url):
        self.__url = url
        self.__dico = collections.OrderedDict()
        with open(self.__url, 'r') as fichier:
            content = fichier.readlines()
            
        content = [line.strip() for line in content]
        for line in content:
            param = line.split(" : ")
            if len(param) == 2:
                nom_param = param[0]
                value_param = param[1]
                self.__dico[nom_param] = value_param
    
    def __write_file(self):
        with open(self.__url, 'w') as fichier:
            for value in self.__dico.items():
                ligne = value[0] + " : " + value[1] + "\n"
                fichier.write(ligne)
    
    def get_param(self, nom_param):
        try:
            return self.__dico[nom_param]
        except KeyError:
            print "ERROR : Parametre " + nom_param + " non existant"
    
    def modif_param(self, nom_param, new_value):
        self.__dico[nom_param] = str(new_value)
        self.__write_file()
        