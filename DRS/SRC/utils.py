import collections
from scipy import sparse
import numpy as np
from datetime import datetime


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


def read_A_matrix(A_matrix_name):
    A_matrix = sparse.load_npz(A_matrix_name).tolil()
    [art, x_min, x_max] = A_matrix[0, -3:].toarray()[0]
    A_matrix[0, -6:] = 0

    return (A_matrix.tocsr(), art, x_min, x_max)


def header_A_matrix(A_matrix_name):
    A_matrix = sparse.load_npz(A_matrix_name).tolil()
    header = A_matrix[0, -5:-3].toarray()[0]
    return header


def date_A_matrix(A_matrix_name):
    A_matrix = sparse.load_npz(A_matrix_name).tolil()
    date_txt = str(int(A_matrix[0, -6]))
    return datetime.strptime(date_txt, "%Y%m%d%H%M")

    
def rint(nombre):
    return int(round(nombre))


def sign(x):
    return int(np.sign(x))
