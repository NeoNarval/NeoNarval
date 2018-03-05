import numpy as np

# Resultat = [ORDRE, LONGUEUR MIN, LONGUEUR MAX]
# Longueurs d'ondes en mm

def max_min(m):
    num_colonnes = 4612
    f2 = 388
    pixel = 0.0135
    alpha = np.radians(63.477) # angle d'arrivee sur le reseau
    gamma = np.radians(0.6) # inclinaison du reseau selon l'axe horizontal
    G = 79 # nombre de traits par mm pour le reseau

    beta_max = (num_colonnes / 2 *pixel) / f2
    
    def longueur(angle): return (np.sin(alpha) + np.sin(alpha + angle))*np.cos(gamma) / (G * m) * 1e6
    
    return [longueur(-beta_max), longueur(beta_max)]